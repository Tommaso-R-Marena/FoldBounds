#!/usr/bin/env python3
"""Fair matched-budget benchmark harness for foldkernel vs foldcheck / foldcheck2 / Python.

    python run_bench.py --quick                 # DEV subset; ~10 min wall budget
    python run_bench.py --full                  # all suites; ~75 min wall budget
    python run_bench.py --full --only S1,S5     # selected suites
    python run_bench.py --summarize-only         # rebuild SUMMARY.md from JSONL rows

Appends one JSON object per timed configuration to ``results/<suite>.jsonl`` and
regenerates ``results/SUMMARY.md``.  Every number in the summary comes from a
JSONL row produced by an actual run.  Standard library + local bench modules
only (``foldbounds`` is used only via the ``pymethods.py`` child process).
"""
from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import statistics
import sys
import time
import traceback
from fractions import Fraction
from pathlib import Path
from typing import Optional

sys.dont_write_bytecode = True

import benchlib as bl
import certgen

PYMETHODS = bl.BENCH_DIR / "pymethods.py"
HPENUM = bl.BENCH_DIR / "hpenum.py"
FINDINGS = bl.FINDINGS_DIR
RESULTS = bl.RESULTS_DIR

RUN_CAP_S = 300.0
BUDGETS_S = (0.01, 0.1, 1.0, 10.0, 60.0)
QUICK_WALL_S = 9.5 * 60.0
FULL_WALL_S = 74.0 * 60.0

SUITES = ("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8")


# ----------------------------------------------------------------------------
# Global run state
# ----------------------------------------------------------------------------

class State:
    def __init__(self, mode: str, only: Optional[set[str]], wall_budget: float):
        self.mode = mode
        self.only = only
        self.t0 = time.perf_counter()
        self.deadline = self.t0 + wall_budget
        self.critical: list[dict] = []
        self.rows_written = 0
        self.skipped: list[str] = []
        self.machine = {
            "hostname": os.uname().nodename,
            "kernel": os.uname().release,
            "cpus": os.cpu_count(),
            "loadavg_start": bl.loadavg(),
            "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "binaries": {k: {"path": str(v), "sha256": bl.sha256_file(v) if v.exists() else None,
                             "bytes": v.stat().st_size if v.exists() else None}
                        for k, v in bl.EXE.items()},
        }

    def remaining(self) -> float:
        return max(0.0, self.deadline - time.perf_counter())

    def can_start(self, need_s: float = 1.0) -> bool:
        return self.remaining() >= need_s

    def clip_timeout(self, want: float) -> float:
        return max(0.5, min(want, RUN_CAP_S, max(0.5, self.remaining() - 0.25)))


STATE: State


# ----------------------------------------------------------------------------
# Row I/O
# ----------------------------------------------------------------------------

def suite_path(suite: str) -> Path:
    return RESULTS / f"{suite}.jsonl"


def append_row(suite: str, row: dict) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    row = dict(row)
    row.setdefault("suite", suite)
    row.setdefault("mode", STATE.mode)
    row.setdefault("loadavg", bl.loadavg())
    row.setdefault("recorded_utc", datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"))
    with open(suite_path(suite), "a", encoding="ascii") as fh:
        fh.write(json.dumps(row, sort_keys=True) + "\n")
    STATE.rows_written += 1


def load_rows(suite: str) -> list[dict]:
    path = suite_path(suite)
    if not path.exists():
        return []
    return [json.loads(ln) for ln in path.read_text(encoding="ascii").splitlines() if ln.strip()]


def clear_suite(suite: str) -> None:
    path = suite_path(suite)
    if path.exists():
        path.unlink()


# ----------------------------------------------------------------------------
# Critical findings
# ----------------------------------------------------------------------------

def save_finding(tag: str, meta: dict, cert: Optional[bytes] = None) -> Path:
    FINDINGS.mkdir(parents=True, exist_ok=True)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    base = FINDINGS / f"{stamp}_{tag}"
    meta = dict(meta)
    meta["tag"] = tag
    (base.with_suffix(".json")).write_text(json.dumps(meta, indent=1) + "\n", encoding="ascii")
    if cert is not None and len(cert) <= 2_000_000:
        base.with_suffix(".cert").write_bytes(cert)
    STATE.critical.append(meta)
    return base


def check_soundness(row: dict, end: Optional[dict], ref: Optional[dict],
                    cert: Optional[bytes], meta: dict) -> None:
    """Stop on accepted-but-unsound or unexpected reject of a matched cert."""
    if row.get("timed_out"):
        return
    if row.get("accepted") is True and ref is not None and end is not None:
        if not bl.encloses(end, ref):
            path = save_finding("unsound_enclosure", {**meta, "row": row, "end": end, "ref": ref}, cert)
            raise SystemExit(f"CRITICAL: accepted enclosure misses reference; saved {path}")
    if row.get("expected_accept") and row.get("accepted") is False and not row.get("timed_out"):
        # matched generator certificates must be accepted
        path = save_finding("unexpected_reject", {**meta, "row": row}, cert)
        raise SystemExit(f"CRITICAL: matched certificate rejected; saved {path}")


# ----------------------------------------------------------------------------
# Timing helpers
# ----------------------------------------------------------------------------

def _n_repeats(wall_s: float) -> int:
    if wall_s < 0.2:
        return 5
    if wall_s < 2.0:
        return 3
    return 1


def _agg(walls: list[float]) -> dict:
    return {
        "wall_s_median": statistics.median(walls),
        "wall_s_min": min(walls),
        "wall_s_all": walls,
        "n_repeats": len(walls),
    }


def run_checker_once(checker: str, cert: bytes, timeout_s: float,
                     scratch_name: str = "run.cert") -> dict:
    path = bl.write_scratch(scratch_name, cert)
    res = bl.timed_run([str(bl.EXE[checker]), "check", str(path)], timeout_s)
    rep = bl.parse_report(checker, res["stdout"])
    out = {
        "checker": checker,
        "timed_out": res["timed_out"],
        "returncode": res["returncode"],
        "wall_s": res["wall_s"],
        "rss_kb": res["rss_kb"],
        "rss_source": res["rss_source"],
        "cpu_user_s": res["cpu_user_s"],
        "cpu_sys_s": res["cpu_sys_s"],
        "accepted": rep["accepted"],
        "verify_s": rep["verify_s"],
        "results": rep["results"],
        "raw": rep["raw"],
        "error": rep["error"],
        "stderr_tail": res["stderr"][-500:] if res["stderr"] else "",
        "cert_bytes": len(cert),
    }
    return out


def measure_checker(checker: str, cert: bytes, timeout_s: float,
                    scratch_name: str = "run.cert") -> dict:
    """One or more timed runs; returns median/min wall and the last parse."""
    timeout_s = STATE.clip_timeout(timeout_s)
    first = run_checker_once(checker, cert, timeout_s, scratch_name)
    walls = [first["wall_s"]]
    best = first
    if not first["timed_out"] and first["returncode"] == 0:
        nrep = _n_repeats(first["wall_s"])
        for i in range(1, nrep):
            if not STATE.can_start(first["wall_s"] + 0.5):
                break
            r = run_checker_once(checker, cert, timeout_s, f"{scratch_name}.{i}")
            walls.append(r["wall_s"])
            if r["wall_s"] < best["wall_s"]:
                best = r
            if r["timed_out"]:
                break
    out = dict(best)
    out.update(_agg(walls))
    # keep the richest parse (prefer an accepted one)
    for r_wall, r in zip(walls, [first] + ([] if len(walls) == 1 else [])):
        pass
    return out


def measure_python(spec: dict, timeout_s: float) -> dict:
    timeout_s = STATE.clip_timeout(timeout_s)
    cmd = [sys.executable, "-B", str(PYMETHODS), json.dumps(spec)]
    first = bl.timed_run(cmd, timeout_s)
    walls = [first["wall_s"]]
    parsed = _parse_pymethods(first)
    if not first["timed_out"] and parsed.get("accepted"):
        nrep = _n_repeats(first["wall_s"])
        for i in range(1, nrep):
            if not STATE.can_start(first["wall_s"] + 0.5):
                break
            r = bl.timed_run(cmd, timeout_s)
            walls.append(r["wall_s"])
            p = _parse_pymethods(r)
            if p.get("accepted") and (parsed.get("verify_s") is None
                                      or (p.get("verify_s") or 0) < (parsed.get("verify_s") or 1e9)):
                parsed = p
                first = r
            if r["timed_out"]:
                break
    out = {
        "checker": "python",
        "timed_out": first["timed_out"],
        "returncode": first["returncode"],
        "wall_s": first["wall_s"],
        "rss_kb": first["rss_kb"],
        "rss_source": first["rss_source"],
        "cpu_user_s": first["cpu_user_s"],
        "cpu_sys_s": first["cpu_sys_s"],
        "stderr_tail": first["stderr"][-500:] if first["stderr"] else "",
        "cert_bytes": None,
    }
    out.update(parsed)
    out.update(_agg(walls))
    return out


def _parse_pymethods(res: dict) -> dict:
    if res["timed_out"] or not res["stdout"].strip():
        return {"accepted": None, "error": "timeout" if res["timed_out"] else "no output",
                "results": {}, "verify_s": None, "producer_s": None}
    try:
        obj = json.loads(res["stdout"].strip().splitlines()[-1])
    except ValueError as exc:
        return {"accepted": None, "error": f"unparsable: {exc}", "results": {},
                "verify_s": None, "producer_s": None}
    if not obj.get("accepted"):
        return {"accepted": obj.get("accepted"), "error": obj.get("error"), "results": {},
                "verify_s": obj.get("verify_s"), "producer_s": obj.get("producer_s")}
    end = obj["endpoints"]
    return {
        "accepted": True,
        "error": None,
        "verify_s": obj.get("verify_s"),
        "producer_s": obj.get("producer_s"),
        "results": {"_": end},
        "exact_paths": obj.get("exact_paths"),
        "bounded_subtrees": obj.get("bounded_subtrees"),
    }


def measure_hpenum(seq: str, timeout_s: float) -> dict:
    timeout_s = STATE.clip_timeout(timeout_s)
    cmd = [sys.executable, "-B", str(HPENUM), seq]
    first = bl.timed_run(cmd, timeout_s)
    walls = [first["wall_s"]]
    parsed: dict = {"accepted": None, "error": None, "histogram": None, "enum_s": None,
                    "conformations": None}
    if not first["timed_out"] and first["stdout"].strip():
        try:
            obj = json.loads(first["stdout"].strip().splitlines()[-1])
            parsed.update({
                "accepted": True,
                "histogram": obj["histogram"],
                "enum_s": obj["enum_s"],
                "conformations": obj["conformations"],
            })
        except Exception as exc:
            parsed["error"] = repr(exc)
        nrep = _n_repeats(first["wall_s"])
        for i in range(1, nrep):
            if not STATE.can_start(first["wall_s"] + 0.5):
                break
            r = bl.timed_run(cmd, timeout_s)
            walls.append(r["wall_s"])
            if r["timed_out"]:
                break
    out = {
        "checker": "hpenum",
        "timed_out": first["timed_out"],
        "returncode": first["returncode"],
        "wall_s": first["wall_s"],
        "rss_kb": first["rss_kb"],
        "rss_source": first["rss_source"],
        "cpu_user_s": first["cpu_user_s"],
        "cpu_sys_s": first["cpu_sys_s"],
        "stderr_tail": first["stderr"][-500:] if first["stderr"] else "",
        "cert_bytes": None,
        "verify_s": parsed.get("enum_s"),
        "producer_s": 0.0,
    }
    out.update(parsed)
    out.update(_agg(walls))
    return out


def endpoints_of(meas: dict, q: Optional[Fraction] = None, t: Optional[int] = None) -> Optional[dict]:
    results = meas.get("results") or {}
    if not results:
        return None
    if q is not None and t is not None:
        key = bl.qkey(q, t)
        if key in results:
            return results[key]
    # python uses "_" ; single-query checkers have one entry
    return next(iter(results.values()))


def accuracy_of(end: Optional[dict], ref: Optional[dict]) -> dict:
    if end is None:
        return {"p_width": None, "p_width_float": None, "log10_zHi_over_zLo": None,
                "log10_zHi_over_Z": None, "log10_Z_over_zLo": None, "encloses_reference": None,
                "exact_match": None}
    acc = bl.accuracy_fields(end, ref)
    if ref is not None:
        acc["exact_match"] = bl.equals_reference(end, ref)
    else:
        acc["exact_match"] = bl.is_exact(end)
    acc["endpoints"] = {e: end[e] for e in bl.ENDPOINTS}
    return acc


def method_label(checker: str, family: str, bounded: str = "-", kappa=None, **params) -> dict:
    return {
        "method": {
            "checker": checker,
            "family": family,
            "bounded": bounded,
            "kappa": kappa,
            **params,
        }
    }


# ----------------------------------------------------------------------------
# Instance helpers
# ----------------------------------------------------------------------------

def load_instances(names: list[str]) -> list[dict]:
    out = []
    for name in names:
        man = bl.load_manifest(name)
        for inst in man["instances"]:
            inst = dict(inst)
            inst["_manifest"] = name
            out.append(inst)
    return out


def ref_for(inst: dict, q: Fraction, t: int) -> Optional[dict]:
    ref = inst.get("reference")
    if not ref:
        return None
    return ref.get(bl.qkey(q, t))


def mid_target(inst: dict) -> tuple[Fraction, int]:
    thr = inst["thresholds"]
    t = thr.get("mid", thr.get("t_b", 1))
    return Fraction(2), int(t)


def iter_qt(inst: dict, quick: bool) -> list[tuple[Fraction, int]]:
    if inst.get("reference") is None:
        thr = inst["thresholds"]
        ts = [thr.get("t_a", 1), thr.get("t_b", 4), thr.get("t_c", 8)]
        return [(Fraction(2), int(t)) for t in ts]
    keys = list(inst["reference"].keys())
    pairs = []
    for k in keys:
        qs, ts = k.split("|")
        pairs.append((Fraction(qs), int(ts)))
    if quick:
        # one q=2 mid / high / t1 only
        q2 = [p for p in pairs if p[0] == Fraction(2)]
        return q2 or pairs[:3]
    return pairs


# ----------------------------------------------------------------------------
# Suite S1 — exactness
# ----------------------------------------------------------------------------

def suite_s1(quick: bool) -> None:
    """Kernel root-X equals reference; equals foldcheck2 root-X when feasible.

    foldcheck2 root-``X`` enumerates the whole conformation tree inside a
    Mathlib-linked binary, so a dense (instance × q × t) grid is not feasible
    under the suite wall budget.  Protocol:

    * kernel: every reference (q, t) on every n<=18 instance (cheap);
    * foldcheck2: one mid-threshold query at q=2 per instance with n<=14
      (quick: n<=12), hard-capped at 60 s.
    """
    if quick:
        insts = [i for i in load_instances(["dev"]) if i["n"] <= 14]
        v2_n_max = 12
    else:
        insts = [i for i in load_instances(["dev", "heldout"]) if i["n"] <= 18 and i.get("reference")]
        v2_n_max = 14
    for inst in insts:
        if not STATE.can_start(2.0):
            STATE.skipped.append(f"S1:{inst['id']}:deadline")
            break
        seq = inst["sequence"]
        end_k_mid = None
        mk_mid = None
        q_mid, t_mid = mid_target(inst)
        for q, t in iter_qt(inst, quick):
            if not STATE.can_start(1.0):
                break
            ref = ref_for(inst, q, t)
            # kernel full
            t_prod = time.perf_counter()
            cert_k = certgen.kernel_full(seq, [(q, t)])
            prod_s = time.perf_counter() - t_prod
            mk = measure_checker("kernel", cert_k, RUN_CAP_S, "s1_kernel.cert")
            end_k = endpoints_of(mk, q, t)
            if q == q_mid and t == t_mid:
                end_k_mid, mk_mid = end_k, mk
            row = {
                "instance_id": inst["id"], "sequence": seq, "n": inst["n"],
                "q": bl.fstr(q), "t": t, "manifest": inst["_manifest"],
                "producer_s": prod_s, "n_leaves": 1, "expected_accept": True,
                **method_label("kernel", "full", "X"),
                **{k: mk[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                       "rss_kb", "rss_source", "cert_bytes", "error",
                                       "wall_s_median", "wall_s_min", "n_repeats")},
                **accuracy_of(end_k, ref),
            }
            check_soundness(row, end_k, ref, cert_k, {"suite": "S1", "id": inst["id"]})
            append_row("S1", row)

        # one foldcheck2 comparison at the mid target
        if inst["n"] > v2_n_max or not STATE.can_start(5.0):
            continue
        q, t = q_mid, t_mid
        ref = ref_for(inst, q, t)
        t_prod = time.perf_counter()
        cert_v2 = certgen.v2_full(seq, q, t)
        prod_s = time.perf_counter() - t_prod
        mv = measure_checker("v2", cert_v2, 60.0, "s1_v2.cert")
        end_v = endpoints_of(mv, q, t)
        agree = None
        if end_k_mid and end_v and mk_mid and mk_mid.get("accepted") and mv.get("accepted"):
            agree = all(bl.frac(end_k_mid[e]) == bl.frac(end_v[e]) for e in bl.ENDPOINTS)
        row = {
            "instance_id": inst["id"], "sequence": seq, "n": inst["n"],
            "q": bl.fstr(q), "t": t, "manifest": inst["_manifest"],
            "producer_s": prod_s, "n_leaves": 1, "expected_accept": True,
            "agrees_with_kernel": agree,
            **method_label("v2", "full", "X"),
            **{k: mv[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                   "rss_kb", "rss_source", "cert_bytes", "error",
                                   "wall_s_median", "wall_s_min", "n_repeats")},
            **accuracy_of(end_v, ref),
        }
        check_soundness(row, end_v, ref, cert_v2, {"suite": "S1", "id": inst["id"], "v2": True})
        append_row("S1", row)


# ----------------------------------------------------------------------------
# Suite S2 — exhaustive speed
# ----------------------------------------------------------------------------

def suite_s2(quick: bool) -> None:
    mixed = "HPHPPHHPHPPHPHHPPHPH"  # length 20; truncate per n
    if quick:
        ns_kernel = list(range(8, 15))
        ns_alt = list(range(8, 13))
    else:
        ns_kernel = list(range(8, 22))
        ns_alt = list(range(8, 18))

    for kind, mkseq in (("allH", lambda n: "H" * n),
                        ("mixed", lambda n: (mixed * 3)[:n])):
        for n in (ns_kernel if kind == "allH" else ns_alt):
            if not STATE.can_start(2.0):
                STATE.skipped.append(f"S2:{kind}:{n}:deadline")
                return
            seq = mkseq(n)
            q, t = Fraction(2), 1
            # kernel
            cert = certgen.kernel_full(seq, [(q, t)])
            mk = measure_checker("kernel", cert, RUN_CAP_S, "s2_kernel.cert")
            conf, exact = bl.conformations(n)
            us_per = None
            if mk.get("accepted") and not mk["timed_out"] and conf > 0 and mk.get("verify_s"):
                us_per = (mk["verify_s"] * 1e6) / conf
            append_row("S2", {
                "instance_id": f"s2-{kind}-{n}", "sequence": seq, "n": n,
                "q": "2/1", "t": 1, "kind": kind,
                "conformations": conf, "conformations_exact": exact,
                "us_per_conformation": us_per,
                "producer_s": 0.0, "n_leaves": 1, "expected_accept": True,
                **method_label("kernel", "full", "X"),
                **{k: mk[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                       "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                       "n_repeats", "error")},
                **accuracy_of(endpoints_of(mk, q, t), None),
            })
            if mk["timed_out"] or (mk.get("wall_s_median") or 0) > RUN_CAP_S:
                break

            # foldcheck2 — stop at first over-cap (keep n modest; binary is slow)
            if n <= (12 if quick else 14) and STATE.can_start(5.0):
                cert2 = certgen.v2_full(seq, q, t)
                mv = measure_checker("v2", cert2, 60.0, "s2_v2.cert")
                us_per2 = None
                if mv.get("accepted") and not mv["timed_out"] and conf > 0 and mv.get("verify_s"):
                    us_per2 = (mv["verify_s"] * 1e6) / conf
                append_row("S2", {
                    "instance_id": f"s2-{kind}-{n}", "sequence": seq, "n": n,
                    "q": "2/1", "t": 1, "kind": kind,
                    "conformations": conf, "conformations_exact": exact,
                    "us_per_conformation": us_per2,
                    "producer_s": 0.0, "n_leaves": 1, "expected_accept": True,
                    **method_label("v2", "full", "X"),
                    **{k: mv[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                           "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                           "n_repeats", "error")},
                    **accuracy_of(endpoints_of(mv, q, t), None),
                })
                if mv["timed_out"]:
                    # record feasibility limit; continue kernel-only
                    pass

            # pure Python enumerator
            if n <= (12 if quick else 16) and STATE.can_start(5.0):
                mh = measure_hpenum(seq, RUN_CAP_S)
                us_perh = None
                if mh.get("accepted") and not mh["timed_out"] and conf > 0 and mh.get("enum_s"):
                    us_perh = (mh["enum_s"] * 1e6) / conf
                append_row("S2", {
                    "instance_id": f"s2-{kind}-{n}", "sequence": seq, "n": n,
                    "q": "2/1", "t": 1, "kind": kind,
                    "conformations": conf, "conformations_exact": exact,
                    "us_per_conformation": us_perh,
                    "producer_s": 0.0, "n_leaves": None, "expected_accept": True,
                    **method_label("hpenum", "full", "-"),
                    **{k: mh[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                           "rss_kb", "wall_s_median", "wall_s_min",
                                           "n_repeats", "error")},
                    "histogram": mh.get("histogram"),
                })


# ----------------------------------------------------------------------------
# Suites S3 / S4 — per-leaf cost and tightness at identical frontiers
# ----------------------------------------------------------------------------

S3_KINDS = [
    ("kernel", "B", None),
    ("v1", "B", None),
    ("v1", "T", None),
    ("v2", "M0", None),
    ("v2", "M2", None),
    ("v2", "M4", None),
]


def suite_s3_s4(quick: bool) -> None:
    if quick:
        ids = ["dev-fixed-12", "dev-allH-10"]
        depths = (4, 6, 8)
        kinds = S3_KINDS[:5]  # drop M4 in quick
    else:
        ids = ["dev-fixed-12", "dev-rand-n14-p0.7",
               "ho-rand-n13-p0.6", "ho-rand-n13-p0.8"]
        depths = (4, 6, 8, 10)
        kinds = S3_KINDS
    by_id = {i["id"]: i for i in load_instances(["dev", "heldout"])}
    for iid in ids:
        inst = by_id[iid]
        seq = inst["sequence"]
        q, t = mid_target(inst)
        ref = ref_for(inst, q, t)
        for d in depths:
            if d > inst["n"] or not STATE.can_start(3.0):
                continue
            t_prod = time.perf_counter()
            fr = certgen.build_frontier(seq, d, want_words=True)
            prod_s = time.perf_counter() - t_prod
            for checker, bounded, kappa in kinds:
                if not STATE.can_start(2.0):
                    STATE.skipped.append(f"S3:{iid}:d{d}:{checker}{bounded}")
                    return
                if checker == "kernel":
                    cert = certgen.kernel_text(seq, [(q, t)], certgen.kernel_tree(fr, kappa))
                else:
                    ver = 1 if checker == "v1" else 2
                    cert = certgen.flat_text(fr, q, t, version=ver, bounded=bounded, kappa=kappa)
                m = measure_checker(checker, cert, RUN_CAP_S, f"s3_{checker}.cert")
                end = endpoints_of(m, q, t)
                n_leaves = fr.n_depth_leaves + fr.n_exact
                us_per_leaf = None
                if m.get("verify_s") and n_leaves > 0 and m.get("accepted"):
                    us_per_leaf = (m["verify_s"] * 1e6) / n_leaves
                base = {
                    "instance_id": iid, "sequence": seq, "n": inst["n"],
                    "q": bl.fstr(q), "t": t, "manifest": inst["_manifest"],
                    "producer_s": prod_s, "n_leaves": n_leaves,
                    "n_depth_leaves": fr.n_depth_leaves, "n_exact": fr.n_exact,
                    "us_per_leaf": us_per_leaf, "expected_accept": True,
                    "depth": d,
                    **method_label(checker, f"uniform({d})", bounded, kappa, depth=d),
                    **{k: m[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                           "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                           "n_repeats", "error")},
                    **accuracy_of(end, ref),
                }
                check_soundness(base, end, ref, cert,
                                {"suite": "S3", "id": iid, "d": d, "checker": checker, "bounded": bounded})
                append_row("S3", base)
                append_row("S4", dict(base))  # same rows feed tightness tables
                if m["timed_out"]:
                    break


# ----------------------------------------------------------------------------
# Suite S5 — accuracy per unit of checker time (main fairness test)
# ----------------------------------------------------------------------------

def _kernel_s5_configs(n: int, quick: bool) -> list[dict]:
    ds = list(range(4, n + 1, 2))
    if quick:
        ds = [d for d in ds if d <= min(n, 10)]
        kappas: list[Optional[int]] = [None, 1, 2]
    else:
        kappas = [None, 1, 2, 3]
    out = [{"family": "full", "depth": None, "kappa": None, "bounded": "X"}]
    for d in ds:
        for kap in kappas:
            out.append({"family": "hybrid", "depth": d, "kappa": kap, "bounded": "B"})
    # adaptive budgets via producer if present
    if (bl.PRODUCER_PATH.exists()) and not quick:
        for bud in (1000, 10000, 100000, 1_000_000):
            out.append({"family": "adaptive", "depth": None, "kappa": None,
                        "bounded": "mixed", "budget": bud})
    elif bl.PRODUCER_PATH.exists() and quick:
        out.append({"family": "adaptive", "depth": None, "kappa": None,
                    "bounded": "mixed", "budget": 10000})
    return out


def _v2_s5_configs(quick: bool) -> list[dict]:
    ds = (4, 6, 8) if quick else (4, 6, 8, 10)
    out = [{"family": "full", "depth": None, "kappa": None, "bounded": "X"}]
    for d in ds:
        out.append({"family": "hybrid", "depth": d, "kappa": None, "bounded": "M2"})
        out.append({"family": "hybrid", "depth": d, "kappa": 1, "bounded": "M2"})
        if not quick:
            out.append({"family": "hybrid", "depth": d, "kappa": None, "bounded": "M0"})
    return out


def _v1_s5_configs(quick: bool) -> list[dict]:
    ds = (4, 6, 8) if quick else (4, 6, 8, 10)
    return [{"family": "uniform", "depth": d, "kappa": None, "bounded": "T"} for d in ds] + \
           [{"family": "uniform", "depth": d, "kappa": None, "bounded": "B"} for d in ds]


def _py_s5_configs(n: int, quick: bool) -> list[dict]:
    out = []
    for d in ((4, 8) if quick else (4, 8, min(12, n))):
        if d <= n:
            out.append({"method": "basic", "params": {"depth": d}})
    LA = dict(use_matching=True, use_transfer=True)
    for k in ((0, 2, 4) if quick else (0, 2, 4, 6)):
        out.append({"method": "lookahead", "params": {"depth": 2, "lookahead": k}})
        out.append({"method": "lookahead", "params": {"depth": 2, "lookahead": k, "use_matching": True}})
    out.append({"method": "lookahead", "params": {"depth": 2, "lookahead": 0, **LA}})
    out.append({"method": "lookahead",
                "params": {"depth": 2, "lookahead": 0, "use_monotone": True,
                           "transfer_memory": 4, **LA}})
    if not quick:
        out.append({"method": "lookahead",
                    "params": {"depth": 2, "lookahead": 0, "use_monotone": True,
                               "use_turn_family": True, "transfer_memory": 4,
                               "terminal_box": True, **LA}})
        for pol in ("breadth", "mass", "probability"):
            for s in (4, 8, 16, 32):
                out.append({"method": "adaptive",
                            "params": {"depth": 2, "steps": s, "policy": pol}})
    else:
        out.append({"method": "adaptive",
                    "params": {"depth": 2, "steps": 8, "policy": "mass"}})
    return out


def _make_kernel_cert(seq: str, q: Fraction, t: int, cfg: dict) -> tuple[bytes, float, int]:
    t0 = time.perf_counter()
    if cfg["family"] == "full":
        cert = certgen.kernel_full(seq, [(q, t)])
        n_leaves = 1
    elif cfg["family"] == "adaptive":
        sys.path.insert(0, str(bl.REPO / "experiments" / "kernel"))
        import producer  # type: ignore
        tree = producer.adaptive(seq, q, t, int(cfg["budget"]))
        cert = certgen.kernel_text(seq, [(q, t)], tree.encode("ascii"))
        n_leaves = tree.count("B") + tree.count("X") + tree.count("E")
    else:
        fr = certgen.build_frontier(seq, int(cfg["depth"]), want_words=False)
        cert = certgen.kernel_text(seq, [(q, t)],
                                   certgen.kernel_tree(fr, cfg["kappa"]))
        n_leaves = fr.n_depth_leaves + fr.n_exact
    return cert, time.perf_counter() - t0, n_leaves


def _make_flat_cert(seq: str, q: Fraction, t: int, checker: str, cfg: dict) -> tuple[bytes, float, int]:
    t0 = time.perf_counter()
    if cfg["family"] == "full":
        cert = certgen.v2_full(seq, q, t)
        return cert, time.perf_counter() - t0, 1
    fr = certgen.build_frontier(seq, int(cfg["depth"]), want_words=True)
    ver = 1 if checker == "v1" else 2
    cert = certgen.flat_text(fr, q, t, version=ver, bounded=cfg["bounded"], kappa=cfg["kappa"])
    return cert, time.perf_counter() - t0, fr.n_depth_leaves + fr.n_exact


def suite_s5(quick: bool) -> None:
    by_id = {i["id"]: i for i in load_instances(["dev", "heldout"])}
    if quick:
        ids = ["dev-allH-10", "dev-fixed-12"]
    else:
        ids = [
            "dev-fixed-12", "dev-rand-n14-p0.7", "dev-rand-n16-p0.7",
            "ho-rand-n13-p0.6", "ho-rand-n13-p0.8", "ho-rand-n15-p0.6",
            "ho-rand-n15-p0.8",
        ]
        # n=17 only if ≥12 minutes remain when we reach it
        if STATE.remaining() >= 12 * 60:
            ids.append("ho-rand-n17-p0.6")
    for iid in ids:
        if iid not in by_id:
            continue
        inst = by_id[iid]
        if not STATE.can_start(10.0):
            STATE.skipped.append(f"S5:{iid}:deadline")
            break
        seq = inst["sequence"]
        q, t = mid_target(inst)
        ref = ref_for(inst, q, t)

        # kernel grid
        for cfg in _kernel_s5_configs(inst["n"], quick):
            if not STATE.can_start(2.0):
                break
            try:
                cert, prod_s, n_leaves = _make_kernel_cert(seq, q, t, cfg)
            except Exception as exc:
                append_row("S5", {
                    "instance_id": iid, "sequence": seq, "n": inst["n"],
                    "q": bl.fstr(q), "t": t, "error": f"producer: {exc}",
                    "timed_out": False, "accepted": None,
                    **method_label("kernel", cfg["family"], cfg.get("bounded", "B"),
                                   cfg.get("kappa"), **{k: cfg[k] for k in cfg
                                                        if k not in ("family", "bounded", "kappa")}),
                })
                continue
            m = measure_checker("kernel", cert, RUN_CAP_S, "s5_kernel.cert")
            end = endpoints_of(m, q, t)
            row = {
                "instance_id": iid, "sequence": seq, "n": inst["n"],
                "q": bl.fstr(q), "t": t, "manifest": inst["_manifest"],
                "producer_s": prod_s, "n_leaves": n_leaves, "expected_accept": True,
                "headline_split": "heldout" if inst["_manifest"] == "heldout" else "dev",
                **method_label("kernel", cfg["family"], cfg.get("bounded", "B"),
                               cfg.get("kappa"),
                               **{k: cfg[k] for k in ("depth", "budget") if k in cfg}),
                **{k: m[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                       "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                       "n_repeats", "error")},
                **accuracy_of(end, ref),
            }
            check_soundness(row, end, ref, cert, {"suite": "S5", "id": iid, "cfg": cfg})
            append_row("S5", row)
            if m["timed_out"] and cfg["family"] == "full":
                break

        # v2 grid (shorter hard cap on larger n so one slow M2 cannot eat the suite)
        alt_cap = 60.0 if inst["n"] >= 15 else RUN_CAP_S
        for cfg in _v2_s5_configs(quick):
            if not STATE.can_start(5.0):
                break
            if inst["n"] >= 15 and cfg.get("bounded") == "M0" and cfg.get("depth", 0) and cfg["depth"] >= 8:
                continue
            cert, prod_s, n_leaves = _make_flat_cert(seq, q, t, "v2", cfg)
            m = measure_checker("v2", cert, alt_cap, "s5_v2.cert")
            end = endpoints_of(m, q, t)
            row = {
                "instance_id": iid, "sequence": seq, "n": inst["n"],
                "q": bl.fstr(q), "t": t, "manifest": inst["_manifest"],
                "producer_s": prod_s, "n_leaves": n_leaves, "expected_accept": True,
                "headline_split": "heldout" if inst["_manifest"] == "heldout" else "dev",
                **method_label("v2", cfg["family"], cfg["bounded"], cfg.get("kappa"),
                               **{k: cfg[k] for k in ("depth",) if k in cfg and cfg[k] is not None}),
                **{k: m[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                       "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                       "n_repeats", "error")},
                **accuracy_of(end, ref),
            }
            check_soundness(row, end, ref, cert, {"suite": "S5", "id": iid, "cfg": cfg})
            append_row("S5", row)
            if m["timed_out"]:
                break

        # v1 grid
        for cfg in _v1_s5_configs(quick):
            if not STATE.can_start(5.0):
                break
            cert, prod_s, n_leaves = _make_flat_cert(seq, q, t, "v1", cfg)
            m = measure_checker("v1", cert, alt_cap, "s5_v1.cert")
            end = endpoints_of(m, q, t)
            row = {
                "instance_id": iid, "sequence": seq, "n": inst["n"],
                "q": bl.fstr(q), "t": t, "manifest": inst["_manifest"],
                "producer_s": prod_s, "n_leaves": n_leaves, "expected_accept": True,
                "headline_split": "heldout" if inst["_manifest"] == "heldout" else "dev",
                **method_label("v1", cfg["family"], cfg["bounded"], None, depth=cfg["depth"]),
                **{k: m[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                       "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                       "n_repeats", "error")},
                **accuracy_of(end, ref),
            }
            check_soundness(row, end, ref, cert, {"suite": "S5", "id": iid, "cfg": cfg})
            append_row("S5", row)
            if m["timed_out"]:
                break

        # Python methods
        for cfg in _py_s5_configs(inst["n"], quick):
            if not STATE.can_start(3.0):
                break
            spec = {"sequence": seq, "q": bl.fstr(q), "t": t,
                    "method": cfg["method"], "params": cfg["params"]}
            m = measure_python(spec, RUN_CAP_S)
            end = endpoints_of(m)
            row = {
                "instance_id": iid, "sequence": seq, "n": inst["n"],
                "q": bl.fstr(q), "t": t, "manifest": inst["_manifest"],
                "producer_s": m.get("producer_s"), "n_leaves": m.get("bounded_subtrees"),
                "expected_accept": True,
                "headline_split": "heldout" if inst["_manifest"] == "heldout" else "dev",
                "unverified": True,
                **method_label("python", cfg["method"], "-", None, **cfg["params"]),
                **{k: m[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                       "rss_kb", "wall_s_median", "wall_s_min",
                                       "n_repeats", "error")},
                **accuracy_of(end, ref),
            }
            # Python is unverified: enclosure failure is a finding but not a Lean bug
            if row.get("accepted") and ref is not None and end is not None and not bl.encloses(end, ref):
                save_finding("python_misses_reference", {"suite": "S5", "id": iid, "row": row})
            append_row("S5", row)
            if m["timed_out"]:
                break


# ----------------------------------------------------------------------------
# Suite S6 — non-enumeration / scaling
# ----------------------------------------------------------------------------

def suite_s6(quick: bool) -> None:
    by_id = {i["id"]: i for i in load_instances(["heldout"])}
    seqs = []
    if "ho-long-30" in by_id:
        seqs.append(("ho-long-30", by_id["ho-long-30"]["sequence"], 1))
    if not quick and "ho-long-48" in by_id:
        seqs.append(("ho-long-48", by_id["ho-long-48"]["sequence"], 1))
    for n in (range(18, 20) if quick else range(18, 23)):
        seqs.append((f"allH-{n}", "H" * n, n // 4))
    depths = (2, 4, 6, 8) if quick else (2, 4, 6, 8, 10, 12)
    for iid, seq, t in seqs:
        n = len(seq)
        q = Fraction(2)
        conf, exact = bl.conformations(n)
        for d in depths:
            if d > n or not STATE.can_start(2.0):
                continue
            t0 = time.perf_counter()
            fr = certgen.build_frontier(seq, d, want_words=True)
            prod_s = time.perf_counter() - t0
            cert_k = certgen.kernel_text(seq, [(q, t)], certgen.kernel_tree(fr, None))
            mk = measure_checker("kernel", cert_k, RUN_CAP_S, "s6_kernel.cert")
            steps = None
            if mk.get("accepted") and mk.get("raw"):
                splits = mk["raw"].get("split_nodes") or 0
                leaves = (mk["raw"].get("bound_leaves") or 0) + (mk["raw"].get("enum_leaves") or 0) \
                         + (mk["raw"].get("exact_leaves") or 0)
                steps = 4 * splits + leaves
            append_row("S6", {
                "instance_id": iid, "sequence": seq, "n": n, "q": "2/1", "t": t,
                "conformations": conf, "conformations_exact": exact,
                "checker_steps": steps,
                "coverage_ratio": (None if steps in (None, 0) else conf / steps),
                "producer_s": prod_s, "n_leaves": fr.n_depth_leaves + fr.n_exact,
                "depth": d, "expected_accept": True,
                **method_label("kernel", f"uniform({d})", "B", None, depth=d),
                **{k: mk[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                       "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                       "n_repeats", "error", "raw")},
                **accuracy_of(endpoints_of(mk, q, t), None),
            })
            # matched B for v1/v2 only on moderate n/depth — v2 M0 on n=30 already
            # burned a full 300 s cap in probing; keep a short feasibility probe.
            if d <= (4 if quick else 6) and n <= 22 and STATE.can_start(5.0):
                for checker, bounded in (("v1", "B"), ("v2", "M0")):
                    if not STATE.can_start(3.0):
                        break
                    ver = 1 if checker == "v1" else 2
                    cert = certgen.flat_text(fr, q, t, version=ver, bounded=bounded, kappa=None)
                    alt_cap = 30.0 if n >= 18 else RUN_CAP_S
                    m = measure_checker(checker, cert, alt_cap, f"s6_{checker}.cert")
                    append_row("S6", {
                        "instance_id": iid, "sequence": seq, "n": n, "q": "2/1", "t": t,
                        "conformations": conf, "conformations_exact": exact,
                        "producer_s": prod_s, "n_leaves": fr.n_depth_leaves + fr.n_exact,
                        "depth": d, "expected_accept": True,
                        **method_label(checker, f"uniform({d})", bounded, None, depth=d),
                        **{k: m[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                               "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                               "n_repeats", "error")},
                        **accuracy_of(endpoints_of(m, q, t), None),
                    })
                    if m["timed_out"]:
                        break


# ----------------------------------------------------------------------------
# Suite S7 — one tree, many queries
# ----------------------------------------------------------------------------

def suite_s7(quick: bool) -> None:
    by_id = {i["id"]: i for i in load_instances(["dev"])}
    iid = "dev-fixed-12"
    inst = by_id[iid]
    seq = inst["sequence"]
    # build one hybrid frontier
    d = 8 if not quick else 6
    fr = certgen.build_frontier(seq, d, want_words=True)
    tree = certgen.kernel_tree(fr, 1)
    # query sets
    all_q = [Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3)]
    all_t = sorted(set(inst["thresholds"].values()))
    universe = [(q, t) for q in all_q for t in all_t]
    counts = (1, 4) if quick else (1, 4, 16, min(64, len(universe)))
    # pad universe if needed
    while len(universe) < 64:
        universe.append((Fraction(2), 1 + (len(universe) % 5)))
    for nq in counts:
        if not STATE.can_start(2.0):
            break
        queries = universe[:nq]
        cert = certgen.kernel_text(seq, queries, tree)
        mk = measure_checker("kernel", cert, RUN_CAP_S, "s7_kernel.cert")
        append_row("S7", {
            "instance_id": iid, "sequence": seq, "n": inst["n"],
            "n_queries": nq, "depth": d, "producer_s": 0.0,
            "n_leaves": fr.n_depth_leaves + fr.n_exact, "expected_accept": True,
            **method_label("kernel", f"hybrid({d},1)", "B", 1, n_queries=nq, depth=d),
            **{k: mk[k] for k in ("timed_out", "returncode", "accepted", "verify_s",
                                   "rss_kb", "cert_bytes", "wall_s_median", "wall_s_min",
                                   "n_repeats", "error")},
            "per_query_verify_s": (None if not mk.get("verify_s") else mk["verify_s"] / nq),
        })
        # foldcheck2 once per query on matched cert
        if nq <= (4 if quick else 16) and STATE.can_start(nq * 2.0):
            walls = []
            ok = True
            for q, t in queries:
                if not STATE.can_start(2.0):
                    ok = False
                    break
                cert2 = certgen.flat_text(fr, q, t, version=2, bounded="M2", kappa=1)
                mv = measure_checker("v2", cert2, RUN_CAP_S, "s7_v2.cert")
                walls.append(mv.get("wall_s_median") or mv["wall_s"])
                if mv["timed_out"] or not mv.get("accepted"):
                    ok = False
                    break
            append_row("S7", {
                "instance_id": iid, "sequence": seq, "n": inst["n"],
                "n_queries": nq, "depth": d, "producer_s": 0.0,
                "n_leaves": fr.n_depth_leaves + fr.n_exact,
                "expected_accept": True, "all_accepted": ok,
                "wall_s_sum": sum(walls) if walls else None,
                "wall_s_median": (statistics.median(walls) if walls else None),
                **method_label("v2", f"hybrid({d},1)", "M2", 1, n_queries=nq, depth=d,
                               aggregate="sum_of_per_query_runs"),
                "per_query_walls": walls,
            })


# ----------------------------------------------------------------------------
# Suite S8 — footprint
# ----------------------------------------------------------------------------

def suite_s8(_quick: bool) -> None:
    row = {
        "instance_id": "binaries",
        "binaries": STATE.machine["binaries"],
        "build_notes_coordinator": {
            "foldcheck_plus_foldcheck2_wall_s": 8 * 60 + 4,
            "foldcheck_plus_foldcheck2_cpu_min": 25,
            "foldkernel_wall_s": 1.4,
            "source": "coordinator-provided; not re-measured by this harness",
        },
        **method_label("meta", "footprint", "-"),
    }
    # size ratios
    bk = STATE.machine["binaries"]["kernel"]["bytes"]
    b1 = STATE.machine["binaries"]["v1"]["bytes"]
    b2 = STATE.machine["binaries"]["v2"]["bytes"]
    if bk and b1 and b2:
        row["size_ratio_v1_over_kernel"] = b1 / bk
        row["size_ratio_v2_over_kernel"] = b2 / bk
    append_row("S8", row)


# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------

def _method_key(row: dict) -> str:
    m = row.get("method") or {}
    parts = [m.get("checker", "?"), m.get("family", "?"), str(m.get("bounded", "-"))]
    if m.get("depth") is not None:
        parts.append(f"d{m['depth']}")
    if m.get("kappa") is not None:
        parts.append(f"k{m['kappa']}")
    if m.get("budget") is not None:
        parts.append(f"b{m['budget']}")
    if m.get("method") and m["checker"] == "python":
        # pymethods puts name in family
        pass
    for k in ("lookahead", "steps", "policy", "transfer_memory", "witness_nodes"):
        if k in m and m[k] is not None:
            parts.append(f"{k}={m[k]}")
    return "/".join(parts)


def _checker_family(row: dict) -> str:
    m = row.get("method") or {}
    c = m.get("checker", "?")
    if c == "python":
        return f"python:{m.get('family')}"
    if c == "kernel":
        return "kernel"
    if c == "v1":
        return f"v1:{m.get('bounded')}"
    if c == "v2":
        return f"v2:{m.get('bounded')}"
    return c


def summarize() -> str:
    lines: list[str] = []
    lines.append("# FoldBounds certified frontier — matched-budget benchmark SUMMARY")
    lines.append("")
    lines.append(f"Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("## Machine / binaries")
    lines.append("")
    lines.append(f"- Host: `{STATE.machine['hostname']}` · Linux `{STATE.machine['kernel']}` · "
                 f"{STATE.machine['cpus']} CPUs")
    lines.append(f"- Load average at harness start: {STATE.machine['loadavg_start']}")
    lines.append(f"- Mode: `{STATE.mode}` · rows written this run: {STATE.rows_written}")
    for k, info in STATE.machine["binaries"].items():
        if info["bytes"] is None:
            lines.append(f"- `{k}`: MISSING at `{info['path']}`")
        else:
            lines.append(f"- `{k}`: {info['bytes']:,} bytes · sha256 `{info['sha256'][:16]}…`")
    lines.append("")
    lines.append("## Commands")
    lines.append("")
    lines.append("```bash")
    lines.append('export PATH="$HOME/.elan/bin:$PATH"')
    lines.append("cd /workspace && . .venv/bin/activate")
    lines.append("lake build foldkernel   # foldcheck / foldcheck2 already built")
    lines.append("python experiments/kernel/bench/run_bench.py --quick")
    lines.append("python experiments/kernel/bench/run_bench.py --full")
    lines.append("```")
    lines.append("")
    lines.append("## Findings")
    lines.append("")
    findings_kernel: list[str] = []
    findings_tie: list[str] = []
    findings_alt: list[str] = []
    findings_inconclusive: list[str] = []

    # ---- S1 ----
    s1 = load_rows("S1")
    lines.append("## S1 — Exactness (just as much accuracy)")
    lines.append("")
    lines.append("Protocol: root-`X` certificates; kernel must match the frozen exact "
                 "reference on every endpoint; foldcheck2 root-`X` compared once per "
                 "instance (q=2, mid threshold, n≤14) with a 60s hard cap. "
                 "Repeats: median/min wall.")
    lines.append("")
    k_rows = [r for r in s1 if (r.get("method") or {}).get("checker") == "kernel"]
    v_rows = [r for r in s1 if (r.get("method") or {}).get("checker") == "v2"]
    k_exact = sum(1 for r in k_rows if r.get("exact_match") is True)
    k_encl = sum(1 for r in k_rows if r.get("encloses_reference") is True)
    k_fail = sum(1 for r in k_rows if r.get("accepted") is not True or r.get("timed_out"))
    agree = sum(1 for r in v_rows if r.get("agrees_with_kernel") is True)
    disagree = sum(1 for r in v_rows if r.get("agrees_with_kernel") is False)
    v_to = sum(1 for r in v_rows if r.get("timed_out"))
    lines.append(f"| metric | value |")
    lines.append(f"|---|---|")
    lines.append(f"| kernel rows | {len(k_rows)} |")
    lines.append(f"| kernel exact_match | {k_exact} |")
    lines.append(f"| kernel encloses_reference | {k_encl} |")
    lines.append(f"| kernel failures/timeouts | {k_fail} |")
    lines.append(f"| foldcheck2 agrees_with_kernel | {agree} |")
    lines.append(f"| foldcheck2 disagrees | {disagree} |")
    lines.append(f"| foldcheck2 timeouts | {v_to} |")
    lines.append("")
    if k_fail == 0 and k_exact == len(k_rows) and len(k_rows) > 0:
        findings_tie.append(
            f"S1: kernel root-X exact on {k_exact}/{len(k_rows)} reference queries "
            f"(disagreements with reference = 0).")
    if disagree == 0 and agree > 0:
        findings_tie.append(f"S1: foldcheck2 root-X agreed with kernel on {agree} finished comparisons "
                            f"(disagreements = 0; timeouts = {v_to}).")
    elif disagree > 0:
        findings_alt.append(f"S1: foldcheck2 disagreed with kernel on {disagree} comparisons.")

    # ---- S2 ----
    s2 = load_rows("S2")
    lines.append("## S2 — Exhaustive speed")
    lines.append("")
    lines.append("Protocol: root-`X` (kernel, foldcheck2) and pure-Python `hpenum` on all-H and "
                 "one mixed family; stop each alternative at first over-cap run.")
    lines.append("")
    lines.append("| kind | n | method | wall_s_median | verify_s | µs/conf | timed_out |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in sorted(s2, key=lambda x: (x.get("kind") or "", x.get("n") or 0,
                                       (x.get("method") or {}).get("checker") or "")):
        m = r.get("method") or {}
        lines.append(
            f"| {r.get('kind')} | {r.get('n')} | {m.get('checker')} | "
            f"{_fmt(r.get('wall_s_median'))} | {_fmt(r.get('verify_s'))} | "
            f"{_fmt(r.get('us_per_conformation'))} | {r.get('timed_out')} |")
    lines.append("")
    # aggregate speedups at shared n (report min/median/max, not every row)
    for other in ("v2", "hpenum"):
        ratios = []
        detail = []
        for kind in ("allH", "mixed"):
            for n in sorted({r["n"] for r in s2 if r.get("kind") == kind}):
                sub = [r for r in s2 if r.get("kind") == kind and r["n"] == n and not r.get("timed_out")]
                by_c = {(r.get("method") or {}).get("checker"): r for r in sub}
                if "kernel" in by_c and other in by_c and by_c["kernel"].get("wall_s_median"):
                    kw = by_c["kernel"]["wall_s_median"]
                    if kw > 0:
                        sp = by_c[other]["wall_s_median"] / kw
                        ratios.append(sp)
                        detail.append((kind, n, sp, by_c[other]["wall_s_median"], kw))
        if ratios:
            kind_n, n, sp, ow, kw = max(detail, key=lambda x: x[2])
            label = "foldcheck2" if other == "v2" else "hpenum"
            findings_kernel.append(
                f"S2: {label}/kernel wall speedup over {len(ratios)} matched (kind,n) pairs: "
                f"min {min(ratios):.1f}× · median {statistics.median(ratios):.1f}× · "
                f"max {max(ratios):.1f}× (max at {kind_n}/n={n}: {ow:.4g}s / {kw:.4g}s).")

    # ---- S3 ----
    s3 = load_rows("S3")
    lines.append("## S3 — Per-leaf cost at identical frontiers")
    lines.append("")
    lines.append("Protocol: matched `uniform(d)` frontiers; µs/leaf = verify_s / n_leaves.")
    lines.append("")
    lines.append("| instance | d | method | leaves | verify_s | µs/leaf | wall_s_median |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in sorted(s3, key=lambda x: (x.get("instance_id") or "", x.get("depth") or 0,
                                       _method_key(x))):
        if r.get("timed_out"):
            continue
        m = r.get("method") or {}
        lines.append(
            f"| {r.get('instance_id')} | {r.get('depth')} | "
            f"{m.get('checker')}:{m.get('bounded')} | {r.get('n_leaves')} | "
            f"{_fmt(r.get('verify_s'))} | {_fmt(r.get('us_per_leaf'))} | "
            f"{_fmt(r.get('wall_s_median'))} |")
    lines.append("")
    # compare kernel B vs others at same instance/depth
    from collections import defaultdict
    groups = defaultdict(list)
    for r in s3:
        if r.get("accepted") and not r.get("timed_out") and r.get("us_per_leaf"):
            groups[(r["instance_id"], r["depth"])].append(r)
    s3_slower, s3_faster = [], []
    for key, rows in groups.items():
        kern = [r for r in rows if (r.get("method") or {}).get("checker") == "kernel"]
        if not kern:
            continue
        ku = kern[0]["us_per_leaf"]
        for r in rows:
            if (r.get("method") or {}).get("checker") == "kernel":
                continue
            ratio = r["us_per_leaf"] / ku if ku else None
            label = (f"{(r.get('method') or {}).get('checker')}:"
                     f"{(r.get('method') or {}).get('bounded')}")
            if ratio and ratio > 1.05:
                s3_slower.append((ratio, key, label, r["us_per_leaf"], ku))
            elif ratio and ratio < 0.95:
                s3_faster.append((ratio, key, label, r["us_per_leaf"], ku))
    if s3_slower:
        s3_slower.sort(reverse=True)
        r0 = s3_slower[0]
        findings_kernel.append(
            f"S3: in {len(s3_slower)} matched (instance,depth,alt) comparisons the alternative "
            f"was slower per leaf than kernel B; largest {r0[0]:.1f}× "
            f"({r0[2]} at {r0[1][0]}/d={r0[1][1]}: {r0[3]:.3g} vs {r0[4]:.3g} µs).")
    if s3_faster:
        findings_alt.append(
            f"S3: {len(s3_faster)} matched comparisons where an alternative was cheaper per leaf "
            f"than kernel B (see S3 table).")

    # ---- S4 ----
    s4 = load_rows("S4")
    lines.append("## S4 — Tightness at identical frontiers")
    lines.append("")
    lines.append("Protocol: same frontiers as S3; compare log10(zHi/Z), log10(Z/zLo), p_width.")
    lines.append("")
    lines.append("| instance | d | method | p_width_float | log10(zHi/Z) | log10(Z/zLo) | wall_s |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in sorted(s4, key=lambda x: (x.get("instance_id") or "", x.get("depth") or 0,
                                       _method_key(x))):
        if not r.get("accepted") or r.get("timed_out"):
            continue
        m = r.get("method") or {}
        lines.append(
            f"| {r.get('instance_id')} | {r.get('depth')} | "
            f"{m.get('checker')}:{m.get('bounded')} | {_fmt(r.get('p_width_float'))} | "
            f"{_fmt(r.get('log10_zHi_over_Z'))} | {_fmt(r.get('log10_Z_over_zLo'))} | "
            f"{_fmt(r.get('wall_s_median'))} |")
    lines.append("")
    s4_alt_tight, s4_kern_tight = 0, 0
    s4_examples = []
    for key, rows in groups.items():
        scored = [(r.get("p_width_float"), r) for r in rows
                  if r.get("p_width_float") is not None and r.get("accepted")]
        if len(scored) < 2:
            continue
        scored.sort(key=lambda x: x[0])
        best_w, best_r = scored[0]
        bm = best_r.get("method") or {}
        kern = [r for r in rows if (r.get("method") or {}).get("checker") == "kernel"
                and r.get("p_width_float") is not None]
        if not kern:
            continue
        kw = kern[0]["p_width_float"]
        if bm.get("checker") != "kernel" and best_w + 1e-15 < kw:
            s4_alt_tight += 1
            s4_examples.append(
                f"{key[0]}/d={key[1]}: {bm.get('checker')}:{bm.get('bounded')} "
                f"= {best_w:.6g} vs kernel B = {kw:.6g} "
                f"(wall {_fmt(best_r.get('wall_s_median'))}s vs {_fmt(kern[0].get('wall_s_median'))}s)")
        elif bm.get("checker") == "kernel" or abs(best_w - kw) <= 1e-15 * max(1.0, kw):
            s4_kern_tight += 1
    if s4_kern_tight:
        findings_kernel.append(
            f"S4: kernel B is tightest or tied on p_width in {s4_kern_tight} matched frontiers.")
    if s4_alt_tight:
        findings_alt.append(
            f"S4: an alternative leaf kind is strictly tighter on p_width in {s4_alt_tight} "
            f"matched frontiers (transfer/memory DP can beat kernel B per leaf). "
            f"Example: {s4_examples[0]}.")

    # ---- S5 ----
    s5 = load_rows("S5")
    lines.append("## S5 — Accuracy per unit of checker time (headline fairness test)")
    lines.append("")
    lines.append("Protocol: fixed configuration grids per method; for each wall-clock budget in "
                 f"{list(BUDGETS_S)}s, pick the configuration with the smallest certified "
                 "`p_width` whose **checker** `wall_s_median` is ≤ budget (producer time reported "
                 "separately, not charged to the budget). Held-out instances are the headline.")
    lines.append("")
    for split in ("heldout", "dev"):
        sub = [r for r in s5 if r.get("headline_split") == split
               and r.get("accepted") and not r.get("timed_out")
               and r.get("p_width_float") is not None
               and r.get("wall_s_median") is not None]
        if not sub:
            lines.append(f"### {split}: no completed rows")
            lines.append("")
            continue
        lines.append(f"### {split}")
        lines.append("")
        # per instance, per budget
        inst_ids = sorted({r["instance_id"] for r in sub})
        families = sorted({_checker_family(r) for r in sub})
        wins = {f: 0 for f in families}
        lines.append("| instance | budget_s | best_family | best_width | best_wall_s | "
                     "kernel_width | kernel_wall_s |")
        lines.append("|---|---|---|---|---|---|---|")
        for iid in inst_ids:
            irows = [r for r in sub if r["instance_id"] == iid]
            for B in BUDGETS_S:
                # best overall
                cand = [r for r in irows if r["wall_s_median"] <= B]
                if not cand:
                    lines.append(f"| {iid} | {B} | — | — | — | — | — |")
                    continue
                best = min(cand, key=lambda r: (r["p_width_float"], r["wall_s_median"]))
                bf = _checker_family(best)
                wins[bf] = wins.get(bf, 0) + 1
                kcand = [r for r in cand if (r.get("method") or {}).get("checker") == "kernel"]
                if kcand:
                    kb = min(kcand, key=lambda r: (r["p_width_float"], r["wall_s_median"]))
                    kw, kwall = kb["p_width_float"], kb["wall_s_median"]
                else:
                    kw = kwall = None
                lines.append(
                    f"| {iid} | {B} | {bf} | {best['p_width_float']:.6g} | "
                    f"{best['wall_s_median']:.4g} | {_fmt(kw)} | {_fmt(kwall)} |")
                if split == "heldout" and kw is not None and not bf.startswith("kernel"):
                    if best["p_width_float"] + 1e-15 < kw:
                        findings_alt.append(
                            f"S5/heldout/{iid}@budget={B}s: {bf} width "
                            f"{best['p_width_float']:.6g} beats kernel {kw:.6g}.")
                    elif abs(best["p_width_float"] - kw) <= 1e-15 * max(1.0, kw):
                        findings_tie.append(
                            f"S5/heldout/{iid}@budget={B}s: tie width {kw:.6g} "
                            f"(winner by wall: {bf}).")
        lines.append("")
        lines.append("Win counts (best width within budget, ties go to first min):")
        lines.append("")
        for f, c in sorted(wins.items(), key=lambda x: -x[1]):
            lines.append(f"- `{f}`: {c}")
        lines.append("")
        k_wins = sum(c for f, c in wins.items() if f.startswith("kernel"))
        tot = sum(wins.values()) or 1
        if split == "heldout":
            findings_kernel.append(
                f"S5 held-out: kernel family wins {k_wins}/{tot} budgeted (instance,budget) slots "
                f"on certified p_width.")
        # accuracy loss check: any accepted row that misses reference?
        lean = [r for r in sub if (r.get("method") or {}).get("checker") in ("kernel", "v1", "v2")]
        misses = [r for r in lean if r.get("encloses_reference") is False]
        lines.append(f"Accepted Lean rows missing the exact reference on {split}: **{len(misses)}** "
                     f"(should be 0).")
        lines.append("")
        if split == "heldout" and not misses:
            findings_tie.append(
                "S5 held-out: every accepted Lean row with a reference encloses it "
                "(accuracy loss vs oracle = 0; disagreements = 0).")

    # ---- S6 ----
    s6 = load_rows("S6")
    lines.append("## S6 — Non-enumeration / scaling")
    lines.append("")
    lines.append("Protocol: kernel `uniform(d)` with `B` leaves on long / all-H sequences; "
                 "coverage = conformations / checker_steps. Matched v1/v2 B/M0 where feasible.")
    lines.append("")
    lines.append("| id | n | d | method | cert_bytes | leaves | verify_s | wall_s | "
                 "conf | steps | coverage |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(s6, key=lambda x: (x.get("n") or 0, x.get("depth") or 0, _method_key(x))):
        m = r.get("method") or {}
        lines.append(
            f"| {r.get('instance_id')} | {r.get('n')} | {r.get('depth')} | "
            f"{m.get('checker')}:{m.get('bounded')} | {r.get('cert_bytes')} | "
            f"{r.get('n_leaves')} | {_fmt(r.get('verify_s'))} | {_fmt(r.get('wall_s_median'))} | "
            f"{r.get('conformations')} | {r.get('checker_steps')} | {_fmt(r.get('coverage_ratio'))} |")
    lines.append("")
    big = [r for r in s6 if (r.get("method") or {}).get("checker") == "kernel"
           and r.get("coverage_ratio") and r["coverage_ratio"] > 10
           and r.get("accepted")]
    if big:
        r = max(big, key=lambda x: x["coverage_ratio"])
        findings_kernel.append(
            f"S6: kernel non-enumerative check covers ~{r['coverage_ratio']:.3g}× more "
            f"conformations than checker steps "
            f"({r.get('instance_id')}, d={r.get('depth')}, conf={r.get('conformations')}, "
            f"steps={r.get('checker_steps')}).")

    # ---- S7 ----
    s7 = load_rows("S7")
    lines.append("## S7 — One tree, many queries")
    lines.append("")
    lines.append("Protocol: one hybrid frontier; kernel carries N query lines in one check; "
                 "foldcheck2 is run once per (q,t) on the matched certificate; compare sum of walls.")
    lines.append("")
    lines.append("| n_queries | method | wall_s | verify_s | per_query |")
    lines.append("|---|---|---|---|---|")
    for r in sorted(s7, key=lambda x: (x.get("n_queries") or 0, _method_key(x))):
        m = r.get("method") or {}
        wall = r.get("wall_s_sum") if m.get("aggregate") else r.get("wall_s_median")
        lines.append(
            f"| {r.get('n_queries')} | {m.get('checker')} | {_fmt(wall)} | "
            f"{_fmt(r.get('verify_s'))} | {_fmt(r.get('per_query_verify_s'))} |")
    lines.append("")
    for nq in sorted({r.get("n_queries") for r in s7 if r.get("n_queries")}):
        kr = [r for r in s7 if r.get("n_queries") == nq
              and (r.get("method") or {}).get("checker") == "kernel" and r.get("wall_s_median")]
        vr = [r for r in s7 if r.get("n_queries") == nq
              and (r.get("method") or {}).get("checker") == "v2" and r.get("wall_s_sum")]
        if kr and vr:
            sp = vr[0]["wall_s_sum"] / kr[0]["wall_s_median"]
            findings_kernel.append(
                f"S7/n_queries={nq}: foldcheck2 sum-of-runs / kernel multi-query wall = {sp:.1f}× "
                f"({vr[0]['wall_s_sum']:.4g}s / {kr[0]['wall_s_median']:.4g}s).")

    # ---- S8 ----
    s8 = load_rows("S8")
    lines.append("## S8 — Footprint")
    lines.append("")
    if s8:
        r = s8[-1]
        lines.append("| binary | bytes | ratio vs foldkernel |")
        lines.append("|---|---|---|")
        bins = r.get("binaries") or {}
        bk = (bins.get("kernel") or {}).get("bytes") or 1
        for name in ("kernel", "v1", "v2"):
            b = (bins.get(name) or {}).get("bytes")
            lines.append(f"| fold{ {'kernel':'kernel','v1':'check','v2':'check2'}[name] } | "
                         f"{b:,} | {b / bk:.1f}× |")
        notes = r.get("build_notes_coordinator") or {}
        lines.append("")
        lines.append(f"Coordinator-provided build cost (not re-measured): "
                     f"foldcheck+foldcheck2 wall {notes.get('foldcheck_plus_foldcheck2_wall_s')}s / "
                     f"{notes.get('foldcheck_plus_foldcheck2_cpu_min')} CPU-min; "
                     f"foldkernel wall {notes.get('foldkernel_wall_s')}s.")
        if r.get("size_ratio_v1_over_kernel"):
            findings_kernel.append(
                f"S8: foldcheck binary is {r['size_ratio_v1_over_kernel']:.0f}× larger than "
                f"foldkernel; foldcheck2 is {r['size_ratio_v2_over_kernel']:.0f}× larger.")
    lines.append("")

    # ---- Findings block (placed after collection; rewrite header section) ----
    # We already started Findings earlier as a placeholder — splice properly:
    find_block = ["## Findings", ""]
    find_block.append("### (i) Kernel better")
    find_block += ([f"- {x}" for x in findings_kernel] or ["- (none in this run)"])
    find_block.append("")
    find_block.append("### (ii) Ties / no accuracy loss")
    find_block += ([f"- {x}" for x in findings_tie] or ["- (none in this run)"])
    find_block.append("")
    find_block.append("### (iii) An alternative better")
    find_block += ([f"- {x}" for x in findings_alt] or ["- (none in this run)"])
    find_block.append("")
    find_block.append("### (iv) Inconclusive / timeouts / skipped")
    find_block += ([f"- {x}" for x in findings_inconclusive] or [])
    if STATE.skipped:
        find_block += [f"- skipped: `{s}`" for s in STATE.skipped]
    if STATE.critical:
        find_block += [f"- CRITICAL finding saved: `{c.get('tag')}`" for c in STATE.critical]
    if not STATE.skipped and not STATE.critical and not findings_inconclusive:
        find_block.append("- (none in this run)")
    find_block.append("")

    # Replace the early stub "## Findings\n\n" with the full block
    text = "\n".join(lines)
    stub = "## Findings\n\n"
    if stub in text:
        text = text.replace(stub, "\n".join(find_block) + "\n", 1)
    else:
        text = "\n".join(find_block) + "\n" + text

    text += "\n## Critical findings\n\n"
    if STATE.critical:
        for c in STATE.critical:
            text += f"- `{c.get('tag')}` — see `bench/findings/`\n"
    else:
        text += "None. No accepted Lean certificate missed a frozen reference in this run.\n"
    text += "\n## Load / remaining work\n\n"
    text += f"- Load average at summary time: {bl.loadavg()}\n"
    text += f"- Skipped entries: {len(STATE.skipped)}\n"
    text += f"- Result JSONL paths: " + ", ".join(
        f"`experiments/kernel/bench/results/{s}.jsonl`" for s in SUITES) + "\n"
    return text


def _fmt(x) -> str:
    if x is None:
        return "—"
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        if abs(x) >= 100 or abs(x) == 0:
            return f"{x:.4g}"
        return f"{x:.6g}"
    return str(x)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main(argv: list[str]) -> int:
    global STATE
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--quick", action="store_true", help="DEV subset; ~10 min wall budget")
    g.add_argument("--full", action="store_true", help="all suites; ~75 min wall budget")
    g.add_argument("--summarize-only", action="store_true",
                   help="rebuild SUMMARY.md from existing JSONL rows")
    ap.add_argument("--only", default="",
                    help="comma-separated suite ids (S1..S8); default = all")
    ap.add_argument("--keep-rows", action="store_true",
                    help="append to existing JSONL files instead of clearing selected suites")
    args = ap.parse_args(argv)

    only = {s.strip().upper() for s in args.only.split(",") if s.strip()} or None
    if only:
        bad = only - set(SUITES)
        if bad:
            ap.error(f"unknown suites: {sorted(bad)}")

    if args.summarize_only:
        STATE = State("summarize", only, wall_budget=1e9)
        # reload binary metadata only
        text = summarize()
        RESULTS.mkdir(parents=True, exist_ok=True)
        (RESULTS / "SUMMARY.md").write_text(text, encoding="utf-8")
        print(f"wrote {RESULTS / 'SUMMARY.md'}")
        return 0

    mode = "quick" if args.quick else "full"
    wall = QUICK_WALL_S if args.quick else FULL_WALL_S
    STATE = State(mode, only, wall)

    for name, path in bl.EXE.items():
        if not path.exists():
            raise SystemExit(f"missing binary: {path}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    FINDINGS.mkdir(parents=True, exist_ok=True)
    bl.SCRATCH.mkdir(parents=True, exist_ok=True)

    plan = [
        ("S1", suite_s1),
        ("S2", suite_s2),
        ("S3", suite_s3_s4),  # also writes S4
        ("S5", suite_s5),
        ("S6", suite_s6),
        ("S7", suite_s7),
        ("S8", suite_s8),
    ]
    # S3 handler also fills S4; clear both when running S3
    for suite, fn in plan:
        if only and suite not in only:
            if not (suite == "S3" and ("S3" in only or "S4" in only)):
                continue
        if not STATE.can_start(1.0):
            STATE.skipped.append(f"{suite}:deadline_before_start")
            break
        if not args.keep_rows:
            clear_suite(suite)
            if suite == "S3":
                clear_suite("S4")
        print(f"=== {suite}  remaining={STATE.remaining():.0f}s  load={bl.loadavg()} ===",
              flush=True)
        try:
            fn(args.quick)
        except SystemExit:
            raise
        except Exception:
            traceback.print_exc()
            STATE.skipped.append(f"{suite}:exception")
            save_finding(f"{suite}_exception", {"traceback": traceback.format_exc()})

    text = summarize()
    (RESULTS / "SUMMARY.md").write_text(text, encoding="utf-8")
    # also dump a compact machine JSON
    (RESULTS / "run_meta.json").write_text(
        json.dumps({
            "mode": STATE.mode,
            "machine": STATE.machine,
            "rows_written": STATE.rows_written,
            "skipped": STATE.skipped,
            "critical": STATE.critical,
            "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "wall_s": time.perf_counter() - STATE.t0,
            "loadavg_end": bl.loadavg(),
        }, indent=1) + "\n", encoding="ascii")
    print(f"done: rows={STATE.rows_written} wall={time.perf_counter() - STATE.t0:.1f}s "
          f"summary={RESULTS / 'SUMMARY.md'}", flush=True)
    return 1 if STATE.critical else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
