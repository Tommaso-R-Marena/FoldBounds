#!/usr/bin/env python3
"""Generate the frozen instance manifests (run ONCE; refuses to overwrite).

    python make_manifests.py            # writes manifest_dev.json, manifest_heldout.json
    python make_manifests.py --verify   # recompute everything and compare with the files

DEV       all-H n in {10,12,14,16}; alternating HP n in {12,16}; HPHHPHHPHHHH;
          HPHPPHHPHPPHPH; 6 random sequences from random.Random(20261117),
          one per (n, pH) in {12,14,16} x {0.5, 0.7}.
HELD-OUT  12 random sequences from random.Random(915273), one per (n, pH) in
          {13,15,17,18} x {0.4,0.6,0.8}; all-H n=18; alternating n=18; two long
          random sequences n=30 and n=48 (pH 0.6, same generator, drawn after
          the 12; no exact reference, scaling suite only).

q in {1, 3/2, 2, 3}.  Thresholds per sequence (n <= 18), from the exact
contact histogram at q = 2:
  t1   = 1
  mid  = the t in 1..kmax whose P(k >= t) is closest to 1/2 (ties: smaller t)
  high = the largest t >= 1 with P(k >= t) >= 1/100 (1 if there is none)
Long sequences: t in {1, 4, 8}.

Exact references (Z, Z_A, Z_notA, P as rational strings) are computed with
``foldkernel`` root ``X`` and cross-checked against the pure-Python enumerator
``hpenum.py`` (required for n <= 14; done up to --xcheck-max, default 18) and,
when the other agent's ``independent_oracle.py`` is importable, against that
too.  Any disagreement aborts.
"""
from __future__ import annotations

import argparse
import datetime
import json
import random
import sys
from fractions import Fraction

sys.dont_write_bytecode = True

import benchlib as bl
import certgen
import hpenum

Q_VALUES = [Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3)]
DEV_SEED = 20261117
HELDOUT_SEED = 915273
LONG_THRESHOLDS = [1, 4, 8]
KERNEL_TIMEOUT_S = 300.0


def random_sequence(rng: random.Random, n: int, ph: float) -> str:
    return "".join("H" if rng.random() < ph else "P" for _ in range(n))


def alternating(n: int) -> str:
    return ("HP" * n)[:n]


def dev_instances() -> list[dict]:
    out = []
    for n in (10, 12, 14, 16):
        out.append({"id": f"dev-allH-{n}", "sequence": "H" * n, "kind": "all-H"})
    for n in (12, 16):
        out.append({"id": f"dev-alt-{n}", "sequence": alternating(n), "kind": "alternating"})
    out.append({"id": "dev-fixed-12", "sequence": "HPHHPHHPHHHH", "kind": "fixed"})
    out.append({"id": "dev-fixed-14", "sequence": "HPHPPHHPHPPHPH", "kind": "fixed"})
    rng = random.Random(DEV_SEED)
    for n in (12, 14, 16):
        for ph in (0.5, 0.7):
            out.append({"id": f"dev-rand-n{n}-p{ph}", "sequence": random_sequence(rng, n, ph),
                        "kind": "random", "seed": DEV_SEED, "h_probability": ph})
    return out


def heldout_instances() -> list[dict]:
    out = []
    rng = random.Random(HELDOUT_SEED)
    for n in (13, 15, 17, 18):
        for ph in (0.4, 0.6, 0.8):
            out.append({"id": f"ho-rand-n{n}-p{ph}", "sequence": random_sequence(rng, n, ph),
                        "kind": "random", "seed": HELDOUT_SEED, "h_probability": ph})
    out.append({"id": "ho-allH-18", "sequence": "H" * 18, "kind": "all-H"})
    out.append({"id": "ho-alt-18", "sequence": alternating(18), "kind": "alternating"})
    for n in (30, 48):
        out.append({"id": f"ho-long-{n}", "sequence": random_sequence(rng, n, 0.6),
                    "kind": "long-random", "seed": HELDOUT_SEED, "h_probability": 0.6})
    return out


def run_kernel_full(seq: str, queries) -> dict:
    path = bl.write_scratch("manifest_ref.cert", certgen.kernel_full(seq, queries))
    res = bl.timed_run([str(bl.EXE["kernel"]), "check", str(path)], KERNEL_TIMEOUT_S)
    if res["timed_out"] or res["returncode"] != 0:
        raise SystemExit(f"foldkernel failed on {seq}: rc={res['returncode']} timeout={res['timed_out']}")
    rep = bl.parse_report("kernel", res["stdout"])
    if not rep["accepted"]:
        raise SystemExit(f"foldkernel rejected root X for {seq}: {rep['error']}")
    return rep


def choose_thresholds(hist: list[int]) -> tuple[dict, dict]:
    kmax = len(hist) - 1
    q2 = Fraction(2)
    p_by_t = {t: hpenum.masses(hist, q2, t)["P"] for t in range(1, kmax + 1)}
    half = Fraction(1, 2)
    mid = min(p_by_t, key=lambda t: (abs(p_by_t[t] - half), t)) if p_by_t else 1
    ok = [t for t, p in p_by_t.items() if p >= Fraction(1, 100)]
    high = max(ok) if ok else 1
    return ({"t1": 1, "mid": mid, "high": high},
            {str(t): bl.fstr(p) for t, p in p_by_t.items()})


def try_independent_oracle():
    """density_of_states from the other agent's oracle, if present (no cache use)."""
    path = bl.REPO / "experiments" / "kernel" / "independent_oracle.py"
    if not path.exists():
        return None
    import importlib.util
    try:
        spec = importlib.util.spec_from_file_location("independent_oracle", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.density_of_states
    except Exception as exc:                      # pragma: no cover - optional dependency
        print(f"independent_oracle.py present but unusable: {exc!r}", file=sys.stderr)
        return None


def complete(inst: dict, xcheck_max: int, oracle) -> dict:
    seq = inst["sequence"]
    n = len(seq)
    inst["n"] = n
    inst["h_count"] = seq.count("H")
    count, _ = bl.conformations(n)
    if inst["kind"] == "long-random":
        inst["thresholds"] = {"t_a": LONG_THRESHOLDS[0], "t_b": LONG_THRESHOLDS[1],
                              "t_c": LONG_THRESHOLDS[2]}
        inst["threshold_set"] = LONG_THRESHOLDS
        inst["reference"] = None
        inst["reference_source"] = "none (scaling suite only)"
        return inst
    # 1. exact histogram from the kernel
    rep = run_kernel_full(seq, [(Fraction(2), 1)])
    hist = rep["raw"]["exact_histogram"] or [0]
    if rep["raw"]["exact_count"] != count or sum(hist) != count:
        raise SystemExit(f"{inst['id']}: kernel count {rep['raw']['exact_count']} != c_(n-1)/4 = {count}")
    # 2. cross-checks
    checks = {"python_enumerator": None, "independent_oracle": None}
    if n <= xcheck_max:
        mine = hpenum.contact_histogram(seq)
        if mine != hist:
            raise SystemExit(f"{inst['id']}: hpenum histogram differs from foldkernel: {mine} vs {hist}")
        checks["python_enumerator"] = "histogram equal"
        if oracle is not None:
            theirs = oracle(seq)
            if list(theirs) != hist:
                raise SystemExit(f"{inst['id']}: independent_oracle differs from foldkernel")
            checks["independent_oracle"] = "histogram equal"
    if n <= 14 and checks["python_enumerator"] is None:
        raise SystemExit("cross-check with the Python enumerator is mandatory for n <= 14")
    # 3. thresholds and references
    thresholds, p_by_t = choose_thresholds(hist)
    tset = sorted(set(thresholds.values()))
    queries = [(q, t) for q in Q_VALUES for t in tset]
    rep = run_kernel_full(seq, queries)
    reference = {}
    for q, t in queries:
        key = bl.qkey(q, t)
        end = rep["results"][key]
        if not bl.is_exact(end):
            raise SystemExit(f"{inst['id']}: root X is not exact for {key}")
        m = hpenum.masses(hist, q, t)
        ref = {k: bl.fstr(v) for k, v in m.items()}
        if not bl.equals_reference(end, ref):
            raise SystemExit(f"{inst['id']}: kernel query result differs from histogram arithmetic at {key}")
        reference[key] = ref
    inst.update({
        "conformations": count,
        "histogram": hist,
        "kmax": len(hist) - 1,
        "thresholds": thresholds,
        "threshold_set": tset,
        "p_at_q2_by_threshold": p_by_t,
        "reference": reference,
        "reference_source": "foldkernel root X (exact histogram and query results)",
        "crosscheck": checks,
    })
    return inst


def build(name: str, xcheck_max: int) -> dict:
    oracle = try_independent_oracle()
    raw = dev_instances() if name == "dev" else heldout_instances()
    instances = []
    for inst in raw:
        print(f"[{name}] {inst['id']:<22} {inst['sequence']}", flush=True)
        instances.append(complete(inst, xcheck_max, oracle))
    return {
        "manifest": name,
        "model": "anchored 2D HP lattice; first bond (0,0)->(1,0); reflections not identified",
        "q_values": [bl.fstr(q) for q in Q_VALUES],
        "seed": DEV_SEED if name == "dev" else HELDOUT_SEED,
        "threshold_rule": "t1=1; mid=argmin_{1<=t<=kmax}|P_{q=2}(k>=t)-1/2| (ties: smaller t); "
                          "high=max{t>=1: P_{q=2}(k>=t)>=1/100} (1 if none); long: {1,4,8}",
        "rational_format": "num/den in lowest terms",
        "reference_key_format": "qnum/qden|t",
        "xcheck_max_n": xcheck_max,
        "independent_oracle_used": oracle is not None,
        "foldkernel_sha256": bl.sha256_file(bl.EXE["kernel"]),
        "instances": instances,
    }


def comparable(m: dict) -> dict:
    return {k: v for k, v in m.items() if k not in ("created_utc", "independent_oracle_used")} | {
        "instances": [{k: v for k, v in i.items() if k != "crosscheck"} for i in m["instances"]]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--verify", action="store_true", help="recompute and compare with the frozen files")
    ap.add_argument("--xcheck-max", type=int, default=18,
                    help="largest n cross-checked with the Python enumerator (>= 14)")
    args = ap.parse_args()
    if args.xcheck_max < 14:
        ap.error("--xcheck-max must be at least 14")
    status = 0
    for name in ("dev", "heldout"):
        path = bl.BENCH_DIR / f"manifest_{name}.json"
        if path.exists() and not args.verify:
            print(f"{path.name} exists: frozen, not regenerated (use --verify to re-check it)")
            continue
        manifest = build(name, args.xcheck_max)
        if args.verify:
            frozen = json.loads(path.read_text(encoding="ascii"))
            same = comparable(frozen) == comparable(manifest)
            print(f"{path.name}: {'identical to a fresh recomputation' if same else 'DIFFERS'}")
            status |= 0 if same else 1
        else:
            manifest["created_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
            path.write_text(json.dumps(manifest, indent=1) + "\n", encoding="ascii")
            print(f"wrote {path}")
    return status


if __name__ == "__main__":
    sys.exit(main())
