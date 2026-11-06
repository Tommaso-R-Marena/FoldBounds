#!/usr/bin/env python3
"""Before/after benchmark: version-1 verified leaves vs version-2 leaves.

For every case the SAME untrusted frontier (Python ``build_frontier`` at a
given depth) and the SAME witness set are serialised in several leaf modes and
checked by the compiled Lean runners:

  v1-B   legacy basic bound            (foldcheck,  version-1 format)
  v1-T   legacy transfer bound         (foldcheck,  version-1 format)
  v1-TW  legacy transfer + witnesses   (foldcheck,  version-1 format)
  M0/M2/M4   memory-k transfer + geometric caps + PD family   (foldcheck2)
  M2W/M4W    the same with the witnesses used by v1-TW        (foldcheck2)

The version-1 texts are also re-checked by ``foldcheck2`` (legacy leaves) and
the results must be bit-identical (backward compatibility).  When an exact
oracle is affordable (independent Python enumeration ``exact_spectrum``) every
accepted interval is checked to enclose the exact Z, event mass, complement
mass and probability.

THIS IS MEASUREMENT, NOT PROOF.  Soundness of every accepted interval is the
Lean theorem ``FoldBounds.FrontierNext.checkFlat2_sound`` (resp.
``FrontierCheck.checkFlat_sound`` for version 1).

Usage (Lake project root):
    lake build foldcheck foldcheck2
    python3 FoldBounds-main/experiments/lean_v2_benchmark.py --lake-root . \
        --exe1 .lake/build/bin/foldcheck --exe2 .lake/build/bin/foldcheck2
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from fractions import Fraction
from pathlib import Path as FsPath

HERE = FsPath(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

from foldbounds.lattice import build_frontier, exact_spectrum  # noqa: E402
from foldbounds.lean_export import certificate_to_lean_text  # noqa: E402
from foldbounds.lean_export_v2 import certificate_to_lean_text_v2  # noqa: E402
from foldbounds.witnesses import propose_dfs_witnesses  # noqa: E402


def frac(s: str) -> Fraction:
    a, b = s.split("/")
    return Fraction(int(a), int(b))


def _limit():
    gb = 12 * 1024 ** 3
    resource.setrlimit(resource.RLIMIT_AS, (gb, gb))


def run_lean(exe: str, lake_root: str, text: str, timeout: float) -> dict:
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "cert.txt")
        FsPath(f).write_text(text)
        t0 = time.time()
        try:
            res = subprocess.run([exe, "check", f], cwd=lake_root, capture_output=True, text=True,
                                 timeout=timeout, preexec_fn=_limit)
        except subprocess.TimeoutExpired:
            return {"accepted": None, "error": "timeout", "wall_s": time.time() - t0}
        wall = time.time() - t0
    if res.returncode != 0 or not res.stdout.strip():
        return {"accepted": None, "error": f"exit {res.returncode}", "wall_s": wall}
    row = json.loads(res.stdout.strip().splitlines()[-1])
    row["wall_s"] = wall
    return row


def exact_values(seq: str, q: Fraction, thr: int):
    spec = exact_spectrum(seq)
    z = sum((c * q ** k for k, c in spec.items()), Fraction(0))
    za = sum((c * q ** k for k, c in spec.items() if k >= thr), Fraction(0))
    return {"Z": z, "A": za, "N": z - za, "P": za / z, "count": sum(spec.values())}


def encloses(r: dict, ex: dict) -> bool:
    g = {k: frac(r[k]) for k in ("zLo", "zHi", "eLo", "eHi", "nLo", "nHi", "pLo", "pHi")}
    return (g["zLo"] <= ex["Z"] <= g["zHi"] and g["eLo"] <= ex["A"] <= g["eHi"]
            and g["nLo"] <= ex["N"] <= g["nHi"] and g["pLo"] <= ex["P"] <= g["pHi"])


CASES_DEFAULT = (
    [("H" * 10, Fraction(2), 3, d) for d in (4, 6)]
    + [("H" * 12, Fraction(2), 4, d) for d in (4, 6)]
    + [("H" * 14, Fraction(2), 4, d) for d in (6, 8)]
    + [("H" * 16, Fraction(2), 5, d) for d in (6, 8)]
    + [("H" * 18, Fraction(2), 6, d) for d in (8, 10)]
    + [("HPHHPHHPHHHH", q, thr, d) for q in (Fraction(1), Fraction(3, 2), Fraction(2))
       for thr in (2, 4) for d in (4, 6)]
    + [("HPHPPHHPHPPHPH", Fraction(2), 2, d) for d in (6, 8)]
    + [("HPHPPHHPHPPHPHHPPH", Fraction(2), 2, d) for d in (8, 10)]
)

CASES_MIXED = (
    [("HPHHPHHPHHHH", q, thr, d) for q in (Fraction(1), Fraction(3, 2), Fraction(2))
     for thr in (2, 4) for d in (4, 6)]
    + [("HPHPPHHPHPPHPH", Fraction(2), 2, d) for d in (6, 8)]
    + [("HPHPPHHPHPPHPHHPPH", Fraction(2), 2, 10)]
)

CASES_QUICK = [("H" * 10, Fraction(2), 3, 4), ("HPHHPHHPHHHH", Fraction(3, 2), 2, 4)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lake-root", default=".")
    ap.add_argument("--exe1", default=".lake/build/bin/foldcheck")
    ap.add_argument("--exe2", default=".lake/build/bin/foldcheck2")
    ap.add_argument("--out", default=str(HERE / "results" / "v2"))
    ap.add_argument("--timeout", type=float, default=900)
    ap.add_argument("--exact-max", type=int, default=16)
    ap.add_argument("--cases", default="default")
    ap.add_argument("--json-name", default="lean_v2_benchmark.json")
    ap.add_argument("--modes", default="v1-B,v1-T,v1-TW,M0,M2,M4,M2W,M4W")
    args = ap.parse_args()
    cases = {"quick": CASES_QUICK, "mixed": CASES_MIXED}.get(args.cases, CASES_DEFAULT)
    out = FsPath(args.out)
    out.mkdir(parents=True, exist_ok=True)
    out_name = args.json_name
    modes = args.modes.split(",")
    rows = []
    exact_cache: dict = {}
    for seq, q, thr, depth in cases:
        t0 = time.time()
        cert = build_frontier(seq, q, thr, depth)
        t_build = time.time() - t0
        t0 = time.time()
        wit = {node.path: propose_dfs_witnesses(node.path, len(seq), node_budget=60)
               for node in cert.frontier if node.kind == "bound"}
        t_wit = time.time() - t0
        ex = None
        if len(seq) <= args.exact_max:
            key = (seq, q, thr)
            if key not in exact_cache:
                exact_cache[key] = exact_values(seq, q, thr)
            ex = exact_cache[key]
        for mode in modes:
            if mode.startswith("v1-"):
                m = mode[3:]
                text = certificate_to_lean_text(cert, transfer=(m != "B"),
                                                witnesses=wit if m == "TW" else None)
                r = run_lean(args.exe1, args.lake_root, text, args.timeout)
                r2 = run_lean(args.exe2, args.lake_root, text, args.timeout)
                keys = ("accepted", "zLo", "zHi", "eLo", "eHi", "nLo", "nHi", "pLo", "pHi")
                compat = all(r.get(k) == r2.get(k) for k in keys)
            else:
                k = int(mode[1])
                w = wit if mode.endswith("W") else None
                text = certificate_to_lean_text_v2(cert, leaf="M", memory=k, witnesses=w)
                r = run_lean(args.exe2, args.lake_root, text, args.timeout)
                compat = None
            row = dict(sequence=seq, n=len(seq), q=str(q), threshold=thr, depth=depth, mode=mode,
                       nodes=len(cert.frontier), cert_bytes=len(text.encode()),
                       producer_s=t_build + (t_wit if mode.endswith("W") else 0.0),
                       v1_v2_compatible=compat, lean=r)
            if ex:
                row.update(exact_Z=str(ex["Z"]), exact_P=str(ex["P"]), conformations=ex["count"])
            if r.get("accepted"):
                zlo, zhi, plo, phi = (frac(r[x]) for x in ("zLo", "zHi", "pLo", "pHi"))
                row.update(Z_ratio=float(zhi / zlo) if zlo > 0 else None, P_width=float(phi - plo))
                if ex:
                    row["encloses_exact"] = encloses(r, ex)
            rows.append(row)
            print(f"{seq:>18} q={str(q):>3} thr={thr} d={depth} {mode:>5} nodes={len(cert.frontier):>5} "
                  f"bytes={row['cert_bytes']:>7} ms={r.get('verify_ms')} "
                  + (f"Z=[{float(frac(r['zLo'])):.4g},{float(frac(r['zHi'])):.4g}] "
                     f"P=[{float(frac(r['pLo'])):.4g},{float(frac(r['pHi'])):.4g}]"
                     if r.get("accepted") else f"FAILED {r.get('error')}")
                  + (f" exactZ={float(ex['Z']):.4g} P={float(ex['P']):.4g} ok={row.get('encloses_exact')}"
                     if ex else "")
                  + ("" if compat in (None, True) else " V1/V2 MISMATCH"), flush=True)
            (out / out_name).write_text(json.dumps(rows, indent=1))
    (out / out_name).write_text(json.dumps(rows, indent=1))
    bad = [r for r in rows if r.get("encloses_exact") is False or r.get("v1_v2_compatible") is False
           or r["lean"].get("accepted") is False]
    print("rows:", len(rows), "problems:", len(bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
