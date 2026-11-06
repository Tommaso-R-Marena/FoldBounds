#!/usr/bin/env python3
"""Benchmark of the Lean compact frontier checker (B / T / TW modes).

For each (sequence, frontier depth) the untrusted Python producer builds a
depth-limited frontier certificate, optionally proposes witnesses, and the
*Lean* checker (compiled ``foldcheck`` executable, or the interpreter) verifies
it.  Reported per run:

* frontier size (number of nodes) and certificate size in bytes;
* producer time (Python, untrusted) and verifier time (Lean, measured inside
  the Lean process around the verified ``checkFlat`` call) plus peak RSS of the
  verifier process;
* the certified intervals for Z and P(k >= threshold) and their widths;
* the exact values when an independent exhaustive enumeration is affordable.

Modes:  B  = uniform cap  3^(n-m) * q^Kmax  per bound node (Python-identical)
        T  = verified transfer relaxation (fixed-prefix-avoiding walks,
             cross-contact DP, degree-based future caps)
        TW = T plus DFS witnesses (lower bounds) for every bound node

THIS IS MEASUREMENT, NOT PROOF.  Soundness of every accepted interval follows
from the Lean theorem ``FoldBounds.FrontierCheck.checkFlat_sound``; the timings
are empirical.

Usage (from the Lake project root):
    lake build foldcheck
    python3 FoldBounds-main/experiments/lean_frontier_benchmark.py \
        --lake-root . --exe .lake/build/bin/foldcheck --out FoldBounds-main/experiments/results
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

from foldbounds.lattice import build_frontier, exact_spectrum, verify_certificate  # noqa: E402
from foldbounds.lean_export import certificate_to_lean_text  # noqa: E402
from foldbounds.witnesses import propose_dfs_witnesses  # noqa: E402


def frac(s: str) -> Fraction:
    a, b = s.split("/")
    return Fraction(int(a), int(b))


def run_lean(args, text: str):
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "cert.txt")
        FsPath(f).write_text(text)
        before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        t0 = time.time()
        if args.exe:
            cmd = [args.exe, "check", f]
        else:
            cmd = ["lake", "env", "lean", "--run", args.runner, "check", f]
        res = subprocess.run(cmd, cwd=args.lake_root, capture_output=True, text=True,
                             check=True, timeout=args.timeout)
        wall = time.time() - t0
        rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    row = json.loads(res.stdout.strip().splitlines()[-1])
    row["wall_s"] = wall
    row["max_rss_kb_upper"] = max(rss, before)  # children high-water mark (upper bound)
    return row


def exact_Z_P(seq: str, q: Fraction, thr: int):
    spec = exact_spectrum(seq)
    z = sum((c * q ** k for k, c in spec.items()), Fraction(0))
    za = sum((c * q ** k for k, c in spec.items() if k >= thr), Fraction(0))
    return z, za / z, sum(spec.values())


def fmt(x: Fraction) -> str:
    return f"{float(x):.4g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lake-root", default=".")
    ap.add_argument("--exe", default="")
    ap.add_argument("--runner", default="FoldBounds-main/scripts/foldcheck.lean")
    ap.add_argument("--out", default=str(HERE / "results"))
    ap.add_argument("--timeout", type=float, default=1800)
    ap.add_argument("--exact-max", type=int, default=14)
    ap.add_argument("--cases", default="default")
    args = ap.parse_args()

    q = Fraction(2)
    if args.cases == "quick":
        cases = [("H" * 10, 4, 3), ("HPHPPHHPHH", 4, 2)]
    else:
        cases = []
        for n, depths in [(10, [2, 4, 6]), (12, [2, 4, 6]), (14, [2, 4, 6, 8]),
                          (16, [2, 4, 6, 8]), (18, [2, 4, 6, 8])]:
            for d in depths:
                cases.append(("H" * n, d, n // 4))
        for seq in ["HPHPPHHPHPPHPHHPPHPH"[:14], "HPHPPHHPHPPHPHHPPH"]:
            for d in [2, 4, 6, 8]:
                cases.append((seq, d, 2))

    rows = []
    for seq, depth, thr in cases:
        t0 = time.time()
        cert = build_frontier(seq, q, thr, depth)
        t_build = time.time() - t0
        t0 = time.time()
        py = verify_certificate(cert)
        t_pyverify = time.time() - t0
        t0 = time.time()
        wit = {node.path: propose_dfs_witnesses(node.path, len(seq), node_budget=60)
               for node in cert.frontier if node.kind == "bound"}
        t_wit = time.time() - t0
        exact = exact_Z_P(seq, q, thr) if len(seq) <= args.exact_max else None
        kinds = {k: sum(1 for x in cert.frontier if x.kind == k) for k in ("exact", "bound", "dead")}
        for mode in ("B", "T", "TW"):
            text = certificate_to_lean_text(cert, transfer=(mode != "B"),
                                            witnesses=wit if mode == "TW" else None)
            try:
                r = run_lean(args, text)
            except subprocess.TimeoutExpired:
                r = {"accepted": None, "timeout": True}
            row = dict(sequence=seq, n=len(seq), depth=depth, threshold=thr, q=str(q), mode=mode,
                       nodes=len(cert.frontier), kinds=kinds, cert_bytes=len(text.encode()),
                       producer_s=t_build + (t_wit if mode == "TW" else 0.0),
                       python_verify_s=t_pyverify, lean=r)
            if exact:
                row["exact_Z"], row["exact_P"], row["num_conformations"] = (
                    str(exact[0]), str(exact[1]), exact[2])
            if r.get("accepted"):
                zlo, zhi, plo, phi = (frac(r[k]) for k in ("zLo", "zHi", "pLo", "pHi"))
                row.update(zLo=str(zlo), zHi=str(zhi), pLo=str(plo), pHi=str(phi),
                           Z_ratio=float(zhi / zlo) if zlo > 0 else None,
                           P_width=float(phi - plo))
                if exact:
                    ok = zlo <= exact[0] <= zhi and plo <= exact[1] <= phi
                    row["encloses_exact"] = ok
                if mode == "B":
                    row["python_match"] = (zlo, zhi, plo, phi) == (
                        py.partition_lower, py.partition_upper,
                        py.probability_lower, py.probability_upper)
            rows.append(row)
            print(f"{seq:>20} d={depth} thr={thr} {mode:>2} nodes={len(cert.frontier):>5} "
                  f"bytes={row['cert_bytes']:>7} verify_ms={r.get('verify_ms')} "
                  f"Z=[{fmt(frac(r['zLo'])) if r.get('accepted') else '-'}, "
                  f"{fmt(frac(r['zHi'])) if r.get('accepted') else '-'}] "
                  f"P=[{fmt(frac(r['pLo'])) if r.get('accepted') else '-'}, "
                  f"{fmt(frac(r['pHi'])) if r.get('accepted') else '-'}]"
                  + (f" exactZ={fmt(exact[0])} encl={row.get('encloses_exact')}" if exact else ""),
                  flush=True)
    out = FsPath(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "lean_frontier_benchmark.json").write_text(json.dumps(rows, indent=1))
    bad = [r for r in rows if r.get("encloses_exact") is False or r.get("python_match") is False
           or r["lean"].get("accepted") is False]
    print("rows:", len(rows), "problems:", len(bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
