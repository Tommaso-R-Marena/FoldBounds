#!/usr/bin/env python3
"""Render benchmark JSON (lean_v2_benchmark*.json) as Markdown before/after tables.

    python3 experiments/lean_v2_tables.py experiments/results/v2/lean_v2_benchmark_mixed.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction


def f(s):
    return float(Fraction(s))


def main(path: str) -> None:
    rows = json.load(open(path))
    print("| sequence | q | thr | d | mode | nodes | bytes | verify ms | Z interval | P interval | P width | exact Z | exact P | encloses |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        L = r["lean"]
        if L.get("accepted"):
            zi = f"[{f(L['zLo']):.4g}, {f(L['zHi']):.4g}]"
            pi = f"[{f(L['pLo']):.4g}, {f(L['pHi']):.4g}]"
            pw = f"{r['P_width']:.4g}"
            ms = L.get("verify_ms")
        else:
            zi = pi = pw = "—"
            ms = f"FAILED ({L.get('error')})"
        ez = f"{f(r['exact_Z']):.6g}" if "exact_Z" in r else "n/a"
        ep = f"{f(r['exact_P']):.4g}" if "exact_P" in r else "n/a"
        print(f"| {r['sequence']} | {r['q']} | {r['threshold']} | {r['depth']} | {r['mode']} | "
              f"{r['nodes']} | {r['cert_bytes']} | {ms} | {zi} | {pi} | {pw} | {ez} | {ep} | "
              f"{r.get('encloses_exact', 'n/a')} |")


if __name__ == "__main__":
    main(sys.argv[1])
