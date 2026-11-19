#!/usr/bin/env python3
"""Plain exhaustive enumeration of the anchored 2D HP lattice model.

UNVERIFIED exact baseline of the benchmark harness.  Standard library only;
imports nothing from the repository (no checker, no ``foldbounds``).

Model (docs/research/CERTIFIED_FRONTIER_KERNEL_SPEC.md, section 1):
self-avoiding walks of n sites on Z^2 starting (0,0),(1,0); reflections are
not identified; k(c) = number of non-consecutive H-H nearest-neighbour pairs.

CLI (used by run_bench.py as a separately timed child process):

    python hpenum.py SEQUENCE

prints one JSON object {"sequence", "n", "histogram", "conformations",
"enum_s"} where ``enum_s`` is the time spent inside the enumeration only.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction


def contact_histogram(seq: str) -> list[int]:
    """hist[k] = number of anchored conformations of ``seq`` with k contacts."""
    n = len(seq)
    if n < 2 or any(ch not in "HP" for ch in seq):
        raise ValueError("sequence must be an H/P string of length >= 2")
    is_h = [1 if ch == "H" else 0 for ch in seq]
    width = 2 * n + 3
    occ = bytearray(width * width)        # 1 = site occupied
    hh = bytearray(width * width)         # 1 = site occupied by an H residue
    hist = [0] * (2 * n + 2)
    start = (n + 1) * width + (n + 1)
    occ[start] = 1
    hh[start] = is_h[0]
    occ[start + 1] = 1
    hh[start + 1] = is_h[1]
    if n == 2:
        return [1]
    steps = (1, -1, width, -width)        # +x, -x, +y, -y
    last = n - 1

    def place(prev: int, j: int, k: int) -> None:
        """Put residue j on every free neighbour of ``prev`` (site of j-1)."""
        h = is_h[j]
        for st in steps:
            p = prev + st
            if occ[p]:
                continue
            if h:
                # H neighbours of p among placed residues, minus the chain
                # neighbour j-1 (consecutive residues never form a contact).
                kk = k + hh[p + 1] + hh[p - 1] + hh[p + width] + hh[p - width] - hh[prev]
            else:
                kk = k
            if j == last:
                hist[kk] += 1
            else:
                occ[p] = 1
                hh[p] = h
                place(p, j + 1, kk)
                occ[p] = 0
                hh[p] = 0

    limit = sys.getrecursionlimit()
    if limit < n + 100:
        sys.setrecursionlimit(n + 100)
    place(start + 1, 2, 0)
    while len(hist) > 1 and hist[-1] == 0:
        hist.pop()
    return hist


def masses(hist: list[int], q: Fraction, t: int) -> dict[str, Fraction]:
    """Exact Z, Z_A, Z_notA and P = Z_A / Z for weight q and event k >= t."""
    z = za = Fraction(0)
    w = Fraction(1)
    for k, c in enumerate(hist):
        term = c * w
        z += term
        if k >= t:
            za += term
        w *= q
    return {"Z": z, "Z_A": za, "Z_notA": z - za, "P": za / z}


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: hpenum.py SEQUENCE", file=sys.stderr)
        return 2
    t0 = time.perf_counter()
    hist = contact_histogram(argv[1])
    dt = time.perf_counter() - t0
    print(json.dumps({"sequence": argv[1], "n": len(argv[1]), "histogram": hist,
                      "conformations": sum(hist), "enum_s": dt}))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
