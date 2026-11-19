#!/usr/bin/env python3
"""Child-process worker for the UNVERIFIED Python bound methods of ``foldbounds``.

run_bench.py starts one fresh interpreter per configuration:

    python pymethods.py '<JSON spec>'

so that every configuration gets a hard timeout, its own peak-RSS reading and
the same process-start overhead as the compiled checkers.  The spec is

    {"sequence": "...", "q": "3/2", "t": 2, "method": "...", "params": {...}}

and the worker prints one JSON object with the eight exact endpoints, the time
spent in the frontier/search producer (``producer_s``) and the time spent in
the bound-checking function itself (``verify_s``).  Only documented entry
points are used (``python -m foldbounds --help``, docs/BOUND_HIERARCHY.md,
docs/ADAPTIVE_CERTIFICATES.md, experiments/compare_*.py):

method      entry point
----------  ------------------------------------------------------------------
basic       lattice.build_frontier(depth) + lattice.verify_certificate
lookahead   lookahead.verify_lookahead(cert, lookahead=k, use_matching=...,
            use_transfer=..., transfer_memory=..., use_monotone=...,
            use_turn_family=..., witness_nodes=..., terminal_box=...)
adaptive    adaptive.adaptive_refine(cert, steps=..., policy=..., lookahead=0,
            use_matching=True) then adaptive.verify_refinement_run (the replay
            is the checker; the greedy search is the producer)

Same model and the same eight quantities as the Lean checkers, but NOT
machine-checked.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def _endpoints(b) -> dict:
    def s(x: Fraction) -> str:
        return f"{x.numerator}/{x.denominator}"
    return {"zLo": s(b.partition_lower), "zHi": s(b.partition_upper),
            "eLo": s(b.event_lower), "eHi": s(b.event_upper),
            "nLo": s(b.complement_lower), "nHi": s(b.complement_upper),
            "pLo": s(b.probability_lower), "pHi": s(b.probability_upper)}


def run(spec: dict) -> dict:
    from foldbounds.lattice import build_frontier, verify_certificate
    from foldbounds.lookahead import verify_lookahead
    from foldbounds.adaptive import adaptive_refine, verify_refinement_run

    seq = spec["sequence"]
    q = Fraction(spec["q"])
    t = int(spec["t"])
    method = spec["method"]
    p = dict(spec.get("params", {}))
    depth = int(p.pop("depth", 2))

    t0 = time.perf_counter()
    cert = build_frontier(seq, q, t, depth)
    producer_s = time.perf_counter() - t0

    if method == "basic":
        t0 = time.perf_counter()
        bounds = verify_certificate(cert)
        verify_s = time.perf_counter() - t0
    elif method == "lookahead":
        t0 = time.perf_counter()
        bounds = verify_lookahead(cert, **p)
        verify_s = time.perf_counter() - t0
    elif method == "adaptive":
        steps = int(p.pop("steps"))
        policy = p.pop("policy")
        t0 = time.perf_counter()
        result = adaptive_refine(cert, steps=steps, policy=policy, lookahead=0,
                                 use_matching=True, **p)
        producer_s += time.perf_counter() - t0
        t0 = time.perf_counter()
        bounds = verify_refinement_run(result)
        verify_s = time.perf_counter() - t0
    else:
        raise ValueError(f"unknown method {method!r}")
    return {"accepted": True, "endpoints": _endpoints(bounds),
            "producer_s": producer_s, "verify_s": verify_s,
            "exact_paths": bounds.exact_paths, "bounded_subtrees": bounds.bounded_subtrees}


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: pymethods.py '<json spec>'", file=sys.stderr)
        return 2
    spec = json.loads(argv[1])
    try:
        out = run(spec)
    except Exception as exc:                       # reported as a row, never hidden
        out = {"accepted": None, "error": f"{type(exc).__name__}: {exc}"}
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
