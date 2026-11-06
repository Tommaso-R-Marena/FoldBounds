#!/usr/bin/env python3
"""Cost/precision trade-off of hybrid version-2 certificates.

The untrusted proposer builds one frontier and chooses, leaf by leaf, between
a non-enumerative memory-k bound (``M<k>``) and exhaustive local evaluation
inside the verifier (``X``).  Policy (importance by certified prefix
contacts): a leaf is evaluated exhaustively iff its prefix already has at
least ``kappa`` H-H contacts (these subtrees carry most of the event mass).
Sweeping ``kappa`` from "never" to "always" traces the verifier-cost versus
probability-width curve between the pure bound certificate and complete
enumeration (``kappa = 0``: every leaf exhaustive = full enumeration).

Measurement, not proof (soundness: ``FrontierNext.checkFlat2_sound``).
"""
from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path as FsPath

HERE = FsPath(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))
sys.path.insert(0, str(HERE))

from foldbounds.lattice import _contacts, _possible_hh_pairs, build_frontier  # noqa: E402
from foldbounds.lean_export import path_to_word  # noqa: E402
from foldbounds.lean_export_v2 import HEADER_V2, _fraction, _suffixes  # noqa: E402
from foldbounds.witnesses import propose_dfs_witnesses  # noqa: E402
from lean_v2_benchmark import encloses, exact_values, frac, run_lean  # noqa: E402


def hybrid_text(cert, kappa: int | None, memory: int, wit) -> tuple[str, int]:
    pairs = _possible_hh_pairs(cert.sequence)
    lines = [HEADER_V2, f"sequence {cert.sequence}", f"q {_fraction(cert.q)}",
             f"threshold {cert.event_at_least}"]
    nx = 0
    for node in cert.frontier:
        w = path_to_word(node.path)
        if node.kind == "exact":
            lines.append(f"E {w}")
        elif node.kind == "dead":
            lines.append(f"D {w}")
        elif kappa is not None and _contacts(node.path, pairs) >= kappa:
            lines.append(f"X {w}")
            nx += 1
        else:
            lines.append(" ".join([f"M{memory}", w, *_suffixes(w, wit.get(node.path, ()))]))
    lines.append("end")
    return "\n".join(lines) + "\n", nx


CASES = [("H" * 12, Fraction(2), 4, 6), ("HPHHPHHPHHHH", Fraction(2), 4, 6),
         ("HPHHPHHPHHHH", Fraction(3, 2), 2, 6), ("HPHPPHHPHPPHPH", Fraction(2), 2, 8),
         ("H" * 14, Fraction(2), 4, 8)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lake-root", default=".")
    ap.add_argument("--exe2", default=".lake/build/bin/foldcheck2")
    ap.add_argument("--out", default=str(HERE / "results" / "v2"))
    ap.add_argument("--memory", type=int, default=4)
    ap.add_argument("--timeout", type=float, default=1800)
    args = ap.parse_args()
    rows = []
    for seq, q, thr, depth in CASES:
        cert = build_frontier(seq, q, thr, depth)
        wit = {nd.path: propose_dfs_witnesses(nd.path, len(seq), node_budget=60)
               for nd in cert.frontier if nd.kind == "bound"}
        ex = exact_values(seq, q, thr)
        nb = sum(1 for nd in cert.frontier if nd.kind == "bound")
        for kappa in [None, 4, 3, 2, 1, 0]:
            text, nx = hybrid_text(cert, kappa, args.memory, wit)
            r = run_lean(args.exe2, args.lake_root, text, args.timeout)
            row = dict(sequence=seq, q=str(q), threshold=thr, depth=depth, memory=args.memory,
                       kappa=kappa, bound_leaves=nb, exhaustive_leaves=nx,
                       cert_bytes=len(text.encode()), lean=r, exact_P=str(ex["P"]))
            if r.get("accepted"):
                plo, phi = frac(r["pLo"]), frac(r["pHi"])
                row.update(P_width=float(phi - plo), encloses_exact=encloses(r, ex))
            rows.append(row)
            print(f"{seq:>16} q={str(q):>3} thr={thr} d={depth} kappa={kappa} X={nx}/{nb} "
                  f"ms={r.get('verify_ms')} P=[{float(frac(r['pLo'])):.4f},{float(frac(r['pHi'])):.4f}] "
                  f"exactP={float(ex['P']):.4f} ok={row.get('encloses_exact')}" if r.get("accepted")
                  else f"{seq} kappa={kappa} FAILED {r}", flush=True)
    out = FsPath(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "lean_v2_hybrid.json").write_text(json.dumps(rows, indent=1))
    bad = [r for r in rows if r.get("encloses_exact") is False or not r["lean"].get("accepted")]
    print("rows:", len(rows), "problems:", len(bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
