"""Reproduce exact all-q cubic 3D HP sequence design and its full geometry audit.

Usage:
  python experiments/cubic_3d_phase_demo.py --length 11
  python experiments/cubic_3d_phase_demo.py --length 12

The n=12 independent audit compiles the separate unsymmetrized C++ scorer
when g++ is available; otherwise it performs the slower full Python audit.
All root isolation/sequence selection uses exact fractions, not floats.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

from foldbounds.contact_map_atlas import build_atlas
from foldbounds.cubic_contact_atlas import (
    audit_cubic_atlas, audit_cubic_landscape,
    build_cubic_atlas, build_cubic_landscape
)
from foldbounds.exact_phase import certify_sequence_comparison
from foldbounds.full_skyline_phase import certify_complete_skyline
from foldbounds.universal_landscape import build_landscape


def independent_cpp_audit(atlas):
    """Compare EVERY contact-mask multiplicity against a second C++ enumerator."""
    if not shutil.which("g++"):
        return False
    source = Path(__file__).with_name("independent_cubic_audit.cpp")
    with TemporaryDirectory() as directory:
        executable = Path(directory) / "cubic_audit"
        subprocess.run(
            ["g++", "-O3", "-std=c++17", str(source), "-o", str(executable)],
            check=True, capture_output=True, text=True, timeout=45
        )
        run = subprocess.run(
            [str(executable), str(atlas.length)], check=True,
            capture_output=True, text=True, timeout=60
        )
        histogram = tuple(tuple(int(v) for v in line.split())
                          for line in run.stdout.splitlines())
    return histogram == atlas.multiplicities


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, choices=[11, 12], default=11)
    args = parser.parse_args()
    n = args.length
    atlas = build_cubic_atlas(n, max_nodes=3_000_000)
    print(f"Cubic 3D HP; n={n}, walks={atlas.conformations:,}, "
          f"contact maps={len(atlas.multiplicities):,}", flush=True)

    model = build_cubic_landscape(
        atlas, require_audit=False, max_length=n,
        max_edge_subsets=3_000_000
    )
    if n == 11:
        if not audit_cubic_landscape(model, atlas, max_length=11,
                                     max_comparisons=8_000_000):
            raise RuntimeError("independent full 3D walk/polynomial audit failed")
        print("Independent full 3D Python walk and all-2048-polynomial audit: PASS")
    else:
        if shutil.which("g++"):
            if not independent_cpp_audit(atlas):
                raise RuntimeError("independent C++ full contact histogram mismatch")
            print("Independent full 3D C++ contact histogram audit: PASS")
        elif not audit_cubic_atlas(atlas, max_length=12,
                                   max_leaves=8_000_000):
            raise RuntimeError("full independent Python geometry audit failed")
        else:
            print("Independent full 3D Python contact histogram audit: PASS")

        # After independent geometry validation, check all polynomial
        # distributions against direct contact-map scores, not the
        # positive graph-subset transform used by the builder.
        for mask in range(1 << n):
            seq = "".join("H" if mask & (1 << i) else "P"
                          for i in range(n))
            if model.contact_spectrum(seq) != atlas.spectrum(seq):
                raise RuntimeError("independent direct-spectrum audit mismatch")
        print("Independent direct-score spectra for all 4096 sequences: PASS")

    certificates = tuple(certify_complete_skyline(model, k)
                         for k in range(n + 1))
    transitions = sum(len(c.switches) for c in certificates)
    expected = 11 if n == 11 else 14
    assert transitions == expected
    assert sum(c.population for c in certificates) == (1 << n)
    print(f"All-q exact reference: {1 << n:,} sequences in {n + 1} H-compositions, "
          f"{transitions} global partition-optimum switches.")
    for cert in certificates:
        print(f"  H={cert.hydrophobic_count:2d} population={cert.population:5d} "
              f"skyline classes={len(cert.skyline_classes):2d} "
              f"all positive pair roots={len(cert.events):2d} "
              f"global switches={len(cert.switches):2d}")
    if n == 11:
        A, B = "HPPHPPHPPHP", "HPPHPHPPHPP"
        c3 = certify_sequence_comparison(
            model, A, B, width=Fraction(1, 1_000_000)
        )
        l2 = build_landscape(build_atlas(11))
        c2 = certify_sequence_comparison(l2, A, B)
        assert c3.positive_root_count == 2 and c2.positive_root_count == 0
        print("\nDimension-dependent exact A/B ordering:")
        print("  2D: A>B for every q>1")
        print("  3D: A>B, B>A, A>B; exactly two crossings at:")
        for root in c3.roots:
            print(f"    q in [{root.lo + 1}, {root.hi + 1}] "
                  f"(~{float((root.lo+root.hi)/2+1):.9f})")
    print("\nCaveat: exact finite cubic HP toy model; not a real protein prediction.")
    print("Python and C++ audit are independently implemented; not Lean-kernel certified.")


if __name__ == "__main__":
    main()
