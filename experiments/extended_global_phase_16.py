"""Reproduce exact n=16, H=8 global partition phase diagram for all real q>1.

Requires only the FoldBounds package and standard Python, no external CAS.
The result is relative to the anchored 2D HP lattice model, not a real
protein folding prediction. q is dimensionless, not temperature in Kelvin.
"""
from fractions import Fraction

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.extended_global_phase import certify_extended_global_diagram


def main():
    atlas = build_atlas(16)
    assert audit_atlas(atlas, max_length=16)
    model = build_landscape(
        atlas, require_audit=False, max_length=16,
        max_edge_subsets=3_000_000
    )
    cert = certify_extended_global_diagram(
        model, 8, width=Fraction(1, 100_000)
    )
    assert cert.population == 12870
    assert cert.coefficientwise_classes == 5
    assert len(cert.coefficientwise_survivors) == 9
    assert cert.crossover is not None and cert.crossover.positive_root_count == 1

    print("FoldBounds exact all-q global partition design certificate")
    print("Geometry: length 16, 1,604,149 anchored walks, 31,940 contact maps")
    print("Fixed composition: eight H, eight P; sequences:", cert.population)
    print("Coefficientwise survivors:", len(cert.coefficientwise_survivors))
    print("Distinct candidate polynomials:", cert.coefficientwise_classes)
    print("Sturm-undominated final polynomial classes:",
          len(cert.globally_relevant_classes))
    print("Global maximum sequence classes in consecutive open q-cells:")
    for i, seqs in enumerate(cert.cell_winners, 1):
        print("  cell", i, ":", ", ".join(seqs))
    for root in cert.crossover.roots:
        print("The ONLY global switching point has q in:",
              f"[{root.lo + 1}, {root.hi + 1}]")
        print("Approximate interval:", float(root.lo + 1), float(root.hi + 1))
    print("Arithmetic: exact integer/rational Python, full independent walk audit")
    print("Not Lean-refined; not an experimental protein result")


if __name__ == "__main__":
    main()
