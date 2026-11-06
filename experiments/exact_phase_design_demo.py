"""Reproduce two exact phase transitions and a complete n=14 design diagram.

Usage: python experiments/exact_phase_design_demo.py

These are full-partition and event-probability statements in the anchored 2D
HP lattice model only. q is a dimensionless contact weight, not Kelvin.
"""
from fractions import Fraction

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.exact_phase import certify_sequence_comparison
from foldbounds.global_phase_diagram import certify_global_phase_diagram


def main():
    a12, b12 = "PPPPPHHHPPHP", "PHPPHHPPPPHP"
    atlas12 = build_atlas(12)
    assert audit_atlas(atlas12)
    land12 = build_landscape(atlas12)
    result = certify_sequence_comparison(
        land12, a12, b12, event_at_least=1, width=Fraction(1, 100_000)
    )
    assert result.positive_root_count == 2
    print("Length 12: exactly two event-probability crossovers.")
    print("A", a12, "spectrum", land12.contact_spectrum(a12))
    print("B", b12, "spectrum", land12.contact_spectrum(b12))
    for i, root in enumerate(result.roots, 1):
        print(f"  crossover {i}: q in [{root.lo + 1}, {root.hi + 1}]")
    print("signs of P_A - P_B:", result.cell_signs)

    atlas14 = build_atlas(14)
    assert audit_atlas(atlas14, max_length=14)
    model = build_landscape(atlas14, require_audit=False)
    d = certify_global_phase_diagram(model, 6, width=Fraction(1, 100_000))
    assert d.population == 3003 and len(d.candidates) == 2
    assert d.comparison is not None
    assert d.comparison.positive_root_count == 2
    print("\nLength 14: global maximum Z among all 3,003 six-H sequences.")
    for seq in d.candidates:
        print("Candidate", seq, "coefficients", model.polynomial(seq))
    for i, root in enumerate(d.comparison.roots, 1):
        print(f"  global switch {i}: q in [{root.lo + 1}, {root.hi + 1}]")
    print("globally optimal sequences in consecutive root-free q-cells:")
    for cell, winners in enumerate(d.cell_winners, 1):
        print(f"  cell {cell}: {winners}")
    print("Assurance: exact rational Python + independent full-walk geometry audit.")
    print("Not a Lean proof or experimental protein-prediction result.")


if __name__ == "__main__":
    main()
