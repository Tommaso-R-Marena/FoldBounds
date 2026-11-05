"""Deterministic exact mutation scan; run with `python experiments/contact_map_atlas_demo.py`."""
from fractions import Fraction
from foldbounds.contact_map_atlas import audit_atlas, build_atlas


def main():
    sequence = "HHPPHHPHH"
    q = Fraction(2)
    atlas = build_atlas(len(sequence))
    if not audit_atlas(atlas):
        raise RuntimeError("independent full-enumeration audit failed")
    baseline = atlas.event_probability(sequence, q, 1)
    print(f"length={atlas.length} conformations={atlas.conformations} contact_maps={len(atlas.multiplicities)}")
    print(f"sequence={sequence} q={q} P(K>=1)={baseline} mean,var={atlas.contact_moments(sequence, q)}")
    for i, original in enumerate(sequence):
        mutant = sequence[:i] + ("P" if original == "H" else "H") + sequence[i+1:]
        result = atlas.event_probability(mutant, q, 1)
        print(f"residue={i+1} {original}->{mutant[i]} P(K>=1)={result} delta={result-baseline}")


if __name__ == "__main__":
    main()
