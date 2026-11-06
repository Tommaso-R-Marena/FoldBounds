"""Contrast certified compactness with certified lack of exact fold confidence."""
from fractions import Fraction

from foldbounds.contact_vs_fold import certify_contact_vs_fold
from foldbounds.local_rewire_confidence import (
    StructuralModel, hp_pair_weights, make_macro_target
)


def show(dims, sequence, q, threshold):
    model = StructuralModel(
        sequence,
        make_macro_target(*dims),
        hp_pair_weights(Fraction(q)),
        Fraction(1, 2),
        Fraction(2),
    )
    c = certify_contact_vs_fold(model, dims, threshold)
    print(f"Model: 3D simple cubic HP, n={c.residues}, "
          f"H/P={c.hydrophobic_residues}/{c.polar_residues}, q={c.q}")
    print(f"  SAME exact model and per-state corrections [1/2,2]")
    print(f"  Event HH>={c.event_threshold}: "
          f"P(event) >= {float(c.contact_event_probability_lower):.12g}")
    print(f"  Explicit target HH={c.target_hh_contacts}")
    print(f"  Exact indexed target: "
          f"P(target) <= {float(c.target_microstate_probability_upper):.12g}")
    print(f"  Conservative rigidly congruent class: "
          f"P(class) <= {float(c.target_rigid_fold_probability_upper):.12g}")
    print(f"  Independent witness-family size (event): "
          f"{c.contact_certificate.family_paths:,}")
    print(f"  Independently rewired target microstates: "
          f"{c.structure_certificate.constructive_conformations:,}")
    print(f"  Both strict majority AND no unique fold majority: "
          f"{c.proves_event_majority_without_unique_fold_majority}")
    print("  THESE ARE TOY LATTICE PROBABILITIES, not real protein confidence.\n")


def main():
    show((3, 3, 3), "H"*216, 4, 108)
    seq = "".join("P" if i % 11 == 0 else "H" for i in range(216))
    show((3, 3, 3), seq, 5, 54)
    show((4, 4, 4), "H"*512, 4, 256)


if __name__ == "__main__":
    main()
