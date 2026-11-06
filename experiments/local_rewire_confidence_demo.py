"""Reproduce exact structural-confidence ceilings for cubic lattice HP models."""
from fractions import Fraction

from foldbounds.local_rewire_confidence import (
    StructuralModel, check_structural_confidence, hp_pair_weights,
    make_macro_target, minimum_selectivity_for_confidence,
)


def run_case(dims, sequence, q: Fraction) -> None:
    path = make_macro_target(*dims)
    assert len(sequence) == len(path)
    model = StructuralModel(
        sequence, path, hp_pair_weights(q),
        correction_lower=Fraction(1, 2),
        correction_upper=Fraction(2),
    )
    result = check_structural_confidence(model)
    print(f"3D HP conformer, n={result.residues}, grid={dims}, "
          f"polar fraction={sequence.count('P')}/{len(sequence)}, q={q}")
    print(f"  completely verified 2x2x2 candidate windows: {result.candidate_cubes}")
    print(f"  nonadjacent factorized patch certificates: {len(result.selected_cubes)}")
    print(f"  independently constructible coordinate-distinct microstates: "
          f"{result.constructive_conformations:,}")
    print(f"  proven exactly isoenergetic states: "
          f"{result.exactly_isoenergetic_conformations:,}")
    print(f"  relative family Boltzmann partition >= "
          f"{result.relative_partition_lower.numerator}/"
          f"{result.relative_partition_lower.denominator}")
    print(f"  exact TARGET microstate probability <= "
          f"{float(result.microstate_probability_upper):.12g}")
    print(f"  rigid-orbit fold confidence ceiling (even with 48 symmetries) <= "
          f"{float(result.rigid_fold_probability_upper):.12g}")
    print(f"  robustness: arbitrary statewise weight corrections in [1/2,2]")
    print(f"  worst-case uniform target/competitor selectivity "
          f"needed for P>=0.9 among isoenergy states: "
          f"{minimum_selectivity_for_confidence(result.exactly_isoenergetic_conformations, Fraction(9,10)):,}")
    print("  IMPORTANT: this is a lattice-model identifiability limitation, "
          "NOT a real protein equilibrium estimate.\n")


def main() -> None:
    for dims in ((3, 3, 3), (4, 4, 4)):
        path = make_macro_target(*dims)
        run_case(dims, "H" * len(path), Fraction(5))
        sequence = "".join(
            "P" if i % 11 == 0 else "H" for i in range(len(path))
        )
        run_case(dims, sequence, Fraction(5))
        balanced = "".join(
            "P" if i % 2 == 0 else "H" for i in range(len(path))
        )
        run_case(dims, balanced, Fraction(3))


if __name__ == "__main__":
    main()
