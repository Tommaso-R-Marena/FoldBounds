"""Reproduce ground-state non-uniqueness in the 3D all-H HP lattice model."""
from fractions import Fraction

from foldbounds.cubic_groundstate_no_go import certify_cubic_groundstate_no_go


def main() -> None:
    for side in (4, 6, 8):
        result = certify_cubic_groundstate_no_go(
            side, q=Fraction(1000),
            correction_lower=Fraction(1, 2),
            correction_upper=Fraction(2),
        )
        c = result.confidence
        print(f"ALL-H 3D cubic model, L={side}, n={result.residues}")
        print(f"  EXACT GLOBAL maximum contacts: {result.global_contact_upper}")
        print(f"  Independently counted target HH: {result.target_contacts}")
        print(f"  Is the target a global ground state at q=1000? "
              f"{result.globally_optimal_for_attractive_hh}")
        print(f"  Distinct energy-equal global ground states (constructive lower): "
              f"{c.exactly_isoenergetic_conformations:,}")
        print(f"  Upper Gibbs probability for one exact state: "
              f"{float(c.microstate_probability_upper):.12g}")
        print(f"  Upper probability modulo any cubic rigid symmetry: "
              f"{float(c.rigid_fold_probability_upper):.12g}")
        print("  Robust to arbitrary per-state reweighting in [1/2,2]")
        print("  This is a toy model negative control, NOT a real protein prediction.\n")


if __name__ == "__main__":
    main()
