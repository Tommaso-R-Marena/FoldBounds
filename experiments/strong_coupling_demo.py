"""Reproduce two >50-residue model-relative robust majority certificates."""
from fractions import Fraction
from foldbounds.strong_coupling import (
    StrongCouplingCertificate, certify_majority_halfline,
    minimum_integer_majority_q, serpentine_witness,
)


def main() -> None:
    for width, height, threshold in ((10, 5, 20), (10, 10, 40), (50, 20, 500)):
        sequence = "H" * (width * height)
        path = serpentine_witness(width, height)
        cert = StrongCouplingCertificate(
            sequence, Fraction(1), threshold, (path,),
            correction_lower=Fraction(1, 2),
            correction_upper=Fraction(2))
        minimum = minimum_integer_majority_q(cert, max_q=100)
        assert minimum is not None
        result = certify_majority_halfline(
            StrongCouplingCertificate(
                sequence, Fraction(minimum), threshold, (path,),
                correction_lower=Fraction(1, 2),
                correction_upper=Fraction(2)))
        print("residues:", len(sequence), "threshold:", threshold)
        print("verified witnesses:", result.witness_count)
        print("universal anchored walk upper:", result.all_walk_count_upper)
        print("smallest certified INTEGER q:", minimum)
        print("guaranteed P lower for ALL real q>=q0:",
              float(result.probability_lower))
        print("uncertainty corrections: arbitrary independent [1/2, 2]")
        print("status:", "CERTIFIED MODEL-RELATIVE MAJORITY")
        print()


if __name__ == "__main__":
    main()
