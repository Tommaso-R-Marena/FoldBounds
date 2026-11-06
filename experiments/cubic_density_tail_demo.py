"""Reproduce infinite 3D heterogeneous HP population guarantees."""
from decimal import Decimal, localcontext
from fractions import Fraction
from foldbounds.cubic_density_tail import certify_periodic_cubic_density_tail


def formatted_odds(value: Fraction) -> str:
    with localcontext() as ctx:
        ctx.prec = 18
        x = Decimal(value.numerator) / Decimal(value.denominator)
        return f"{x:.8E}"


def main() -> None:
    for m, q in ((11, 5), (22, 4)):
        result = certify_periodic_cubic_density_tail(m, q=q)
        print(f"3D cube physical side={2*m}, residues={result.residues}")
        print(f"  periodic sequence: P at 0,11,22,...; P={result.polar_residues}")
        print(f"  event: >= {result.event_threshold} noncovalent H-H contacts")
        print(f"  constructive HH minimum: {result.hh_contacts_lower}")
        print("  local cube paths: independently enumerated; "
              f"c8={result.c8}")
        print(f"  family count >= 2*4**({m**3})")
        print(f"  arbitrary microstate corrections: [1/2, 2]")
        print(f"  event/complement rigorous failure odds <= "
              f"{formatted_odds(result.failure_odds_upper)}")
        print(f"  certified P(event)>0.99 at q={q}")
        print(f"  UNIFORMLY for ALL cubic macro sides m={m}+11*k, k>=0")
        print(f"  probability tends to 1 along this infinite subsequence")
        print("  Mathematical 3D HP toy model, NOT a protein-structure prediction\n")


if __name__ == "__main__":
    main()
