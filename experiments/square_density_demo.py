"""Reproduce model-relative infinite-family square-lattice population certificates."""
from decimal import Decimal, localcontext
from fractions import Fraction

from foldbounds.square_density import find_first_square_tail


def as_decimal(value: Fraction, precision: int = 30) -> str:
    with localcontext() as ctx:
        ctx.prec = precision
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def main() -> None:
    for q, rho in ((3, Fraction(1, 20)), (4, Fraction(1, 4)),
                   (5, Fraction(1, 3))):
        result = find_first_square_tail(q=q, density=rho,
                                        goal_probability=Fraction(99, 100),
                                        max_side=160)
        assert result and result.success
        print(f"q={q}, event K >= ({rho})*n")
        print(f"  infinite sizes: n=m^2, m={result.side_min}+{result.side_step}*k, k>=0")
        print(f"  first n={result.side_min ** 2}, t={result.event_threshold}, "
              f"witness K={result.witness_contacts}")
        print(f"  independently enumerated c16={result.c16}")
        print("  P(event) >=", as_decimal(result.probability_lower))
        print("  failure odds <=", as_decimal(result.failure_odds_upper))
        print("  next-side odds multiplier <", as_decimal(result.next_side_ratio))
        print("  ratio-of-multipliers <", as_decimal(result.next_ratio_multiplier))
        print("  robust to arbitrary per-state weights in [1/2,2]")
        print("  mathematical 2D all-H model; NOT a real protein prediction\n")


if __name__ == "__main__":
    main()
