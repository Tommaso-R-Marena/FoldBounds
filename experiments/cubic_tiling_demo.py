"""Reproduce strict, adversarially robust 3D HP lattice family bounds.

No external structures are accessed. This is synthetic lattice mathematics,
not validation on natural proteins or evidence of folding accuracy.
"""
from decimal import Decimal, localcontext
from fractions import Fraction

from foldbounds.cubic_tiling import (
    TiledCertificate, certify_majority_halfline,
    minimum_certified_integer_q,
)


def fractional_decimal(x: Fraction, precision: int = 45) -> str:
    with localcontext() as ctx:
        ctx.prec = precision
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def run_case(dims, sequence, threshold, q) -> None:
    cert = TiledCertificate(
        dims, sequence, Fraction(q), threshold,
        correction_lower=Fraction(1, 2),
        correction_upper=Fraction(2)
    )
    result = certify_majority_halfline(cert)
    failure_odds = result.complement_mass_upper / result.event_mass_lower
    single_witness_event = cert.correction_lower * cert.q ** result.family_min_contacts
    single_witness_probability = (
        single_witness_event / (single_witness_event + result.complement_mass_upper)
    )
    print(f"3D HP box {tuple(2*d for d in dims)}, n={result.residues}, "
          f"H={sequence.count('H')}, P={sequence.count('P')}")
    print("Event: nonconsecutive H-H contacts >=", threshold)
    print("Exact verified family count:", result.family_paths)
    print("Guaranteed minimum HH contacts of every family member:",
          result.family_min_contacts)
    print("All complete 3D self-avoiding walks <=:", result.total_walks_upper)
    print("Minimum certified integer contact weight q:",
          minimum_certified_integer_q(cert, max_q=20))
    print("At q =", q)
    print("  rigorously certified event probability >=", fractional_decimal(result.probability_lower))
    print("  rigorous competing/event mass ratio <=", fractional_decimal(failure_odds))
    print("  same-K single-witness probability bound =", fractional_decimal(single_witness_probability))
    print("  for EVERY real q above threshold: strict majority; arbitrary "
          "per-state corrections in [1/2,2].")
    print("  MODEL-RELATIVE RESULT ONLY, not real-protein accuracy.\n")


def main() -> None:
    for dims, period, threshold, q in (
        ((3, 3, 3), 11, 54, 5),
        ((4, 4, 4), 11, 256, 6),
        ((6, 6, 6), 11, 864, 5),
    ):
        n = 8 * dims[0] * dims[1] * dims[2]
        sequence = "".join("P" if i % period == 0 else "H" for i in range(n))
        run_case(dims, sequence, threshold, q)


if __name__ == "__main__":
    main()
