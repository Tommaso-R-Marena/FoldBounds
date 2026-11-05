"""Minimal, deterministic reproduction of the rational q-window guarantee."""
from fractions import Fraction

from foldbounds.lattice import build_frontier
from foldbounds.temperature_transport import (
    certify_temperature_window, transport_certificate,
)


def main() -> None:
    cert = build_frontier("HHHH", Fraction(2), 1, depth=4)
    point = transport_certificate(cert, Fraction(3))
    window = certify_temperature_window((cert,), Fraction(2), Fraction(3))
    print("HHHH, event K>=1, source q=2, target q=3")
    print("Target partition:", point.partition_lower, point.partition_upper)
    print("Target population:", point.probability_lower, point.probability_upper)
    print("Uniform q=[2,3] population:", window.probability_lower,
          window.probability_upper)
    print("Uniform majority:", window.guaranteed_majority)
    assert point.probability_lower == point.probability_upper == Fraction(6, 13)
    assert window.guaranteed_majority == "complement"


if __name__ == "__main__":
    main()
