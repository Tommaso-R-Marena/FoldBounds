"""Reproduce heterogeneous 3D infinite-size certificates at 20% and 25% P."""
from fractions import Fraction
from foldbounds.cubic_general_tail import certify_general_cubic_density_tail


def main() -> None:
    for side, period, density, q in (
        (20, 5, Fraction(1, 4), 5),
        (24, 4, Fraction(1, 4), 8),
    ):
        cert = certify_general_cubic_density_tail(
            side, polar_period=period, density=density, q=q,
            correction_lower=Fraction(1, 2),
            correction_upper=Fraction(2),
            goal_probability=Fraction(99, 100),
        )
        print(f"3D H/P lattice, 1 P every {period} sites, q={q}")
        print(f"  Infinite progression: macro side {side}+{cert.side_step}*k (k>=0)")
        print(f"  First full box: side {2*side}, n={cert.residues}")
        print(f"  polar residues={cert.polar_residues}")
        print(f"  event: K_HH>={cert.event_threshold}, minimum family HH={cert.hh_contacts_lower}")
        print(f"  verified family count >= 2*4**({side**3})")
        print(f"  model-relative failure odds < 10^-9: "
              f"{cert.failure_odds_upper < Fraction(1, 10**9)}")
        print(f"  P(event)>99% for ALL eligible sizes: "
              f"{cert.probability_lower > Fraction(99, 100)}")
        print("  exact arithmetic; statewise corrections in [1/2,2]; "
              "no experimental protein accuracy claims\n")


if __name__ == "__main__":
    main()
