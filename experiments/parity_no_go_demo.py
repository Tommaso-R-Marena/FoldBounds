"""Reproduce flat-HP-energy structural-confidence ceiling for parity-pure sequences."""
from decimal import Decimal, localcontext
from fractions import Fraction

from foldbounds.local_rewire_confidence import make_macro_target
from foldbounds.parity_no_go import certify_hp_parity_no_go


def scientific(value: Fraction) -> str:
    with localcontext() as ctx:
        ctx.prec = 18
        return f"{Decimal(value.numerator)/Decimal(value.denominator):.12E}"


def main() -> None:
    for dimensions in ((2, 2, 2), (3, 3, 3), (4, 4, 4)):
        path = make_macro_target(*dimensions)
        n = len(path)
        seq = "HP" * (n // 2)
        cert = certify_hp_parity_no_go(
            seq, path, q=Fraction(1000),
            correction_lower=Fraction(1, 2),
            correction_upper=Fraction(2)
        )
        print(f"Simple-cubic HP parity obstruction, n={n}, "
              f"H={cert.hydrophobic_residues}, P={n-cert.hydrophobic_residues}")
        print("  every full conformation has K_HH=0 at ALL q>0")
        print("  constructive distinct SAWs >= 3**(n-2)")
        print("  exact target microstate confidence upper:",
              scientific(cert.probability_microstate_upper))
        print("  conservative rigid fold-class confidence upper:",
              scientific(cert.probability_rigid_fold_upper))
        print("  adversarial microstate corrections: [1/2,2]")
        print("  lattice-model limitation ONLY; no real protein inference.\n")


if __name__ == "__main__":
    main()
