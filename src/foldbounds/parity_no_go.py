"""Global HP-model parity obstruction and nonidentifiability certificate.

For any canonical cubic-lattice self-avoiding walk, each unit step flips
the parity of x+y+z. Therefore a noncovalent nearest-neighbor HH contact
requires hydrophobic residues on OPPOSITE sequence-index parities.

If all H positions have the SAME parity (or none are H), the contact count
is exactly ZERO for *every* conformation, independent of the chain
geometry. For the classical HP Hamiltonian with no HP/PP interactions,
all conformations consequently have nominal Boltzmann weight one at
every contact fugacity q>0.

After fixing root sites (0,0,0),(1,0,0), every word of n-2 steps from
{+x,+y,+z} yields a distinct, valid, coordinate-monotone self-avoiding
walk. Thus the conformational ensemble has at least 3**(n-2) states.
This constructive lower bound supplies exact confidence ceilings
without generating the exponentially many conformations.

The parity obstruction is standard bipartite-lattice mathematics, not
a novel discovery, not a statement about all polymer models, and not
a valid conclusion for real protein contact-force fields.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .local_rewire_confidence import StructuralModel, _validate, hp_pair_weights

Point = tuple[int, int, int]
Path = tuple[Point, ...]


@dataclass(frozen=True)
class ParityNoGoResult:
    residues: int
    hydrophobic_residues: int
    compatible_h_index_parity: int | None
    all_conformations_have_zero_hh_contacts: bool
    monotone_walks_lower: int
    probability_microstate_upper: Fraction
    probability_rigid_fold_upper: Fraction
    correction_lower: Fraction
    correction_upper: Fraction

    @property
    def rules_out_exact_state_majority(self) -> bool:
        return self.probability_microstate_upper < Fraction(1, 2)


def certify_hp_parity_no_go(
    sequence: str, target: Path, *, q: Fraction,
    correction_lower: Fraction = Fraction(1),
    correction_upper: Fraction = Fraction(1),
) -> ParityNoGoResult:
    """Either reject insufficient hypotheses or certify the ENTIRE ensemble.

    The target is validated as a full canonical 3D walk, but no special
    geometry of the target is required for the parity obstruction.

    Validity scope: exactly the simple-cubic HP model with positive HH
    contact weight q, HP and PP weight 1, fixed chain length and
    fully self-avoiding rooted 3D walks, with arbitrary positive
    statewise corrections in [lower,upper].
    """
    _validate(StructuralModel(
        sequence, target, hp_pair_weights(q),
        correction_lower, correction_upper
    ))
    h_parities = {i % 2 for i, residue in enumerate(sequence) if residue == "H"}
    if not set(sequence) <= {"H", "P"}:
        raise ValueError("only the explicitly declared H/P alphabet is supported")
    if len(h_parities) > 1:
        raise ValueError("inconclusive: H residues occupy both chain-index parities")
    compatible = next(iter(h_parities)) if h_parities else None
    n = len(sequence)
    # Every continuation of the fixed first bond using positive coordinate
    # steps +x,+y,+z strictly increases x+y+z. Consequently no revisit is
    # possible, and distinct words produce distinct indexed walks.
    # The 3^(n-2) lower count is independent of the proposed target.
    universal_count_lower = 3**(n - 2)
    lo, hi = correction_lower, correction_upper
    exact_cap = hi / (hi + lo * (universal_count_lower - 1))
    if universal_count_lower <= 48:
        orbit_cap = Fraction(1)
    else:
        orbit_cap = (
            48 * hi
            / (48 * hi + lo * (universal_count_lower - 48))
        )
    return ParityNoGoResult(
        n, sum(residue == "H" for residue in sequence),
        compatible, True, universal_count_lower, exact_cap,
        orbit_cap, lo, hi
    )
