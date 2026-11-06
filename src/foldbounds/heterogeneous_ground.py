"""Composition-robust HP event bounds from complete Hamiltonian-square families.

Every P residue removes at most two nonbackbone contacts, or three when
it is a chain endpoint. Reuse exact all-H ground-state *geometrical* counts
to get a guaranteed heterogeneous score floor for each witness.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from .ground_state_transfer import (
    GroundStateCertificate, verify_ground_state_certificate,
)


def compact_hp_score_floor(sequence: str, side: int) -> int:
    if type(side) is not int or not 2 <= side <= 8:
        raise ValueError('side must be an integer in [2,8]')
    if not isinstance(sequence, str) or len(sequence) != side ** 2 or set(sequence) - {'H', 'P'}:
        raise ValueError('expected exactly side^2 H/P residues')
    penalty = sum(3 if i in (0, len(sequence)-1) else 2
                  for i, aa in enumerate(sequence) if aa == 'P')
    return max(0, (side-1) ** 2 - penalty)


@dataclass(frozen=True)
class MixedHPPopulation:
    sequence: str
    side: int
    contact_floor: int
    threshold: int
    q: Fraction
    witness_count: int
    all_walk_count_upper: int
    probability_lower: Fraction
    correction_lower: Fraction
    correction_upper: Fraction


def certify_mixed_hp_population(
    cert: GroundStateCertificate, sequence: str, q: Fraction,
    threshold: int, *, correction_lower: Fraction = Fraction(1),
    correction_upper: Fraction = Fraction(1), block_max: int = 11
) -> MixedHPPopulation:
    """All-q' >= q event lower bound; statewise corrections may vary with q'.

    This does NOT assert that these Hamiltonian square walks are ground
    states of the heterogeneous sequence, only that every one has at least
    the checked score floor.
    """
    from .saw_count_cap import optimal_saw_count_cap, verify_saw_count_cap
    if not verify_ground_state_certificate(cert):
        raise ValueError('ground-state count certificate invalid')
    floor = compact_hp_score_floor(sequence, cert.side)
    if (type(threshold) is not int or not 1 <= threshold <= floor or
        not isinstance(q, Fraction) or q < 1):
        raise ValueError('event not guaranteed for the square family')
    if (not isinstance(correction_lower, Fraction) or
        not isinstance(correction_upper, Fraction) or
        not 0 < correction_lower <= correction_upper):
        raise ValueError('invalid adversarial corrections')
    cap = optimal_saw_count_cap(len(sequence)-2, block_max)
    if not verify_saw_count_cap(cap):
        raise AssertionError('invalid SAW count certificate')
    good = correction_lower * cert.degeneracy * q ** floor
    bad = correction_upper * cap.count_upper * q ** (threshold-1)
    return MixedHPPopulation(sequence, cert.side, floor, threshold, q,
        cert.degeneracy, cap.count_upper, good / (good + bad),
        correction_lower, correction_upper)
