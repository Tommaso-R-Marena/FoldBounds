"""Exact weighted lower masses from low-turn self-avoiding suffix families.

Every counted completion is explicitly checked for self-avoidance, so its
exact rational Boltzmann weight contributes to a valid lower mass. Paths
with <= 2 changes of suffix direction form a polynomial-cardinality
*restricted* family, not an exact solution for all lattice conformations.

Tested Python research code. No Lean equivalence or algorithmic novelty
claim is made.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .lattice import (DIRECTIONS, Path, _check_q, _check_sequence,
                      _contacts, _possible_hh_pairs, _validate_path)


@dataclass(frozen=True)
class TurnFamilyLower:
    event: Fraction
    complement: Fraction
    distinct_paths: int

    @property
    def partition(self) -> Fraction:
        return self.event + self.complement


def turn_family_mass(sequence: str, prefix: Path, q: Fraction,
                     event_at_least: int, *, max_turns: int = 2) -> TurnFamilyLower:
    """Count all valid extensions with at most max_turns *suffix* direction changes.

    The first suffix direction is free; a change is charged only between
    two consecutive suffix moves. This definition is stable for arbitrary
    canonical self-avoiding prefixes, including a bent prefix.

    For the built-in 0..2-turn cases, each direction string is unique,
    yields at most one embedded completion, and there are at most
      4 * [1 + 3(r-1) + 9*(r-1)*(r-2)/2]
    strings for r >= 1 and max_turns=2, independent of the exponential
    population of all self-avoiding walks. Actual self-avoidance and all
    nonconsecutive H-H contacts are recomputed from the full paths.
    """
    _check_sequence(sequence)
    _check_q(q)
    if type(event_at_least) is not int:
        raise ValueError("event threshold must be integer")
    if type(max_turns) is not int or not 0 <= max_turns <= 2:
        raise ValueError("max_turns must be an integer in [0, 2]")
    _validate_path(prefix, len(sequence))
    n = len(sequence)
    pairs = _possible_hh_pairs(sequence)
    event = complement = Fraction(0)
    count = 0
    # The direction before the prefix/suffix boundary is NOT counted as a turn.
    stack = [(prefix, None, 0)]
    while stack:
        path, last_direction, turns = stack.pop()
        if len(path) == n:
            score = _contacts(path, pairs)
            weight = q ** score
            if score >= event_at_least:
                event += weight
            else:
                complement += weight
            count += 1
            continue
        x, y = path[-1]
        for direction in DIRECTIONS:
            turns_next = turns + int(
                last_direction is not None and last_direction != direction)
            if turns_next > max_turns:
                continue
            dx, dy = direction
            point = (x + dx, y + dy)
            if point in path:
                continue
            stack.append((path + (point,), direction, turns_next))
    return TurnFamilyLower(event, complement, count)
