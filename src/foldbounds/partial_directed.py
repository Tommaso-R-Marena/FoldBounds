"""Exact finite-HP contact spectra of partially directed self-avoiding walks.

This is a *lower-mass witness family*, not a model of all 2D self-avoiding
walks. Three pairwise disjoint directional languages are counted:

* E-forward walks: E/N/S moves, with no N<->S reversal unless an E intervenes;
* N-forward walks: N/E/W moves, with no E<->W reversal unless an N intervenes,
  and at least one W;
* S-forward walks: S/E/W moves, with no E<->W reversal unless an S intervenes,
  and at least one W.

Every path is rooted at ((0,0),(1,0)). Each language is self-avoiding
because longitudinal coordinates never decrease and lateral steps within
one longitudinal layer are strictly monotone. The E family has no W;
N and S families with a W have a different vertical sign, hence the
three languages cannot overlap. Prefix restrictions are imposed as
actual step constraints, so each returned path extends its given prefix.

The finite-state DP stores the current layer's signed lateral displacement,
the previous layer's starting lateral displacement relative to the current
start, and the current layer's first sequence index. Because a layer is a
strictly monotone interval, the preceding layer's residue at a given
lateral coordinate, if any, is uniquely reconstructed. All nonconsecutive
contacts in these paths lie across *adjacent* longitudinal layers; contacts
from a new lateral step into the previous layer are scored exactly.

This is a Python mathematical research prototype, not a Lean-verified
implementation or a claim of novel asymptotic self-avoiding-walk theory.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from typing import Literal

from .lattice import Path, _check_sequence, _validate_path

# (longitudinal forward, lateral plus, lateral minus)
_AXES: dict[str, tuple[tuple[int, int], tuple[int, int], tuple[int, int]]] = {
    "E": ((1, 0), (0, 1), (0, -1)),
    "N": ((0, 1), (1, 0), (-1, 0)),
    "S": ((0, -1), (1, 0), (-1, 0)),
}


def _rooted_axis_spectrum(sequence: str, prefix: Path, axis: Literal["E", "N", "S"]) -> dict[int, int]:
    """Histogram in one orientation, restricted by an arbitrary legal prefix.

    A state is (current_offset, previous_start_offset, current_start_index,
    has_west, score). previous_start_offset is None only in the initial
    N/S layer (the initial E bond lies entirely in the same row).

    State counts are nonnegative integers. Each step chooses exactly one
    allowed direction and preserves its corresponding physical walk.
    """
    n = len(sequence)
    forward, side_plus, side_minus = _AXES[axis]
    # The root's first bond lies across columns in E orientation, but within
    # the single starting row in N/S orientations.
    if axis == "E":
        start = (0, 0, 1, False, 0)
    else:
        start = (1, None, 0, False, 0)
    states: dict[tuple[int, int | None, int, bool, int], int] = {start: 1}
    for j in range(2, n):
        constraint = ((prefix[j][0] - prefix[j-1][0],
                       prefix[j][1] - prefix[j-1][1]) if j < len(prefix) else None)
        nxt: dict[tuple[int, int | None, int, bool, int], int] = defaultdict(int)
        for (dc, dp, cs, has_west, k), count in states.items():
            if constraint is None or constraint == forward:
                nxt[(0, -dc, j, has_west, k)] += count
            for sign, delta in ((1, side_plus), (-1, side_minus)):
                if constraint is not None and constraint != delta:
                    continue
                if dc * sign < 0:
                    continue  # would revisit the site immediately behind
                new_dc = dc + sign
                new_score = k
                if dp is not None and min(dp, 0) <= new_dc <= max(dp, 0):
                    prev_start = cs - 1 - abs(dp)
                    i = prev_start + abs(new_dc - dp)
                    # Independently impose the non-backbone pair condition.
                    if j - i >= 2 and sequence[i] == sequence[j] == "H":
                        new_score += 1
                nxt[(new_dc, dp, cs, has_west or (axis != "E" and sign == -1),
                     new_score)] += count
        states = nxt
        if not states:
            break
    out: dict[int, int] = defaultdict(int)
    for (_, _, _, seen_west, k), count in states.items():
        if axis == "E" or seen_west:
            out[k] += count
    return dict(out)


@lru_cache(maxsize=64)
def partially_directed_spectrum(sequence: str, prefix: Path = ((0, 0), (1, 0))) -> tuple[tuple[int, int], ...]:
    """Exact (contact count, multiplicity) of a disjoint genuine SAW subset.

    Calling this for every bounded subtree in a verified prefix frontier
    yields disjoint lower mass; a complete-conformation witness matching
    this language must be removed before crediting the spectrum additionally.
    """
    _check_sequence(sequence)
    _validate_path(prefix, len(sequence))
    out: dict[int, int] = defaultdict(int)
    for axis in ("E", "N", "S"):
        for k, count in _rooted_axis_spectrum(sequence, prefix, axis).items():
            out[k] += count
    return tuple(sorted(out.items()))


def is_partially_directed_path(path: Path) -> bool:
    """Independently classify a full anchored walk (geometry validated by caller).

    Uses step-language acceptance rather than DP score computation to allow
    overlap checks against externally supplied full-length witnesses.
    """
    if len(path) < 2 or path[:2] != ((0, 0), (1, 0)):
        return False
    moves = tuple((b[0] - a[0], b[1] - a[1]) for a, b in zip(path[1:], path[2:]))
    for axis, (forward, plus, minus) in _AXES.items():
        # N/S initial horizontal layer already moved +1 from root index 0.
        lateral = 0 if axis == "E" else 1
        west = False
        valid = True
        for move in moves:
            if move == forward:
                lateral = 0
            elif move == plus and lateral >= 0:
                lateral += 1
            elif move == minus and lateral <= 0:
                lateral -= 1
                if axis != "E":
                    west = True
            else:
                valid = False
                break
        if valid and (axis == "E" or west):
            return True
    return False


def east_directed_path_count(remaining_steps: int) -> int:
    """Closed recurrence for rooted E-progressive SAW subset cardinality.

    Number of E/N/S paths with vertical sign fixed inside each column:
    A(z)=(1+z)/(1-2z-z²), a(0)=1, a(1)=3 and
    a(r)=2*a(r-1)+a(r-2). Runs in O(r) integer operations, with no DP
    over coordinates or contacts. This is classical PDSAW enumeration.
    """
    if type(remaining_steps) is not int or remaining_steps < 0:
        raise ValueError('remaining_steps must be a nonnegative integer')
    if remaining_steps == 0:
        return 1
    previous, current = 1, 3
    for _ in range(2, remaining_steps+1):
        previous, current = current, 2*current+previous
    return current


def universal_hp_partition_lower(sequence: str, q) -> int:
    """Nonenumerative rigorous Z>=a(n-2) for every 2D HP sequence and q>=1.

    For every E-progressive rooted walk, its Boltzmann weight q**k >=1.
    Distinct directed walks are genuine distinct anchored SAWs, therefore
    their elementary integer cardinality lower-bounds the full partition.
    This coarse bound loses all energy information but scales to long n.
    """
    from .lattice import _check_q
    _check_sequence(sequence)
    _check_q(q)
    return east_directed_path_count(len(sequence)-2)
