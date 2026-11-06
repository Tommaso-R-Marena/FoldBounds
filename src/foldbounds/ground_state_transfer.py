"""Exact complete all-H square-lattice ground-state degeneracy (side 2..8).

Connected-frontier integer DP counts canonical first-East Hamiltonian walks
over all compatible square placements. Not Lean-verified or real-protein physics.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache


def _canon(front, left):
    labels = {}
    out = []
    for v in (*front, left):
        if v and v not in labels:
            labels[v] = len(labels) + 1
        out.append(labels.get(v, 0))
    return tuple(out[:-1]), out[-1]


def _step(state, x, y, w, h, right, down, degree, force_east):
    front, left = state
    up = front[x]
    if int(bool(up)) + int(bool(left)) + int(right) + int(down) != degree:
        return None
    if (right and x == w - 1) or (down and y == h - 1):
        return None
    if force_east and (left or up or not right or down):
        return None
    if up and up == left:  # would close a cycle
        return None
    old = set(front) - {0}
    if left:
        old.add(left)
    merged = {v for v in (up, left) if v}
    label = max(old, default=0) + 1
    newfront = tuple((label if down else 0) if i == x
                     else (label if v in merged else v)
                     for i, v in enumerate(front))
    newleft = label if right else 0
    active = set(newfront) - {0}
    if newleft:
        active.add(newleft)
    vanished = bool(merged) and label not in active
    if any(v not in active for v in old - merged):
        return None
    if not merged and not right and not down:
        return None
    if (x, y) == (w - 1, h - 1):
        return (newfront, newleft) if vanished and not active else None
    if vanished or not active:
        return None
    return _canon(newfront, newleft)


@lru_cache(maxsize=200)
def count_from_start(w: int, h: int, sx: int, sy: int) -> int:
    """Count all Hamiltonian paths starting (sx,sy), with first step East."""
    if any(type(v) is not int for v in (w, h, sx, sy)) or not (
        2 <= w <= 10 and h >= 2 and 0 <= sx < w - 1 and 0 <= sy < h
    ):
        raise ValueError('invalid rectangle or starting site')
    states = {(tuple([0] * w), 0, False): 1}
    for y in range(h):
        for x in range(w):
            nxt_states = defaultdict(int)
            for (front, left, used_end), count in states.items():
                if (x, y) == (sx, sy):
                    possibilities = ((1, False),)
                elif used_end:
                    possibilities = ((2, False),)
                else:
                    possibilities = ((2, False), (1, True))
                for degree, is_end in possibilities:
                    for right in (False, True):
                        for down in (False, True):
                            nxt = _step((front, left), x, y, w, h, right,
                                        down, degree, (x, y) == (sx, sy))
                            if nxt is not None:
                                f, l = nxt
                                nxt_states[f, l, used_end or is_end] += count
            states = nxt_states
    terminal = (tuple([0] * w), 0, True)
    if set(states) - {terminal}:
        raise AssertionError('open terminal connectivity state')
    return states.get(terminal, 0)


@dataclass(frozen=True)
class GroundStateCertificate:
    side: int
    max_contacts: int
    per_start: tuple[tuple[int, int, int], ...]
    degeneracy: int


def _per_start(side):
    # Mirror through the square's horizontal midline: East is unchanged.
    counts = {}
    for y in range((side + 1) // 2):
        for x in range(side - 1):
            count = count_from_start(side, side, x, y)
            counts[x, y] = count
            counts[x, side - 1 - y] = count
    return tuple((x, y, counts[x, y])
                 for y in range(side) for x in range(side - 1))


@lru_cache(maxsize=8)
def build_ground_state_certificate(side: int) -> GroundStateCertificate:
    if type(side) is not int or not 2 <= side <= 8:
        raise ValueError('square side must be an integer from 2 to 8')
    rows = _per_start(side)
    return GroundStateCertificate(side, (side - 1) ** 2, rows,
                                  sum(row[2] for row in rows))


def verify_ground_state_certificate(cert: GroundStateCertificate) -> bool:
    """Recompute the recurrence; separate brute/sequence oracles in tests."""
    if (not isinstance(cert, GroundStateCertificate)
        or type(cert.side) is not int or not 2 <= cert.side <= 8):
        return False
    if type(cert.max_contacts) is not int or cert.max_contacts != (cert.side - 1) ** 2:
        return False
    if type(cert.degeneracy) is not int or cert.degeneracy <= 0:
        return False
    if not isinstance(cert.per_start, tuple) or any(
        not isinstance(row, tuple) or len(row) != 3 or
        any(type(v) is not int for v in row) for row in cert.per_start
    ):
        return False
    rows = _per_start(cert.side)
    return cert.per_start == rows and cert.degeneracy == sum(v for _, _, v in rows)


@dataclass(frozen=True)
class GroundPopulation:
    side: int
    q: Fraction
    threshold: int
    ground_states: int
    walk_count_cap: int
    probability_lower: Fraction
    correction_lower: Fraction
    correction_upper: Fraction


def certify_ground_population(cert: GroundStateCertificate, q: Fraction,
                              threshold: int, *,
                              correction_lower: Fraction = Fraction(1),
                              correction_upper: Fraction = Fraction(1),
                              block_max: int = 11) -> GroundPopulation:
    """Uniform q'>=q event guarantee, also under statewise weights [l,u].

    P(K>=t) >= l*C*q^K / (l*C*q^K + u*N*q^(t-1)),
    where N is independently checked short-block rooted SAW count cap.
    """
    from .saw_count_cap import optimal_saw_count_cap, verify_saw_count_cap
    if not verify_ground_state_certificate(cert):
        raise ValueError('invalid ground-state certificate')
    if not isinstance(q, Fraction) or q < 1 or type(threshold) is not int or not (1 <= threshold <= cert.max_contacts):
        raise ValueError('invalid q or event threshold')
    if (not isinstance(correction_lower, Fraction) or
        not isinstance(correction_upper, Fraction) or
        not 0 < correction_lower <= correction_upper):
        raise ValueError('invalid statewise correction interval')
    cap = optimal_saw_count_cap(cert.side ** 2 - 2, block_max)
    if not verify_saw_count_cap(cap):
        raise AssertionError('count certificate rejected')
    favorable = correction_lower * cert.degeneracy * q ** cert.max_contacts
    unfavorable = correction_upper * cap.count_upper * q ** (threshold - 1)
    return GroundPopulation(cert.side, q, threshold, cert.degeneracy,
        cap.count_upper, favorable / (favorable + unfavorable),
        correction_lower, correction_upper)


@dataclass(frozen=True)
class DirectedGroundUnion:
    side: int
    q: Fraction
    ground_states: int
    directed_overlap: int
    directed_lower: Fraction
    union_lower: Fraction


def directed_ground_union(cert: GroundStateCertificate, q: Fraction) -> DirectedGroundUnion:
    """Combine all ground states with PR #56 PDSAW mass, subtracting overlap."""
    from .partial_directed import partially_directed_spectrum
    if not verify_ground_state_certificate(cert) or not isinstance(q, Fraction) or q < 1:
        raise ValueError('invalid certificate or q')
    spectrum = dict(partially_directed_spectrum('H' * (cert.side ** 2)))
    overlap = spectrum.get(cert.max_contacts, 0)
    if not 0 <= overlap <= cert.degeneracy:
        raise AssertionError('overlap count exceeds ground degeneracy')
    directed = sum((count * q ** score for score, count in spectrum.items()), Fraction(0))
    return DirectedGroundUnion(cert.side, q, cert.degeneracy, overlap,
        directed, directed + (cert.degeneracy - overlap) * q ** cert.max_contacts)
