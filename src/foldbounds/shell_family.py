"""Exact mass of extremal-quadrant shells over a disjoint fixed-depth SAW prefix code.

Tested 2D HP research method, NOT currently Lean verified.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from .lattice import (_check_sequence, _check_q, _validate_path, _possible_hh_pairs,
                      _contacts, _neighbors, ROOT, Path)

_QUADRANTS = ((1, 1), (1, -1), (-1, 1), (-1, -1))
_DIRS = ((1, 0), (-1, 0), (0, 1), (0, -1))


@dataclass(frozen=True)
class ShellLower:
    event: Fraction
    complement: Fraction
    distinct_paths: int
    seed_count: int
    accepted_quadrants: int
    @property
    def partition(self) -> Fraction:
        return self.event + self.complement


def extremal_shell_mass(sequence: str, prefix: Path, q: Fraction, event_at_least: int, *,
                         core_length: int) -> ShellLower:
    """Exactly count a prefix-disjoint union of monotone completions.

    Enumerates valid core paths of exactly core_length residues, extending
    the given prefix. Each core contributes every outward monotone quadrant
    for which its endpoint maximizes the signed x+y potential over core.
    These completions are guaranteed SAW and have no nonconsecutive
    suffix/suffix contacts; cross-contacts can only touch the FIRST new site.
    Quadrant pair overlaps are precisely cardinal rays, subtracted once.
    Prefix cores are disjoint by uniqueness of the core_length path prefix.
    No complete length-n paths are enumerated unless core_length == n.
    """
    _check_sequence(sequence)
    _check_q(q)
    n = len(sequence)
    if type(event_at_least) is not int:
        raise ValueError("event_at_least must be an integer")
    _validate_path(prefix, n)
    if type(core_length) is not int or not len(prefix) <= core_length <= n:
        raise ValueError("core_length must be an integer between prefix and full length")
    # Validate BEFORE looking up a cached computation: Python otherwise
    # treats bool and int as equal cache keys, bypassing strict type checks.
    return _exact_shell_mass(sequence, prefix, q, event_at_least, core_length)


@lru_cache(maxsize=128)
def _exact_shell_mass(sequence: str, prefix: Path, q: Fraction,
                      event_at_least: int, core_length: int) -> ShellLower:
    n = len(sequence)
    pairs = _possible_hh_pairs(sequence)
    event = complement = Fraction(0)
    count = seeds = quadrants = 0
    stack = [prefix]
    while stack:
        p = stack.pop()
        m = len(p)
        if m < core_length:
            stack.extend(_neighbors(p))
            continue
        seeds += 1
        fixed = _contacts(p, pairs)
        r = n-m
        if r == 0:
            value = q**fixed
            if fixed >= event_at_least:
                event += value
            else:
                complement += value
            count += 1
            continue
        x,y=p[-1]
        potentials = {sign: max(sx*px+sy*py for px,py in p)
                      for sign in _QUADRANTS for sx,sy in [sign]}
        eligible = tuple(sign for sign in _QUADRANTS if sign[0]*x+sign[1]*y==potentials[sign])
        quadrants += len(eligible)
        pset = set(p)
        # Unique cardinal direction's first appended site can contact only
        # already placed H positions with chain index <= m-2.
        weights = {}
        for dx,dy in _DIRS:
            new=(x+dx,y+dy)
            if new in pset:
                continue
            cross = sum(sequence[m]=='H' and sequence[i]=='H' and i < m-1
                        and abs(new[0]-px)+abs(new[1]-py)==1
                        for i,(px,py) in enumerate(p))
            weights[(dx,dy)] = (fixed+cross, q**(fixed+cross))
        direction_count={d:0 for d in _DIRS}
        for sx,sy in eligible:
            direction_count[(sx,0)] += 1
            direction_count[(0,sy)] += 1
        # For each quadrant: each first-step direction has 2^(r-1)
        # completions. Intersections between adjacent quadrants: one ray.
        for d,multiplicity in direction_count.items():
            if not multiplicity:
                continue
            if d not in weights:
                raise AssertionError("extremal quadrant points into occupied prefix")
            score,w=weights[d]
            copies=multiplicity*(1 << (r-1))-(1 if multiplicity==2 else 0)
            if copies < 0:
                raise AssertionError("quadrant union underflow")
            value=copies*w
            if score >= event_at_least:
                event+=value
            else:
                complement+=value
            count+=copies
    return ShellLower(event,complement,count,seeds,quadrants)