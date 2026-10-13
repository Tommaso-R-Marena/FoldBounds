"""Exact-rational transfer relaxation of weighted 2D HP lattice completions.

The relaxed continuations avoid all FIXED-prefix sites and never reverse
immediately, but may intersect their OWN suffix. Every genuine self-avoiding
completion belongs to this larger finite set. A DP counts contacts crossing
from suffix H residues to placed-prefix H residues exactly. Contacts entirely
within the suffix are bounded by an integral bipartite b-matching.
This is a model-relative tested checker, NOT a Lean-formalized executable.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from .isoperimetric import prefix_contact_cap
from .lattice import (DIRECTIONS, Path, _check_q, _check_sequence,
                      _contacts, _possible_hh_pairs, _validate_path)
from .matching import _maximum_bipartite_bmatching


@dataclass(frozen=True)
class TransferUpper:
    partition: Fraction
    event: Fraction
    complement: Fraction
    future_future_cap: int
    relaxed_walk_count: int
    coefficient_count: int


def future_future_cap(sequence: str, prefix_len: int) -> int:
    """Bipartite b-matching upper bound for H-H contacts with both ends future."""
    n=len(sequence)
    pairs=tuple((i,j) for i,j in _possible_hh_pairs(sequence) if i>=prefix_len)
    caps={i: (3 if i==n-1 else 2) for i in range(prefix_len,n)
          if sequence[i]=='H'}
    return _maximum_bipartite_bmatching(caps,pairs)


@lru_cache(maxsize=128)
def relaxed_cross_spectrum(sequence: str, path: Path) -> tuple[tuple[int,int], ...]:
    """Integer transfer DP indexed by position, direction, cross-contact score.

    Distinct relaxed walks with identical state accumulate integer
    multiplicities. No float arithmetic or learned numerical approximation.
    """
    n,m=len(sequence),len(path)
    fixed=frozenset(path)
    old_dir=(path[-1][0]-path[-2][0],path[-1][1]-path[-2][1])
    states={(path[-1][0],path[-1][1],old_dir[0],old_dir[1],0):1}
    fixed_h=tuple((i,path[i]) for i in range(m) if sequence[i]=='H')
    for j in range(m,n):
        next_states=defaultdict(int)
        for (x,y,dx_old,dy_old,k),multiplicity in states.items():
            for dx,dy in DIRECTIONS:
                if dx==-dx_old and dy==-dy_old:
                    continue
                xx,yy=x+dx,y+dy
                if (xx,yy) in fixed:
                    continue
                cross=0
                if sequence[j]=='H':
                    cross=sum(j-i>=2 and abs(xx-x0)+abs(yy-y0)==1
                              for i,(x0,y0) in fixed_h)
                next_states[(xx,yy,dx,dy,k+cross)]+=multiplicity
        states=next_states
        if not states:
            break
    histogram=defaultdict(int)
    for (*_,score),multiplicity in states.items():
        histogram[score]+=multiplicity
    return tuple(sorted(histogram.items()))


def nonbacktracking_mass_upper(sequence: str, path: Path, q: Fraction,
                               event_at_least: int) -> TransferUpper:
    """Sound (possibly loose) upper bounds on partition and event masses.

    For an actual extension: k_final = fixed + cross + futurefuture.
    A universal cap based on prefix perimeter bounds k_final, and a
    b-matching caps futurefuture. Conditioning the latter on the *observed
    cross score* is safe for genuine paths and can reduce upper weights.
    Relaxed invalid paths are included with nonnegative weights, so
    overcounting them cannot invalidate the enclosure.
    """
    _check_sequence(sequence)
    _check_q(q)
    _validate_path(path,len(sequence))
    if type(event_at_least) is not int:
        raise ValueError('event_at_least must be an integer')
    n,m=len(sequence),len(path)
    pairs=_possible_hh_pairs(sequence)
    fixed_score=_contacts(path,pairs)
    if m==n:
        weight=q**fixed_score
        return TransferUpper(weight,weight if fixed_score>=event_at_least else Fraction(0),
                             weight if fixed_score<event_at_least else Fraction(0),0,1,1)
    flow_cap=future_future_cap(sequence,m)
    iso_cap=prefix_contact_cap(path,n)
    zhi=ehi=nhi=Fraction(0)
    count=0
    spectrum=relaxed_cross_spectrum(sequence,path)
    for cross,multiplicity in spectrum:
        # Every actual continuation has
        # futurefuture <= min(flow_cap, iso_cap-fixed_score-cross).
        # If a relaxed (invalid) walk violates the perimeter cap, its
        # nonnegative overcount remains harmless.
        cap=min(flow_cap,max(0,iso_cap-fixed_score-cross))
        upper_score=fixed_score+cross+cap
        weight=multiplicity*q**upper_score
        count+=multiplicity
        zhi+=weight
        if upper_score>=event_at_least:
            ehi+=weight
        if fixed_score+cross<event_at_least:
            nhi+=weight
    return TransferUpper(zhi,ehi,nhi,flow_cap,count,len(spectrum))
