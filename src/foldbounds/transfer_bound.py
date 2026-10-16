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

from .isoperimetric import prefix_contact_cap, minimum_final_box_sum
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
    suffix_memory: int = 2


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



@lru_cache(maxsize=128)
def recent_memory_contact_spectrum(
        sequence: str, path: Path, suffix_memory: int,
        local_cap: int | None = None,
        with_endpoint: bool = False
) -> tuple[tuple[tuple, int], ...]:
    """Exact relaxed-walk histogram by (prefix contacts, local suffix contacts).

    A relaxed continuation avoids the fixed prefix and the most recent K
    sites. For each suffix residue j, any contact (i,j) with i >= prefix_len,
    j-i >= 2 and j-i <= K is visible while j is appended. Such pairs are
    counted once, even if a longer relaxed trajectory is not self-avoiding.
    For a genuine conformation these observed pairs are a *subset* of its
    future/future H-H contacts, so they give a safe contact-score floor.

    Distinct paths with equal history and both scores are aggregated with
    integer multiplicities. With local_cap, the second score is the exact
    *saturated* value min(observed_local, local_cap). This safely preserves
    threshold classification while limiting state growth; None gives the
    full exact local score. This is a relaxed-walk histogram, not a true
    SAW ensemble histogram. With with_endpoint=True, the histogram key is
    (cross, local, terminal_site), retaining geometry without path enumeration.
    """
    _check_sequence(sequence)
    _validate_path(path, len(sequence))
    if type(suffix_memory) is not int or not 2 <= suffix_memory <= 8:
        raise ValueError("suffix_memory must be an integer in [2, 8]")
    if local_cap is not None and (type(local_cap) is not int or local_cap < 0):
        raise ValueError("local_cap must be a nonnegative integer or None")
    if type(with_endpoint) is not bool:
        raise ValueError("with_endpoint must be a boolean")
    if suffix_memory == 2 and not with_endpoint:
        return tuple(((cross, 0), multiplicity)
                     for cross, multiplicity in relaxed_cross_spectrum(sequence, path))

    n, m = len(sequence), len(path)
    fixed = frozenset(path)
    fixed_h = tuple((i, path[i]) for i in range(m) if sequence[i] == 'H')
    # State: recent sites, exact prefix/suffix contacts, observed suffix contacts.
    states = {(path[-1:], 0, 0): 1}
    for j in range(m, n):
        next_states = defaultdict(int)
        for (history, cross, local), multiplicity in states.items():
            x, y = history[-1]
            first_history_index = j - len(history)
            for dx, dy in DIRECTIONS:
                point = (x + dx, y + dy)
                if point in fixed or point in history:
                    continue
                cross_new = local_new = 0
                if sequence[j] == 'H':
                    cross_new = sum(
                        j - i >= 2 and abs(point[0] - px) + abs(point[1] - py) == 1
                        for i, (px, py) in fixed_h
                    )
                    if local_cap is None or local < local_cap:
                        local_new = sum(
                            i >= m and j - i >= 2 and sequence[i] == 'H'
                            and abs(point[0] - px) + abs(point[1] - py) == 1
                            for i, (px, py) in enumerate(history, first_history_index)
                        )
                next_history = (history + (point,))[-suffix_memory:]
                new_local = local + local_new
                if local_cap is not None:
                    new_local = min(new_local, local_cap)
                next_states[(next_history, cross + cross_new,
                             new_local)] += multiplicity
        states = next_states
        if not states:
            break

    histogram = defaultdict(int)
    for (history, cross, local), multiplicity in states.items():
        key = (cross, local, history[-1]) if with_endpoint else (cross, local)
        histogram[key] += multiplicity
    return tuple(sorted(histogram.items()))


@lru_cache(maxsize=128)
def recent_memory_cross_spectrum(
        sequence: str, path: Path, suffix_memory: int
) -> tuple[tuple[int, int], ...]:
    """Backward-compatible marginal over local suffix contact scores."""
    marginal = defaultdict(int)
    for (cross, _), multiplicity in recent_memory_contact_spectrum(
            sequence, path, suffix_memory, local_cap=0):
        marginal[cross] += multiplicity
    return tuple(sorted(marginal.items()))


def nonbacktracking_mass_upper(sequence: str, path: Path, q: Fraction,
                               event_at_least: int,
                               suffix_memory: int = 2, *,
                               terminal_box: bool = False) -> TransferUpper:
    """Sound (possibly loose) upper bounds on partition and event masses.

    For an actual extension: k_final = fixed + cross + futurefuture.
    A universal cap based on prefix perimeter bounds k_final, and a
    b-matching caps futurefuture. Conditioning the latter on the *observed
    cross score* is safe for genuine paths and can reduce upper weights.
    Relaxed invalid paths are included with nonnegative weights, so
    overcounting them cannot invalidate the enclosure.

    With terminal_box, the transfer DP preserves each relaxed path's final
    endpoint. For a genuine completion, the enclosing bounding box must
    contain both that endpoint and every fixed-prefix site. The square-grid
    perimeter/isoperimetric contact cap therefore strengthens conditionally
    on the known terminal point. This caps both event and total upper mass,
    without assuming anything about invalid relaxed trajectories.
    """
    _check_sequence(sequence)
    _check_q(q)
    _validate_path(path,len(sequence))
    if type(event_at_least) is not int:
        raise ValueError('event_at_least must be an integer')
    if type(suffix_memory) is not int or not 2 <= suffix_memory <= 8:
        raise ValueError('suffix_memory must be an integer in [2, 8]')
    if type(terminal_box) is not bool:
        raise ValueError('terminal_box must be a boolean')
    n,m=len(sequence),len(path)
    pairs=_possible_hh_pairs(sequence)
    fixed_score=_contacts(path,pairs)
    if m==n:
        weight=q**fixed_score
        return TransferUpper(weight,weight if fixed_score>=event_at_least else Fraction(0),
                             weight if fixed_score<event_at_least else Fraction(0),0,1,1,
                             suffix_memory)
    flow_cap=future_future_cap(sequence,m)
    iso_cap=prefix_contact_cap(path,n)
    zhi=ehi=nhi=Fraction(0)
    count=0
    # Threshold saturation retains exactly the local-contact information
    # needed for event/complement classification, avoiding excess DP states.
    local_cap=max(0,event_at_least-fixed_score)
    spectrum=recent_memory_contact_spectrum(
        sequence,path,suffix_memory,local_cap=local_cap,
        with_endpoint=terminal_box)
    xs = [site[0] for site in path]
    ys = [site[1] for site in path]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)

    @lru_cache(maxsize=None)
    def terminal_contact_cap(point: tuple[int, int]) -> int:
        x,y = point
        width = max(xmax, x)-min(xmin, x)+1
        height = max(ymax, y)-min(ymin, y)+1
        return max(0,n+1-minimum_final_box_sum(n,width,height))

    for key,multiplicity in spectrum:
        if terminal_box:
            cross,local_contacts,terminal=key
            geometry_cap=min(iso_cap,terminal_contact_cap(terminal))
        else:
            cross,local_contacts=key
            geometry_cap=iso_cap
        # Every actual continuation has
        # futurefuture <= min(flow_cap, iso_cap-fixed_score-cross).
        # If a relaxed (invalid) walk violates the perimeter cap, its
        # nonnegative overcount remains harmless.
        cap=min(flow_cap,max(0,geometry_cap-fixed_score-cross))
        upper_score=fixed_score+cross+cap
        weight=multiplicity*q**upper_score
        count+=multiplicity
        zhi+=weight
        if upper_score>=event_at_least:
            ehi+=weight
        # Observed nonconsecutive suffix contacts are guaranteed for a
        # genuine walk, even though the relaxed DP does not remember every
        # older site. Prune impossible complement classifications only.
        if fixed_score+cross+local_contacts<event_at_least:
            nhi+=weight
    return TransferUpper(zhi,ehi,nhi,flow_cap,count,len(spectrum),suffix_memory)
