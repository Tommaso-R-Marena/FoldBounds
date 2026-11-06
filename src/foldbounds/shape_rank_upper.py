"""Exact-rational 2D HP upper mass via shape-conditioned relaxed walks.

A finite-memory walk transfer strictly contains the true SAW ensemble. Its
endpoint/bounding-box counts yield a per-relaxed-path isoperimetric maximum
contact score. A separate, verified short-block count cap bounds the total
number of *true* SAW completions, so sorting candidate cap weights and
keeping at most that many is a rigorous uniform-matroid upper bound.

Every inclusion is model-relative; this Python implementation is not a Lean
proof, and transfer memory ignores arbitrary old self-intersections.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from .horizon_certificate import ROOT, Path, _validate, _score_interval
from .saw_count_cap import optimal_saw_count_cap, verify_saw_count_cap

STEPS=((1,0),(0,1),(-1,0),(0,-1))

@dataclass(frozen=True)
class ShapeClass:
    score_upper: int
    relaxed_count: int

@dataclass(frozen=True)
class ShapeRankCertificate:
    sequence: str
    prefix: Path
    suffix_memory: int
    histogram: tuple[ShapeClass,...]
    true_count_upper: int

@dataclass(frozen=True)
class ShapeRankUpper:
    partition_upper: Fraction
    event_upper: Fraction
    complement_upper: Fraction
    relaxed_paths: int
    ranked_paths: int
    count_upper: int
    score_classes: int


def _prefix_directions(path: Path, history_length: int) -> tuple[int,...]:
    return tuple(STEPS.index((b[0]-a[0],b[1]-a[1])) for a,b in zip(path,path[1:]))[-history_length:]


def _contains_recent(new: tuple[int,int], endpoint: tuple[int,int], history: tuple[int,...]) -> bool:
    if new==endpoint:return True
    x,y=endpoint
    for direction in reversed(history):
        dx,dy=STEPS[direction]
        x-=dx;y-=dy
        if new==(x,y):return True
    return False


def _relaxed_shape_histogram(sequence: str, prefix: Path, suffix_memory: int,
                              max_states: int) -> tuple[ShapeClass,...]:
    """Count a superset of all genuine completions by (contact-score cap).

    States track last memory-1 directed steps, absolute tip and extremal
    coordinates, not the full visited set. Fixed prefix sites are obstacles.
    Old suffix sites outside memory may be revisited, making this a safe
    *overcount*. Final shape area smaller than n rules out a genuine SAW.
    """
    _validate(sequence,prefix)
    if type(suffix_memory) is not int or not 2<=suffix_memory<=7:
        raise ValueError('suffix_memory must lie in [2,7]')
    if type(max_states) is not int or max_states<1:
        raise ValueError('max_states must be positive')
    n=len(sequence);blocked=set(prefix)
    xs=[p[0] for p in prefix];ys=[p[1] for p in prefix]
    initial=(*prefix[-1],min(xs),max(xs),min(ys),max(ys),
             _prefix_directions(prefix,suffix_memory-1))
    states={initial:1}
    for _ in range(len(prefix),n):
        nxt=defaultdict(int)
        for (x,y,xmin,xmax,ymin,ymax,history),count in states.items():
            for i,(dx,dy) in enumerate(STEPS):
                p=(x+dx,y+dy)
                if p in blocked or _contains_recent(p,(x,y),history):continue
                entry=(p[0],p[1],min(xmin,p[0]),max(xmax,p[0]),
                       min(ymin,p[1]),max(ymax,p[1]),
                       (history+(i,))[-(suffix_memory-1):])
                nxt[entry]+=count
        states=nxt
        if len(states)>max_states:
            raise OverflowError('relaxed transfer exceeds declared finite-state budget')
        if not states:break
    low, high = _score_interval(sequence,prefix)
    histogram=defaultdict(int)
    for (_,_,xmin,xmax,ymin,ymax,_),count in states.items():
        width=xmax-xmin+1;height=ymax-ymin+1
        if width*height<n:continue  # impossible n distinct lattice vertices
        cap=min(high,max(0,n+1-width-height))
        if cap<low:continue  # cannot retain already established contacts
        histogram[cap]+=count
    return tuple(ShapeClass(k,histogram[k]) for k in sorted(histogram,reverse=True))


def build_shape_rank_certificate(sequence: str, prefix: Path=ROOT, *,
                                 suffix_memory: int=4, max_states: int=1500000) -> ShapeRankCertificate:
    histogram=_relaxed_shape_histogram(sequence,prefix,suffix_memory,max_states)
    count_cert=optimal_saw_count_cap(len(sequence)-len(prefix))
    if not verify_saw_count_cap(count_cert):raise AssertionError('block count invalid')
    return ShapeRankCertificate(sequence,prefix,suffix_memory,histogram,count_cert.count_upper)


def verify_shape_rank_certificate(cert: ShapeRankCertificate, q: Fraction,
                                  event_at_least: int=1,
                                  max_states: int=1500000) -> ShapeRankUpper:
    """Independently replay shape histogram and check greedy largest-C upper.

    The histogram generator and q-weight optimizer are untrusted. Here the
    complete finite-memory transfer is re-executed and each class count/cap
    compared. For a scalable small proof artifact, a later refinement should
    verify per-step transition histograms instead of recomputing them.
    """
    if not isinstance(cert,ShapeRankCertificate):raise ValueError('malformed shape certificate')
    _validate(cert.sequence,cert.prefix)
    if not isinstance(q,Fraction) or q<1 or type(event_at_least) is not int:
        raise ValueError('expected Fraction q>=1 and integer threshold')
    recount=_relaxed_shape_histogram(cert.sequence,cert.prefix,cert.suffix_memory,max_states)
    if cert.histogram!=recount:raise ValueError('mismatched geometric transfer histogram')
    count_cert=optimal_saw_count_cap(len(cert.sequence)-len(cert.prefix))
    if not verify_saw_count_cap(count_cert) or cert.true_count_upper!=count_cert.count_upper:
        raise ValueError('forged SAW total count cap')
    available=count_cert.count_upper
    z=e=c=Fraction(0)
    chosen=0
    for cls in cert.histogram:
        count=min(available,cls.relaxed_count)
        if not count:break
        mass=count*q**cls.score_upper
        z+=mass
        if cls.score_upper>=event_at_least:e+=mass
        c+=mass  # Without a positive lower score, any candidate may be complement.
        available-=count;chosen+=count
    # Complement upper can be tightened to zero if the prefix itself already
    # has >= threshold contacts, independently of future geometry.
    fixed_score,_=_score_interval(cert.sequence,cert.prefix)
    if fixed_score>=event_at_least:c=Fraction(0)
    return ShapeRankUpper(z,e,c,sum(x.relaxed_count for x in cert.histogram),
                          chosen,count_cert.count_upper,len(cert.histogram))
