"""Proof-carrying finite-horizon SAW *upper* envelopes for 2D HP chains.

Each genuine continuation belongs to exactly one certified leaf. Unknown
leaves use independently checkable rooted-SA W block-count and prefix-box
contact caps; complete leaves contribute exact Boltzmann weights. A finite
vacancy-component test can certify nonextendable prefixes even when the
terminal site has free neighbors (insufficient area or parity capacity).

Unlike restricted-language witness methods, this upper bound controls *all*
conformations. This is elementary exact-integer Python assurance, NOT Lean
kernel verification and NOT a new protein-physics theorem.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Literal

from .saw_count_cap import optimal_saw_count_cap, verify_saw_count_cap

Point=tuple[int,int]
Path=tuple[Point,...]
ROOT: Path=((0,0),(1,0))
STEPS=((1,0),(-1,0),(0,1),(0,-1))

@dataclass(frozen=True)
class HorizonLeaf:
    path: Path
    kind: Literal['exact','dead','bound']
    score_lower: int
    score_upper: int
    count_upper: int

@dataclass(frozen=True)
class HorizonCertificate:
    sequence: str
    prefix: Path
    horizon: int
    leaves: tuple[HorizonLeaf,...]

@dataclass(frozen=True)
class HorizonEnvelope:
    partition_upper: Fraction
    event_upper: Fraction
    complement_upper: Fraction
    exact_mass: Fraction
    exact_event: Fraction
    exact_complement: Fraction
    leaves: int
    dead_leaves: int
    bounded_leaves: int


def _validate(seq: str, path: Path) -> None:
    if not isinstance(seq,str) or len(seq)<2 or set(seq)-{'H','P'}:
        raise ValueError('expected a 2D HP sequence with at least two residues')
    if not isinstance(path,tuple) or not 2<=len(path)<=len(seq) or path[:2]!=ROOT:
        raise ValueError('invalid anchored prefix')
    if any(not isinstance(p,tuple) or len(p)!=2 or any(type(z) is not int for z in p) for p in path):
        raise ValueError('integer coordinates required')
    if len(set(path))!=len(path):
        raise ValueError('prefix intersects itself')
    if any(abs(a[0]-b[0])+abs(a[1]-b[1])!=1 for a,b in zip(path,path[1:])):
        raise ValueError('prefix must take lattice steps')


def _children(path: Path) -> tuple[Path,...]:
    used=set(path); x,y=path[-1]
    return tuple(path+((x+dx,y+dy),) for dx,dy in STEPS if (x+dx,y+dy) not in used)


def _score_interval(seq: str, path: Path) -> tuple[int,int]:
    m=len(path); known=pending=0
    for i in range(len(seq)):
        if seq[i]!='H':continue
        for j in range(i+3,len(seq),2):
            if seq[j]!='H':continue
            if j>=m:pending+=1
            elif abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1:known+=1
    xs=[p[0] for p in path];ys=[p[1] for p in path]
    w0=max(xs)-min(xs)+1; h0=max(ys)-min(ys)+1;n=len(seq)
    box=min(w+max(h0,(n+w-1)//w) for w in range(w0,n+1))
    return known, min(known+pending,max(0,n+1-box))


def _component_possible(path: Path, n: int) -> bool:
    """Necessary vacancy parity capacity. False rigorously proves no suffix.

    Each extension after the current tip lies inside one component of unplaced
    grid sites adjacent to the tip. An r-vertex suffix uses ceil(r/2)
    opposite-checkerboard sites and floor(r/2) same-checkerboard sites.
    A component exceeding the finite search budget is conservatively allowed.
    """
    r=n-len(path)
    if r==0:return True
    blocked=set(path); tip=path[-1]
    starts=[(tip[0]+dx,tip[1]+dy) for dx,dy in STEPS
            if (tip[0]+dx,tip[1]+dy) not in blocked]
    if not starts:return False
    seen:set[Point]=set()
    p=(tip[0]+tip[1])&1
    for first in starts:
        if first in seen:continue
        queue=[first];seen.add(first);j=0
        while j<len(queue):
            x,y=queue[j];j+=1
            if len(queue)>2*r+2: return True  # incomplete exploration => fail open
            for dx,dy in STEPS:
                v=(x+dx,y+dy)
                if v not in blocked and v not in seen:
                    seen.add(v);queue.append(v)
        other=sum(((x+y)&1)!=p for x,y in queue)
        same=len(queue)-other
        if other>=(r+1)//2 and same>=r//2:
            return True
    return False


def _count_cap(remaining: int) -> int:
    cert=optimal_saw_count_cap(remaining)
    if not verify_saw_count_cap(cert):
        raise AssertionError('SAW count-cap verification failed')
    return cert.count_upper


def build_horizon_certificate(sequence: str, prefix: Path=ROOT, *,
                              horizon: int=4, max_leaves: int=50000) -> HorizonCertificate:
    _validate(sequence,prefix)
    if type(horizon) is not int or horizon<0 or type(max_leaves) is not int or max_leaves<1:
        raise ValueError('horizon>=0 and max_leaves>=1 required')
    target=min(len(sequence),len(prefix)+horizon)
    result: list[HorizonLeaf]=[]
    def visit(path: Path) -> None:
        low,high=_score_interval(sequence,path)
        if len(path)==len(sequence):
            result.append(HorizonLeaf(path,'exact',low,low,1));return
        if not _component_possible(path,len(sequence)):
            result.append(HorizonLeaf(path,'dead',low,high,0));return
        # Budget is fail-open: stop by bounding instead of dropping subtrees.
        if len(path)>=target or len(result)>=max_leaves:
            result.append(HorizonLeaf(path,'bound',low,high,_count_cap(len(sequence)-len(path))))
            return
        for child in _children(path):visit(child)
    visit(prefix)
    return HorizonCertificate(sequence,prefix,horizon,tuple(result))


def verify_horizon_certificate(cert: HorizonCertificate, q: Fraction,
                               event_at_least: int=1) -> HorizonEnvelope:
    """Independently replay complete leaf partition and all safety caps.

    The producer's leaf ordering, count estimates, score bounds and dead-end
    assertions are untrusted. Every leaf must be accessible by an actual
    lattice prefix, and no legal immediate branch may be omitted.
    """
    if not isinstance(cert,HorizonCertificate):raise ValueError('malformed certificate')
    _validate(cert.sequence,cert.prefix)
    if type(cert.horizon) is not int or cert.horizon<0:
        raise ValueError('bad horizon')
    if not isinstance(q,Fraction) or q<1 or type(event_at_least) is not int:
        raise ValueError('exact Fraction q>=1 and integer threshold required')
    if not isinstance(cert.leaves,tuple) or not cert.leaves:
        raise ValueError('empty leaf frontier')
    submitted:dict[Path,HorizonLeaf]={}
    for node in cert.leaves:
        if not isinstance(node,HorizonLeaf):raise ValueError('malformed leaf')
        _validate(cert.sequence,node.path)
        if node.path[:len(cert.prefix)]!=cert.prefix or node.path in submitted:
            raise ValueError('duplicate or alien prefix')
        if len(node.path)>len(cert.prefix)+cert.horizon:
            raise ValueError('leaf beyond horizon')
        low,high=_score_interval(cert.sequence,node.path)
        if type(node.score_lower) is not int or node.score_lower>low or type(node.score_upper) is not int or node.score_upper<high:
            raise ValueError('forged score enclosure')
        if type(node.count_upper) is not int or node.count_upper<0:
            raise ValueError('bad count cap')
        if node.kind=='exact':
            if (len(node.path)!=len(cert.sequence) or node.count_upper!=1 or
                    node.score_lower!=low or node.score_upper!=low):
                raise ValueError('forged exact node')
        elif node.kind=='dead':
            if len(node.path)==len(cert.sequence) or _component_possible(node.path,len(cert.sequence)) or node.count_upper!=0:
                raise ValueError('unsupported dead-end claim')
        elif node.kind=='bound':
            if len(node.path)==len(cert.sequence) or node.count_upper<_count_cap(len(cert.sequence)-len(node.path)):
                raise ValueError('forged SAW count cap')
        else:raise ValueError('unknown node kind')
        submitted[node.path]=node
    entered:set[Path]=set(); z=e=n=exact=ee=nn=Fraction(0)
    dead=bounded=0
    def visit(path: Path) -> None:
        nonlocal z,e,n,exact,ee,nn,dead,bounded
        node=submitted.get(path)
        if node is not None:
            entered.add(path)
            if node.kind=='dead':dead+=1;return
            if node.kind=='bound':bounded+=1
            mass=node.count_upper*q**node.score_upper
            z+=mass
            if node.score_upper>=event_at_least:e+=mass
            if node.score_lower<event_at_least:n+=mass
            if node.kind=='exact':
                exact+=mass
                if node.score_lower>=event_at_least:ee+=mass
                else:nn+=mass
            return
        if len(path)==len(cert.sequence):raise ValueError('terminal conformation omitted')
        for child in _children(path):
            if child not in prefixes:raise ValueError('missing legal child subtree')
            visit(child)
    prefixes={path[:k] for path in submitted for k in range(len(cert.prefix),len(path)+1)}
    if cert.prefix not in prefixes:raise ValueError('root missing')
    visit(cert.prefix)
    if entered!=set(submitted):raise ValueError('shadowed or disconnected leaf')
    return HorizonEnvelope(z,e,n,exact,ee,nn,len(submitted),dead,bounded)



def select_best_horizon_certificate(
    sequence: str, prefix: Path=ROOT, *, q: Fraction=Fraction(2),
    horizon: int=8, max_leaves: int=50000
) -> HorizonCertificate:
    """Choose the smallest *independently verified* upper bound among depths.

    Geometric restrictions can make a deeper block-product bound worse than
    its parent (block partition is not inherently nested). Do not assume
    monotonically improving partial refinements. Instead verify each depth
    and select its minimum exact-rational partition cap; every candidate
    is separately sound, so this choice preserves rigorous coverage.
    """
    if type(horizon) is not int or horizon<0:
        raise ValueError('horizon must be nonnegative')
    if not isinstance(q,Fraction) or q<1:
        raise ValueError('q must be an exact Fraction >= 1')
    best=None
    for d in range(min(horizon,len(sequence)-len(prefix))+1):
        candidate=build_horizon_certificate(sequence,prefix,horizon=d,
                                            max_leaves=max_leaves)
        score=verify_horizon_certificate(candidate,q).partition_upper
        if best is None or score<best[0]:best=(score,candidate)
    assert best is not None
    return best[1]
