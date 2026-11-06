"""Exact per-pair Cartesian coordinate box geometry + all-network GNM audit.

A stricter source-only alternative to an isotropic distance-shift bound.
For nominal vector r (PDB in milliangstrom) and independent Cartesian
coordinate error ±eps per atom, every component of the difference vector
is independently in [r_k-2eps,r_k+2eps]. Hence EXACT extrema are:

  minimum_squared_distance = sum(max(abs(r_k)-2eps,0)**2)
  maximum_squared_distance = sum((abs(r_k)+2eps)**2)

These bounds are individually *attainable* for each atom pair, unlike
crude triangle bounds. An edge is certainly present only if max_sq<R²;
it is possibly present only if min_sq<R² (strict cutoff).

Different edges share coordinates, so independently choosing all their
statuses remains a graph OVERAPPROXIMATION. The complete finite-model GNM
rank scan uses small exact Woodbury inverse witnesses, no floats or labels.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .gnm_joint_rank import Edge,canonical_edges
from .gnm_reduced_woodbury import ReducedWorldCache


@dataclass(frozen=True)
class ExactBoxGraphs:
    n:int
    certainly_on:tuple[Edge,...]
    possibly_on:tuple[Edge,...]
    per_pair_optional:tuple[Edge,...]


@dataclass(frozen=True)
class ExactBoxRankAudit:
    n:int
    certain_nonlocal:int
    possible_nonlocal:int
    optional_contacts:int
    exact_graph_worlds:int
    globally_invariant:int
    opposite_sign_witnessed:int
    noninvariant_with_ties_only:int
    always_tied:int
    reduced_inverses_checked:int
    # Frozen signed witnesses for independently checkable counterexamples.
    # (i,j,positive_mask,negative_mask), with exact inverse checks upstream.
    flipping_witnesses:tuple[tuple[int,int,int,int],...]


def exact_coordinate_box_graphs(
    coordinates_milliangstrom:tuple[tuple[int,int,int],...],
    radius_angstrom:Fraction,epsilon_angstrom:Fraction,
)->ExactBoxGraphs:
    """Exact necessary/sufficient per-edge contact status for Cartesian boxes.

    Does NOT assert that every combination of optional edges is jointly
    geometrically realizable by a single coordinate perturbation.
    """
    coords=coordinates_milliangstrom
    if not isinstance(coords,tuple) or any(
        not isinstance(pt,tuple) or len(pt)!=3 or any(type(c) is not int for c in pt)
        for pt in coords
    ):
        raise ValueError('expected exact milliangstrom integer CA positions')
    n=len(coords)
    if type(radius_angstrom) is not Fraction or radius_angstrom<=0:
        raise ValueError('radius must be a strictly positive Fraction')
    if type(epsilon_angstrom) is not Fraction or epsilon_angstrom<0:
        raise ValueError('epsilon must be a nonnegative Fraction')
    # Validate n and mandatory chain via canonical_edges. All subsequent
    # extrema use only exact rational operations and integer source coords.
    backbone={(i,i+1) for i in range(n-1)}
    canonical_edges(n,backbone)
    radius_sq=(1000*radius_angstrom)**2
    two_eps=2000*epsilon_angstrom
    certain=set(backbone)
    possible=set(backbone)
    for i in range(n):
        for j in range(i+3,n):
            d=tuple(abs(coords[i][k]-coords[j][k]) for k in range(3))
            d_min_sq=sum(max(Fraction(0),x-two_eps)**2 for x in d)
            d_max_sq=sum((x+two_eps)**2 for x in d)
            if d_max_sq<radius_sq:
                certain.add((i,j))
            if d_min_sq<radius_sq:
                possible.add((i,j))
    if not certain<=possible:
        raise ArithmeticError('pairwise minimum/maximum geometry inconsistent')
    on=canonical_edges(n,certain)
    possibly=canonical_edges(n,possible)
    return ExactBoxGraphs(n,on,possibly,tuple(sorted(possible-certain)))


def audit_exact_box_gnm_ranks(graphs:ExactBoxGraphs,*,
                              max_worlds:int=8192)->ExactBoxRankAudit:
    """Complete exact ranking over every independently selected optional edge.

    Each graph uses ONE verified base full inverse plus a checked reduced
    integer (<=m)x(<=m) inverse; no stochastic sampling, random targets or
    floating-point acceptance. Fail closed on excessive work limits.
    """
    if not isinstance(graphs,ExactBoxGraphs):
        raise TypeError('expected exact pairwise geometry')
    n=graphs.n
    certain=canonical_edges(n,graphs.certainly_on)
    possible=canonical_edges(n,graphs.possibly_on)
    if not set(certain)<=set(possible):
        raise ValueError('malformed contact envelope')
    optional=tuple(sorted(set(possible)-set(certain)))
    if graphs.per_pair_optional!=optional:
        raise ValueError('stored optional contacts do not match graph envelope')
    if type(max_worlds) is not int or max_worlds<1 or (1<<len(optional))>max_worlds:
        raise ValueError('complete enumeration exceeds exact matrix budget')
    cache=ReducedWorldCache(n,certain,possible,max_optional=20)
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    pos=[-1]*len(pairs)
    neg=[-1]*len(pairs)
    zero=[False]*len(pairs)
    for mask in range(1<<len(optional)):
        diagonal=cache.get(mask).diagonal_numerators()
        for k,(i,j) in enumerate(pairs):
            delta=diagonal[i]-diagonal[j]
            if delta>0 and pos[k]<0:pos[k]=mask
            elif delta<0 and neg[k]<0:neg[k]=mask
            elif delta==0:zero[k]=True
    invariant=reverse=tie_only=always_tie=0
    witnesses=[]
    for (i,j),a,b,equal in zip(pairs,pos,neg,zero):
        if a>=0 and b>=0:
            reverse+=1
            witnesses.append((i,j,a,b))
        elif not equal and ((a>=0) ^ (b>=0)):
            invariant+=1
        elif a<0 and b<0:
            always_tie+=1
        else:
            tie_only+=1
    if invariant+reverse+tie_only+always_tie!=len(pairs):
        raise ArithmeticError('internal pair classification incomplete')
    return ExactBoxRankAudit(
        n,len(certain)-(n-1),len(possible)-(n-1),len(optional),
        1<<len(optional),invariant,reverse,tie_only,always_tie,
        cache.count,tuple(witnesses)
    )
