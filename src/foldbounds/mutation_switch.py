"""Exactly certified 3-D HP lattice mutation-dependent compactness populations.

Research model (NOT all-atom or protein folding): fixed 8-bead cubic SAWs,
nonconsecutive geometric contact count >= 2 defines a compact-state event,
Boltzmann weights q**(HH contacts) * b**(chain bends). Sequence P/H mutates
without changing the geometric event; q,b vary over rational rectangles.

The untrusted proposer is checked by complete re-enumeration and independent
exact Bernstein polynomial arithmetic. Worst-case *state-specific* energy
reweighting is permitted without independence assumptions.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from typing import NamedTuple
from .ensemble_crossover_3d import geometry_3d
from .bernstein_population import _monomial_to_bernstein
from .ensemble_uncertainty import lower_probability, upper_probability

N = 8
COUNT = 13565


def _sequence(mask: int) -> str:
    if type(mask) is not int or not 0 <= mask < (1 << N):
        raise ValueError('sequence is an 8-bit H/P mask')
    return ''.join('H' if mask & (1 << i) else 'P' for i in range(N))


def _mask(sequence: str) -> int:
    if type(sequence) is not str or len(sequence) != N or set(sequence) - {'H', 'P'}:
        raise ValueError('required 8-character exact H/P sequence')
    return sum(1 << i for i, aa in enumerate(sequence) if aa == 'H')


@lru_cache(maxsize=1)
def geometric_fingerprints():
    """Exact lossless compression of the full 13,565-state geometric ensemble.

    A fingerprint encodes every nonconsecutive contact and the bend count.
    All states with one fingerprint have identical weights for *every* H/P
    sequence and every contact/bending parameter, so grouping commutes with
    mutation and Gibbs reweighting. Every full 3D state is still enumerated.
    """
    fingerprints=Counter()
    for pairs,bends in geometry_3d(N):
        edges=tuple(sorted(((1 << i) | (1 << j)) for i,j in pairs))
        fingerprints[(edges,bends,int(len(pairs)>=2))]+=1
    if sum(fingerprints.values())!=COUNT:
        raise AssertionError('geometry fingerprint compression is incomplete')
    return tuple(sorted(fingerprints.items()))


@lru_cache(maxsize=256)
def compact_spectrum(mask: int):
    """Count all 3D SAWs by HH energy, bends and geometric compactness.

    Fingerprint compression is *exact* for all mutation sequences, because
    each fingerprint retains every pair in the contact graph.
    """
    _sequence(mask)
    hist=Counter()
    for (edges,bends,is_compact),multiplicity in geometric_fingerprints():
        hh=sum((mask & edge)==edge for edge in edges)
        hist[(hh,bends,is_compact)]+=multiplicity
    if sum(hist.values()) != COUNT:
        raise AssertionError('anchored 3D geometry cardinality mismatch')
    return tuple(sorted(hist.items()))


def nominal_probability(mask: int, q: F, b: F) -> F:
    _sequence(mask)
    if type(q) is not F or type(b) is not F or q<=0 or b<=0:
        raise ValueError('positive rational model parameters required')
    hist=compact_spectrum(mask)
    a=sum((c*q**k*b**r for (k,r,compact),c in hist if compact),F(0))
    z=sum((c*q**k*b**r for (k,r,_),c in hist),F(0))
    return a/z


@dataclass(frozen=True)
class ExactPopulationBox:
    sequence: str
    q_lower: F
    q_upper: F
    bend_lower: F
    bend_upper: F
    event_lower: F
    event_upper: F
    microstates: int
    bernstein_coefficients: int
    refinement_depth: int


def _bernstein_one(hist, qa: F, qb: F, ba: F, bb: F):
    """Prove an exact convex-hull population enclosure for a closed box."""
    kd=max(key[0] for key,_ in hist)
    rd=max(key[1] for key,_ in hist)
    Q=[_monomial_to_bernstein(k,kd,qa,qb) for k in range(kd+1)]
    B=[_monomial_to_bernstein(r,rd,ba,bb) for r in range(rd+1)]
    E=[[F(0)]*(rd+1) for _ in range(kd+1)]
    C=[[F(0)]*(rd+1) for _ in range(kd+1)]
    for (k,r,ev),count in hist:
        dest=E if ev else C
        for i in range(kd+1):
            x=count*Q[k][i]
            if x:
                for j in range(rd+1):
                    dest[i][j]+=x*B[r][j]
    ratios=[]
    for i in range(kd+1):
        for j in range(rd+1):
            e,c=E[i][j],C[i][j]
            if e<0 or c<0 or e+c<=0:
                raise ArithmeticError('positive Bernstein weight cannot be proved')
            ratios.append(e/(e+c))
    return min(ratios), max(ratios), len(ratios)


def certify_population_box(sequence: str, q_lower: F, q_upper: F,
                           bend_lower: F, bend_upper: F,
                           depth: int=0) -> ExactPopulationBox:
    m=_mask(sequence)
    vals=(q_lower,q_upper,bend_lower,bend_upper)
    if not all(type(x) is F for x in vals) or not 0<q_lower<=q_upper or not 0<bend_lower<=bend_upper:
        raise ValueError('exact rational positive nonempty rectangle required')
    if type(depth) is not int or not 0<=depth<=3:
        raise ValueError('refinement depth must be integer 0..3')
    hist=compact_spectrum(m)
    boxes=[(q_lower,q_upper,bend_lower,bend_upper)]
    for _ in range(depth):
        nxt=[]
        for qa,qb,ba,bb in boxes:
            qm=(qa+qb)/2;bm=(ba+bb)/2
            nxt.extend(((qa,qm,ba,bm),(qa,qm,bm,bb),
                        (qm,qb,ba,bm),(qm,qb,bm,bb)))
        boxes=nxt
    bounds=[_bernstein_one(hist,*box) for box in boxes]
    return ExactPopulationBox(sequence,q_lower,q_upper,bend_lower,bend_upper,
                              min(v[0] for v in bounds),max(v[1] for v in bounds),
                              COUNT,bounds[0][2],depth)


def verify_population_box(cert: ExactPopulationBox) -> bool:
    if not isinstance(cert,ExactPopulationBox):return False
    try:
        return cert == certify_population_box(cert.sequence,cert.q_lower,
                    cert.q_upper,cert.bend_lower,cert.bend_upper,
                    cert.refinement_depth)
    except (TypeError,ValueError,ArithmeticError):return False


@dataclass(frozen=True)
class MutationSwitchCertificate:
    parent: ExactPopulationBox
    mutant: ExactPopulationBox
    mutated_residue: int   # one-indexed
    factor_lower: F
    factor_upper: F
    parent_robust_upper: F
    mutant_robust_lower: F
    orientation: str  # 'low-to-high' only (to make interpretation unambiguous)


def certify_mutation_switch(parent_sequence: str,mutant_sequence: str,
             q_lower:F,q_upper:F,bend_lower:F,bend_upper:F,
             factor_lower:F=F(5,6),factor_upper:F=F(6,5),
             depth:int=0)->MutationSwitchCertificate:
    p=_mask(parent_sequence);m=_mask(mutant_sequence)
    d=p^m
    if d==0 or d & (d-1):
        raise ValueError('mutation requires exactly one H/P substitution')
    if not all(type(x) is F for x in (factor_lower,factor_upper)) or not 0<factor_lower<=factor_upper:
        raise ValueError('invalid rational energy-error multipliers')
    before=certify_population_box(parent_sequence,q_lower,q_upper,bend_lower,bend_upper,depth)
    after=certify_population_box(mutant_sequence,q_lower,q_upper,bend_lower,bend_upper,depth)
    ub=upper_probability(before.event_upper,factor_lower,factor_upper)
    lb=lower_probability(after.event_lower,factor_lower,factor_upper)
    if not ub<F(1,2)<lb:
        raise ValueError('a robust minority-to-majority switch was not proved')
    return MutationSwitchCertificate(before,after,d.bit_length(),factor_lower,
                  factor_upper,ub,lb,'low-to-high')


def verify_mutation_switch(c:MutationSwitchCertificate)->bool:
    if not isinstance(c,MutationSwitchCertificate):return False
    try:
        p=c.parent;m=c.mutant
        if (p.q_lower,p.q_upper,p.bend_lower,p.bend_upper,p.refinement_depth)!=(m.q_lower,m.q_upper,m.bend_lower,m.bend_upper,m.refinement_depth):return False
        return c==certify_mutation_switch(p.sequence,m.sequence,p.q_lower,p.q_upper,
                        p.bend_lower,p.bend_upper,c.factor_lower,c.factor_upper,p.refinement_depth)
    except (TypeError,ValueError,ArithmeticError):return False


def enumerate_mutation_landscape(q:int=8,b:int=1,lo:F=F(5,6),hi:F=F(6,5)):
    """Exhaustively scan all 1,024 undirected one-residue sequence edges."""
    if type(q) is not int or q<1 or type(b) is not int or b<1:
        raise ValueError('positive integer q,b for exact sequence landscape scan')
    if type(lo) is not F or type(hi) is not F or not 0<lo<=hi:
        raise ValueError('positive exact error range required')
    p={i:nominal_probability(i,F(q),F(b)) for i in range(256)}
    minority={i:upper_probability(x,lo,hi)<F(1,2) for i,x in p.items()}
    majority={i:lower_probability(x,lo,hi)>F(1,2) for i,x in p.items()}
    flips=[]
    for a in range(256):
        for bit in range(N):
            z=a^(1<<bit)
            if z<=a:continue
            if minority[a] and majority[z]:flips.append((a,z,bit+1))
            if minority[z] and majority[a]:flips.append((z,a,bit+1))
    return tuple(sorted(flips))


@dataclass(frozen=True)
class UniformMutationLandscape:
    """Complete strict-majority classification of every one-site HP mutation.

    The ordered triples are (minority-mask, majority-mask, one-based site).
    Each of the 1,024 undirected edges is inspected exactly once.
    """
    q_lower: F
    q_upper: F
    bend_lower: F
    bend_upper: F
    factor_lower: F
    factor_upper: F
    robust_switch_edges: tuple[tuple[int,int,int],...]
    checked_edges: int
    microstates_per_sequence: int


def certify_uniform_mutation_landscape(q_lower:F,q_upper:F,
                             bend_lower:F,bend_upper:F,
                             factor_lower:F=F(5,6),factor_upper:F=F(6,5)):
    if not all(type(x) is F for x in (factor_lower,factor_upper)) or not 0<factor_lower<=factor_upper:
        raise ValueError('valid exact energy uncertainty factors required')
    prob=[]
    for m in range(1<<N):
        box=certify_population_box(_sequence(m),q_lower,q_upper,
                                   bend_lower,bend_upper)
        prob.append((lower_probability(box.event_lower,factor_lower,factor_upper),
                     upper_probability(box.event_upper,factor_lower,factor_upper)))
    flips=[];checked=0
    for original in range(1<<N):
        for site in range(N):
            mutant=original^(1<<site)
            if mutant<=original:continue
            checked+=1
            if prob[original][1]<F(1,2)<prob[mutant][0]:
                flips.append((original,mutant,site+1))
            elif prob[mutant][1]<F(1,2)<prob[original][0]:
                flips.append((mutant,original,site+1))
    if checked != N*(1<<(N-1)):
        raise AssertionError('full 8-residue sequence hypercube not scanned')
    return UniformMutationLandscape(q_lower,q_upper,bend_lower,bend_upper,
              factor_lower,factor_upper,tuple(sorted(flips)),checked,COUNT)


def verify_uniform_mutation_landscape(c:UniformMutationLandscape)->bool:
    if not isinstance(c,UniformMutationLandscape):return False
    try:
        return c==certify_uniform_mutation_landscape(c.q_lower,c.q_upper,
                  c.bend_lower,c.bend_upper,c.factor_lower,c.factor_upper)
    except (ValueError,TypeError,ArithmeticError):return False