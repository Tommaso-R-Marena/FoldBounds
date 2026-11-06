"""Length-10 independent 3D cubic-lattice mutational compactness validation.

This is a toy model: compact class = >=3 geometric nonconsecutive contacts.
The conformation set is enumerated exactly (308,981 anchored cubic SAWs).
No atomistic force field or biological phenotype is asserted.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from .ensemble_crossover_3d import geometry_3d
from .mutation_switch import _bernstein_one
from .ensemble_uncertainty import lower_probability, upper_probability

N=10
SAWS=308981
ORIGINAL='PHPHPHPHPH'
MUTANT='HHPHPHPHPH'


def _mask10(seq:str)->int:
    if type(seq) is not str or len(seq)!=N or set(seq)-{'H','P'}:
        raise ValueError('10-character H/P sequence required')
    return sum(1<<i for i,x in enumerate(seq) if x=='H')


@lru_cache(maxsize=1)
def geometry_fingerprints10():
    """Every entry encodes all contact pairs and bend count exactly."""
    freq=Counter()
    for pairs,bends in geometry_3d(N):
        code=tuple(sorted((1<<i)|(1<<j) for i,j in pairs))
        freq[(code,bends,int(len(pairs)>=3))]+=1
    if sum(freq.values())!=SAWS:
        raise AssertionError('missing anchored self-avoiding conformations')
    return tuple(sorted(freq.items()))


@lru_cache(maxsize=1024)
def exact_spectrum10(seq:str):
    bits=_mask10(seq)
    hist=Counter()
    for (edges,bends,compact),multiplicity in geometry_fingerprints10():
        hh=sum((bits & edge)==edge for edge in edges)
        hist[(hh,bends,compact)]+=multiplicity
    if sum(hist.values())!=SAWS:
        raise AssertionError('mass spectrum lost conformations')
    return tuple(sorted(hist.items()))


def population10(seq:str,q:F,b:F)->F:
    if type(q) is not F or type(b) is not F or q<=0 or b<=0:
        raise ValueError('positive rational q and b required')
    hist=exact_spectrum10(seq)
    event=sum((v*q**k*b**r for (k,r,c),v in hist if c),F())
    total=sum((v*q**k*b**r for (k,r,c),v in hist),F())
    return event/total


@dataclass(frozen=True)
class TenResidueMutationCertificate:
    original:str
    mutant:str
    q_lower:F
    q_upper:F
    bend_lower:F
    bend_upper:F
    correction_lower:F
    correction_upper:F
    original_nominal_lower:F
    original_nominal_upper:F
    mutant_nominal_lower:F
    mutant_nominal_upper:F
    original_robust_upper:F
    mutant_robust_lower:F
    conformations:int
    fingerprints:int


def certify_ten_mutation(q_lower:F=F(8),q_upper:F=F(12),
       bend_lower:F=F(3,4),bend_upper:F=F(5,4),
       factor_lower:F=F(5,6),factor_upper:F=F(6,5)):
    if not all(type(x) is F for x in (q_lower,q_upper,bend_lower,bend_upper,factor_lower,factor_upper)):
        raise ValueError('exact rational inputs required')
    if not 0<q_lower<=q_upper or not 0<bend_lower<=bend_upper or not 0<factor_lower<=factor_upper:
        raise ValueError('nonempty positive interval required')
    assert _mask10(ORIGINAL)^_mask10(MUTANT)==1
    a0,a1,ca=_bernstein_one(exact_spectrum10(ORIGINAL),q_lower,q_upper,bend_lower,bend_upper)
    b0,b1,cb=_bernstein_one(exact_spectrum10(MUTANT),q_lower,q_upper,bend_lower,bend_upper)
    after_min=lower_probability(b0,factor_lower,factor_upper)
    before_max=upper_probability(a1,factor_lower,factor_upper)
    if not before_max<F(1,2)<after_min:
        raise ValueError('specified parameter/error ranges do not certify a switch')
    return TenResidueMutationCertificate(ORIGINAL,MUTANT,q_lower,q_upper,
      bend_lower,bend_upper,factor_lower,factor_upper,a0,a1,b0,b1,
      before_max,after_min,SAWS,len(geometry_fingerprints10()))


def verify_ten_mutation(c:TenResidueMutationCertificate)->bool:
    if not isinstance(c,TenResidueMutationCertificate):return False
    try:
        return c==certify_ten_mutation(c.q_lower,c.q_upper,c.bend_lower,c.bend_upper,
                                       c.correction_lower,c.correction_upper)
    except (ValueError,TypeError,ArithmeticError):return False