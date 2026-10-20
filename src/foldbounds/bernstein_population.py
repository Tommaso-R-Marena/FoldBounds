"""Exact bivariate Bernstein-ratio enclosures for finite HP ensembles.

Works for non-monotone population response surfaces; does not require a
positive-derivative sign witness. Every bound is exact rational arithmetic.
The finite histogram itself must be recomputed by a trusted finite model or
supplied via a separately certified mass-parameter polynomial.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from .ensemble_uncertainty import _valid_model,lower_probability,upper_probability


@lru_cache(maxsize=4096)
def _monomial_to_bernstein(exp:int,degree:int,l:F,h:F)->tuple[F,...]:
    """Exact coefficients of x**exp on [l,h], in Bernstein degree 'degree'."""
    if not (0<=exp<=degree and 0<l<=h):
        raise ValueError('positive interval and valid polynomial degree required')
    delta=h-l
    a=[F(comb(exp,k))*l**(exp-k)*delta**k for k in range(exp+1)]
    return tuple(sum((a[i]*F(comb(j,i),comb(degree,i))
                      for i in range(min(j,exp)+1)),F(0))
                 for j in range(degree+1))


def _bernstein_ratio_bounds(hist,t:int, qlo:F,qhi:F,blo:F,bhi:F):
    """Population must be inside ratios of positive coefficient pairs."""
    if not hist or any(type(v) is not int or v<=0 for v in hist.values()):
        raise ValueError('histogram coefficients must be positive exact counts')
    qdeg=max(k for k,r in hist)
    bdeg=max(r for k,r in hist)
    Qt={i:_monomial_to_bernstein(i,qdeg,qlo,qhi) for i in range(qdeg+1)}
    Bt={i:_monomial_to_bernstein(i,bdeg,blo,bhi) for i in range(bdeg+1)}
    event=[[F(0)]*(bdeg+1) for _ in range(qdeg+1)]
    comp=[[F(0)]*(bdeg+1) for _ in range(qdeg+1)]
    for (k,r),count in hist.items():
        target=event if k>=t else comp
        for i in range(qdeg+1):
            a=count*Qt[k][i]
            if not a:continue
            for j in range(bdeg+1):target[i][j]+=a*Bt[r][j]
    ratios=[]
    for er,cr in zip(event,comp):
        for e,c in zip(er,cr):
            if e<0 or c<0 or e+c<=0:
                raise AssertionError('positive Bernstein coefficient proof unavailable')
            ratios.append(e/(e+c))
    return min(ratios),max(ratios)


@dataclass(frozen=True)
class BernsteinPopulationCertificate:
    sequence:str
    dimension:int
    threshold:int
    q_lower:F
    q_upper:F
    bend_lower:F
    bend_upper:F
    depth:int
    microstates:int
    subboxes:int
    lower_probability:F
    upper_probability:F
    # Conservative enclosure valid across all q,b in the original rectangle.


def certify_bernstein_population(seq:str,dim:int,t:int,
                                  q_lower:F,q_upper:F,bend_lower:F,bend_upper:F,
                                  depth:int=2)->BernsteinPopulationCertificate:
    if any(type(x) is not F for x in (q_lower,q_upper,bend_lower,bend_upper)):
        raise TypeError('all parameter endpoints must be exact Fractions')
    if not 0<q_lower<=q_upper or not 0<bend_lower<=bend_upper:
        raise ValueError('invalid positive parameter rectangle')
    if type(depth) is not int or depth not in range(0,6):
        raise ValueError('depth must be an integer in 0..5 (exponential subdivision)')
    hist=_valid_model(seq,dim,t)
    boxes=[(q_lower,q_upper,bend_lower,bend_upper)]
    for _ in range(depth):
        nxt=[]
        for q0,q1,b0,b1 in boxes:
            qm=(q0+q1)/2;bm=(b0+b1)/2
            # Degenerate dimensions still produce duplicate but sound boxes.
            nxt.extend(((q0,qm,b0,bm),(q0,qm,bm,b1),
                        (qm,q1,b0,bm),(qm,q1,bm,b1)))
        boxes=nxt
    lowers=[];uppers=[]
    for q0,q1,b0,b1 in boxes:
        lo,hi=_bernstein_ratio_bounds(hist,t,q0,q1,b0,b1)
        lowers.append(lo);uppers.append(hi)
    return BernsteinPopulationCertificate(seq,dim,t,q_lower,q_upper,
        bend_lower,bend_upper,depth,sum(hist.values()),len(boxes),
        min(lowers),max(uppers))


def verify_bernstein_population(cert:BernsteinPopulationCertificate)->bool:
    if not isinstance(cert,BernsteinPopulationCertificate):return False
    try:
        return cert==certify_bernstein_population(cert.sequence,cert.dimension,
           cert.threshold,cert.q_lower,cert.q_upper,cert.bend_lower,
           cert.bend_upper,cert.depth)
    except (ValueError,TypeError,AssertionError,OverflowError,ZeroDivisionError):
        return False


def certified_bernstein_with_errors(cert:BernsteinPopulationCertificate,
                                    factor_lower:F,factor_upper:F):
    if not verify_bernstein_population(cert):
        raise ValueError('invalid Bernstein population certificate')
    return (lower_probability(cert.lower_probability,factor_lower,factor_upper),
            upper_probability(cert.upper_probability,factor_lower,factor_upper))