"""Exact symbolic ensemble phase certificates in a two-parameter 2D HP lattice model.

This is a toy polymer model, not a physical protein force field.
The mathematical certificate uses rational polynomial arithmetic and checks
conformation counts independently. No floating-point assumption is trusted.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from typing import NamedTuple

DIRECTIONS = ((1,0),(-1,0),(0,1),(0,-1))

class Geometry(NamedTuple):
    contacts: tuple[tuple[int,int], ...]
    bends: int

@lru_cache(maxsize=None)
def geometry(n: int) -> tuple[Geometry, ...]:
    if not 2 <= n <= 16:
        raise ValueError('n must be between 2 and 16; enumeration is exponential')
    result=[]
    def rec(path,last,bends):
        if len(path)==n:
            pairs=tuple((i,j) for i in range(n) for j in range(i+2,n)
                        if abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1)
            result.append(Geometry(pairs,bends)); return
        x,y=path[-1]
        for d in DIRECTIONS:
            nxt=(x+d[0],y+d[1])
            if nxt not in path:
                rec(path+(nxt,),d,bends+(d!=last))
    rec(((0,0),(1,0)),(1,0),0)
    return tuple(result)

def spectrum(seq: str) -> Counter:
    if type(seq) is not str or not seq or set(seq)-{'H','P'} or len(seq)<2:
        raise ValueError('sequence must be H/P string length >=2')
    h=Counter()
    for g in geometry(len(seq)):
        k=sum(seq[i]=='H' and seq[j]=='H' for i,j in g.contacts)
        h[(k,g.bends)]+=1
    return h

def add_poly(a,b):
    return tuple((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
                 for i in range(max(len(a),len(b))))

def mul_poly(a,b):
    z=[Fraction(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): z[i+j]+=x*y
    return tuple(z)

def deriv(a):return tuple((i*a[i] for i in range(1,len(a)))) or (Fraction(0),)

def eval_poly(p,b):
    if not isinstance(b,Fraction):raise TypeError('exact rational input required')
    acc=Fraction(0)
    for c in reversed(p):acc=acc*b+c
    return acc

def weighted_polys(h, q, threshold):
    if not isinstance(q,Fraction) or q<=0:raise ValueError('positive rational contact strength required')
    if type(threshold) is not int: raise ValueError('threshold must be integer')
    n=max(r for k,r in h)+1
    e=[Fraction(0)]*n;c=[Fraction(0)]*n
    for (k,r),count in h.items():
        (e if k>=threshold else c)[r]+=count*q**k
    return tuple(e),tuple(c)

def sign_changes(p):
    signs=[1 if x>0 else -1 for x in p if x]
    return sum(signs[i]!=signs[i-1] for i in range(1,len(signs)))

def unique_switch_bracket(h,q,t,b_low,b_high,bits=64):
    if not all(isinstance(v,Fraction) for v in (q,b_low,b_high)) or not 0<b_low<b_high:
        raise ValueError('positive exact fractions, increasing bend bracket')
    e,c=weighted_polys(h,q,t)
    diff=tuple(x-y for x,y in zip(e,c))
    if sign_changes(diff)!=1:
        raise ValueError('uniqueness NOT established by Descartes; cannot certify')
    vlow=eval_poly(diff,b_low);vhigh=eval_poly(diff,b_high)
    if not (vlow<0 and vhigh>0):raise ValueError('must straddle a negative-to-positive crossing')
    for _ in range(bits):
        mid=(b_low+b_high)/2
        if eval_poly(diff,mid)<0:b_low=mid
        else:b_high=mid
    return b_low,b_high

def monotone_in_bend_by_coefficients(h,q,t):
    """Sufficient exact certificate: each coefficient of E'C-EC' nonnegative.

    For any positive b, P'= (E'C-EC')/(E+C)^2. If all coefficients
    of the numerator are nonnegative, P is nondecreasing for b>0.
    """
    e,c=weighted_polys(h,q,t)
    left=mul_poly(deriv(e),c)
    right=mul_poly(e,deriv(c))
    delta=tuple(left[i]-right[i] for i in range(max(len(left),len(right))))
    return all(v>=0 for v in delta),delta

def probability(h,q,b,t):
    if not isinstance(b,Fraction) or b<=0:raise ValueError('positive exact rational bend weight required')
    e,c=weighted_polys(h,q,t)
    x,y=eval_poly(e,b),eval_poly(c,b)
    return x/(x+y)

if __name__=='__main__':
    h=spectrum('HPHHPHHPHHHH')
    for q in (Fraction(3,2),Fraction(2),Fraction(3)):
        e,c=weighted_polys(h,q,2)
        d=tuple(x-y for x,y in zip(e,c))
        ok,cc=monotone_in_bend_by_coefficients(h,q,2)
        print('q=',q,'states=',sum(h.values()),'signchanges=',sign_changes(d),'mono_bend_coeff=',ok,'neg_coeffs=',sum(v<0 for v in cc))
        print('sign at q,b=1/4,4:',eval_poly(d,Fraction(1,4)),eval_poly(d,Fraction(4)))
        if eval_poly(d,Fraction(1,4))<0<eval_poly(d,Fraction(4)):
            lo,hi=unique_switch_bracket(h,q,2,Fraction(1,4),Fraction(4),bits=32)
            print('switch bend in:',lo,hi,'floats=',float(lo),float(hi))

def robust_bend_sign_certificate(h,t,q_origin=Fraction(3,2),b_origin=Fraction(1,4)):
    from math import comb
    polynomial=Counter()
    ea=[(k,r,c) for (k,r),c in h.items() if k>=t]
    ca=[(k,r,c) for (k,r),c in h.items() if k<t]
    for ka,ra,na in ea:
        for kc,rc,nc in ca:
            if ra!=rc:
                polynomial[(ka+kc,ra+rc-1)] += (ra-rc)*na*nc
    shifted=Counter()
    for (powq,powb),coeff in polynomial.items():
        for iq in range(powq+1):
            x=coeff*comb(powq,iq)*q_origin**(powq-iq)
            for ib in range(powb+1):
                shifted[(iq,ib)] += x*comb(powb,ib)*b_origin**(powb-ib)
    negative=[(k,v) for k,v in shifted.items() if v<0]
    return len(shifted),negative[:8],sum(v>0 for v in shifted.values())
if __name__=='__main__':
  print('GLOBAL BEND DERIV:',robust_bend_sign_certificate(spectrum('HPHHPHHPHHHH'),2))

from dataclasses import dataclass

@dataclass(frozen=True)
class RobustPopulationEnvelope:
    """Proof-carrying box enclosure with *attained* probability endpoints.

    Sound only if all exact coefficients of the shifted E'C - EC' polynomial
    are nonnegative. This ensures bending monotonicity; contact monotonicity
    follows from the threshold event and the positive exponential family.
    """
    sequence: str
    threshold: int
    contact_lower: Fraction
    contact_upper: Fraction
    bend_lower: Fraction
    bend_upper: Fraction
    lower_probability: Fraction
    upper_probability: Fraction
    positive_derivative_terms: int
    zero_derivative_terms: int


def certify_rectangular_population(seq: str, threshold: int,
                                    contact_lower: Fraction, contact_upper: Fraction,
                                    bend_lower: Fraction, bend_upper: Fraction
                                    ) -> RobustPopulationEnvelope:
    if not all(isinstance(x,Fraction) for x in
               (contact_lower,contact_upper,bend_lower,bend_upper)):
        raise TypeError('all four parameter bounds must be exact Fractions')
    if not 0<contact_lower<=contact_upper or not 0<bend_lower<=bend_upper:
        raise ValueError('positive nonempty parameter box required')
    h=spectrum(seq)
    terms,negative,positive=robust_bend_sign_certificate(
        h,threshold,contact_lower,bend_lower)
    if negative:
        raise ValueError('nonnegative-coefficient bend certificate not established')
    low=probability(h,contact_lower,bend_lower,threshold)
    high=probability(h,contact_upper,bend_upper,threshold)
    return RobustPopulationEnvelope(seq,threshold,contact_lower,contact_upper,
                                    bend_lower,bend_upper,low,high,positive,terms-positive)


def verify_population_envelope(cert: RobustPopulationEnvelope) -> bool:
    if not isinstance(cert,RobustPopulationEnvelope):return False
    try:
        return certify_rectangular_population(
            cert.sequence,cert.threshold,cert.contact_lower,cert.contact_upper,
            cert.bend_lower,cert.bend_upper)==cert
    except (TypeError,ValueError):
        return False


def naive_decoupled_population_bounds(seq,threshold,qlo,qhi,blo,bhi):
    """Independent extrema of event and complement masses; typically loose."""
    h=spectrum(seq)
    el,cl=weighted_polys(h,qlo,threshold)
    eh,ch=weighted_polys(h,qhi,threshold)
    e_low,c_low=eval_poly(el,blo),eval_poly(cl,blo)
    e_high,c_high=eval_poly(eh,bhi),eval_poly(ch,bhi)
    return e_low/(e_low+c_high),e_high/(e_high+c_low)