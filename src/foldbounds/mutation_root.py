"""Exact Descartes-isolated robust-majority q crossover (fixed bend weight).

Only enumerated 8- and 10-mer cubic HP compactness models are supported.
The certificate recomputes the complete rational event/complement polynomials.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from .mutation_switch import _mask, compact_spectrum
from .mutation_switch_10 import exact_spectrum10,_mask10


def _valid_hist(sequence:str):
    if len(sequence)==8:
        return compact_spectrum(_mask(sequence))
    if len(sequence)==10:
        _mask10(sequence)
        return exact_spectrum10(sequence)
    raise ValueError('only exact cubic 8-mer/10-mer models supported')


def _sign_changes(coeff):
    signs=[1 if x>0 else -1 for x in coeff if x]
    return sum(signs[i]!=signs[i-1] for i in range(1,len(signs)))


def _polynomial(hist,b:F,low:F,high:F):
    kdeg=max(k for (k,r,ev),v in hist)
    E=[F(0)]*(kdeg+1);C=[F(0)]*(kdeg+1)
    for (k,r,ev),count in hist:
        (E if ev else C)[k]+=count*b**r
    return tuple(low*E[k]-high*C[k] for k in range(kdeg+1))


def _eval(coeff,q):
    value=F(0)
    for c in reversed(coeff): value=value*q+c
    return value


@dataclass(frozen=True)
class CriticalContactRoot:
    sequence: str
    bend:F
    factor_lower:F
    factor_upper:F
    bracket_low:F
    bracket_high:F
    polynomial_coefficients:tuple[F,...]
    positive_root_sign_changes:int
    bisection_rounds:int


def certify_contact_root(sequence:str,bend:F,factor_lower:F,factor_upper:F,
                         q_low:F,q_high:F,bits:int=48)->CriticalContactRoot:
    if any(type(v) is not F for v in (bend,factor_lower,factor_upper,q_low,q_high)):
        raise ValueError('all inputs must be exact Fractions')
    if not 0<bend or not 0<factor_lower<=factor_upper or not 0<q_low<q_high:
        raise ValueError('strictly positive rational interval required')
    if type(bits) is not int or not 1<=bits<=60:
        raise ValueError('bisection steps must be integer 1..60')
    h=_valid_hist(sequence)
    D=_polynomial(h,bend,factor_lower,factor_upper)
    changes=_sign_changes(D)
    if changes!=1 or not _eval(D,q_low)<0<_eval(D,q_high):
        raise ValueError('unique positive robust-majority crossover not established')
    for _ in range(bits):
        mid=(q_low+q_high)/2
        if _eval(D,mid)<0:q_low=mid
        else:q_high=mid
    if not _eval(D,q_low)<0<_eval(D,q_high):
        raise AssertionError('rational sign isolator invariant violated')
    return CriticalContactRoot(sequence,bend,factor_lower,factor_upper,
                               q_low,q_high,D,changes,bits)


def verify_contact_root(cert:CriticalContactRoot)->bool:
    if not isinstance(cert,CriticalContactRoot):return False
    try:
        # A root bracket is checked directly rather than reconstructing an
        # assumed original 6..10 interval, permitting any valid isolator.
        if not (type(cert.bracket_low) is F and type(cert.bracket_high) is F and
                0<cert.bracket_low<cert.bracket_high and
                type(cert.bisection_rounds) is int and 1<=cert.bisection_rounds<=60):
            return False
        coeff=_polynomial(_valid_hist(cert.sequence),cert.bend,
                         cert.factor_lower,cert.factor_upper)
        return (cert.polynomial_coefficients==coeff and
                cert.positive_root_sign_changes==1==_sign_changes(coeff) and
                _eval(coeff,cert.bracket_low)<0<_eval(coeff,cert.bracket_high))
    except (TypeError,ValueError,ArithmeticError,ZeroDivisionError):return False