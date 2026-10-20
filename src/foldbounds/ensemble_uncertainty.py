"""Certified sensitivity of finite HP conformational populations to energy errors.

All calculations use exact rational weight factors and an independently
recomputed finite lattice histogram. No assumptions about a physical force
field, fitted parameters, or calibrated error bars are made.

A model weight is w(x)=q**contacts(x)*b**bends(x). An unknown statewise
multiplicative correction f(x) lives in [factor_low,factor_high]. A common
factor cancels; worst-case fold-event odds multiply by factor_high/factor_low.
The relative-weight correction can be correlated arbitrarily across states:
this worst-case interval is valid even without independence assumptions.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from typing import Literal
from .ensemble_crossover import (
    spectrum, probability, robust_bend_sign_certificate, weighted_polys,
    eval_poly, sign_changes,
)
from .ensemble_crossover_3d import spectrum_3d

Dimension = Literal[2,3]


def _strict_fraction(x, name):
    if type(x) is not F:
        raise TypeError(f'{name} must be exact Fraction')
    return x


def _valid_model(seq, dim, threshold):
    if type(dim) is not int or dim not in (2,3):
        raise ValueError('dimension must be 2 or 3')
    if type(threshold) is not int:
        raise ValueError('threshold must be integer')
    # Individual full-state traversal is exponential, intentionally guarded.
    if dim == 2:
        if type(seq) is not str or not 2 <= len(seq) <= 14:
            raise ValueError('2D benchmark length must be 2..14')
        return spectrum(seq)
    if type(seq) is not str or not 2 <= len(seq) <= 9:
        raise ValueError('3D benchmark length must be 2..9')
    return spectrum_3d(seq)


def lower_probability(p: F, factor_low: F, factor_high: F) -> F:
    """Sharp worst-case event probability if relative weights are uncertain."""
    for name,v in (('p',p),('factor_low',factor_low),('factor_high',factor_high)):
        _strict_fraction(v,name)
    if not F(0)<=p<=F(1) or not F(0)<factor_low<=factor_high:
        raise ValueError('invalid population or positive weight correction interval')
    return factor_low*p/(factor_low*p + factor_high*(1-p))


def upper_probability(p: F, factor_low: F, factor_high: F) -> F:
    for name,v in (('p',p),('factor_low',factor_low),('factor_high',factor_high)):
        _strict_fraction(v,name)
    if not F(0)<=p<=F(1) or not F(0)<factor_low<=factor_high:
        raise ValueError('invalid population or positive weight correction interval')
    return factor_high*p/(factor_high*p + factor_low*(1-p))


@dataclass(frozen=True)
class CertifiedErrorBox:
    sequence: str
    dimension: int
    threshold: int
    q_lower: F
    q_upper: F
    bend_lower: F
    bend_upper: F
    factor_lower: F
    factor_upper: F
    microstates: int
    derivative_terms: int
    positive_derivative_terms: int
    nominal_probability_lower: F
    nominal_probability_upper: F
    worst_probability_lower: F
    worst_probability_upper: F
    status: str


def certify_error_box(seq: str, dim: int, threshold: int,
                      q_lower: F,q_upper: F,bend_lower: F,bend_upper: F,
                      factor_lower: F,factor_upper: F) -> CertifiedErrorBox:
    """Provably sharp envelope over *all* allowed statewise corrections.

    If the derivative numerator has no negative coefficients after shifting
    to the box's lower q and bend corner, P is monotone in bending throughout
    the box. Event threshold ensures P is also monotone in contact strength.
    Thus extrema occur at corners. For each corner, choosing all event-state
    factors at one endpoint, all complement factors at the opposite endpoint
    attains the worst-case bound. Fully adversarial corrections are permitted.
    """
    vals=(('q_lower',q_lower),('q_upper',q_upper),
          ('bend_lower',bend_lower),('bend_upper',bend_upper),
          ('factor_lower',factor_lower),('factor_upper',factor_upper))
    for name,v in vals: _strict_fraction(v,name)
    if not F(0)<q_lower<=q_upper or not F(0)<bend_lower<=bend_upper:
        raise ValueError('nonempty positive q/bend box required')
    if not F(0)<factor_lower<=factor_upper:
        raise ValueError('nonempty positive multiplicative error bounds required')
    hist=_valid_model(seq,dim,threshold)
    terms,negative,positive=robust_bend_sign_certificate(
        hist,threshold,q_lower,bend_lower)
    if negative:
        raise ValueError('bend-monotonicity coefficient proof not established')
    p_lo=probability(hist,q_lower,bend_lower,threshold)
    p_hi=probability(hist,q_upper,bend_upper,threshold)
    assert p_lo<=p_hi
    worst_lo=lower_probability(p_lo,factor_lower,factor_upper)
    worst_hi=upper_probability(p_hi,factor_lower,factor_upper)
    status=('robust-majority-event' if worst_lo>F(1,2) else
            'robust-majority-complement' if worst_hi<F(1,2) else
            'undetermined')
    return CertifiedErrorBox(seq,dim,threshold,q_lower,q_upper,
        bend_lower,bend_upper,factor_lower,factor_upper,sum(hist.values()),
        terms,positive,p_lo,p_hi,worst_lo,worst_hi,status)


def verify_error_box(obj: CertifiedErrorBox) -> bool:
    if not isinstance(obj,CertifiedErrorBox):return False
    try:
        expected=certify_error_box(obj.sequence,obj.dimension,obj.threshold,
            obj.q_lower,obj.q_upper,obj.bend_lower,obj.bend_upper,
            obj.factor_lower,obj.factor_upper)
        return expected==obj
    except (TypeError,ValueError,AssertionError,ZeroDivisionError,OverflowError):
        return False


def _polynomial_difference(h: Counter, q: F, threshold: int,
                           event_factor: F, complement_factor: F):
    e,c=weighted_polys(h,q,threshold)
    return tuple(event_factor*(e[i] if i<len(e) else 0) -
                 complement_factor*(c[i] if i<len(c) else 0)
                 for i in range(max(len(e),len(c))))


def _unique_positive_root_bracket(poly, lo: F, hi: F, steps: int=65):
    if type(steps) is not int or not 1<=steps<=256:
        raise ValueError('bounded positive number of bisection steps required')
    if sign_changes(poly)!=1:
        raise ValueError('unique positive root is not certified by Descartes')
    if not eval_poly(poly,lo)<0<eval_poly(poly,hi):
        raise ValueError('positive root not strictly bracketed')
    a,b=lo,hi
    for _ in range(steps):
        m=(a+b)/2
        v=eval_poly(poly,m)
        if v<0:a=m
        else:b=m
    return a,b


@dataclass(frozen=True)
class RobustSwitchBand:
    sequence: str
    dimension: int
    threshold: int
    contact_strength: F
    factor_lower: F
    factor_upper: F
    bend_domain_lower: F
    bend_domain_upper: F
    guaranteed_complement_until: tuple[F,F]
    guaranteed_event_from: tuple[F,F]
    nominal_crossover: tuple[F,F]
    bits: int
    # All brackets are exact rational pairs, each containing a unique root.


def certify_switch_band(seq: str, dim: int, threshold: int, q: F,
                        bend_lower: F, bend_upper: F,
                        factor_lower: F, factor_upper: F,
                        bits: int=64) -> RobustSwitchBand:
    """Bracket both ends of adversarial-error-induced fold-switch ambiguity.

    Below the lower root, even maximizing event weights cannot give
    a majority. Above the upper root, even minimizing event weights cannot
    prevent a majority. Inside band both are possible at some allowed factors.
    Contact/bend model and correction factors are *hypotheses*, not measured.
    """
    for name,v in (('q',q),('bend_lower',bend_lower),('bend_upper',bend_upper),
                   ('factor_lower',factor_lower),('factor_upper',factor_upper)):
        _strict_fraction(v,name)
    if not (q>0 and 0<bend_lower<bend_upper and 0<factor_lower<=factor_upper):
        raise ValueError('positive ordered rational inputs required')
    hist=_valid_model(seq,dim,threshold)
    terms,negative,positive=robust_bend_sign_certificate(
        hist,threshold,q,bend_lower)
    if negative:
        raise ValueError('cannot certify monotonicity throughout bend range')
    # All errors in [l,u]. To *rule out* event majority require uE<lC;
    # to *guarantee* event majority require lE>uC.
    minority_poly=_polynomial_difference(hist,q,threshold,
                                          factor_upper,factor_lower)
    majority_poly=_polynomial_difference(hist,q,threshold,
                                          factor_lower,factor_upper)
    nominal_poly=_polynomial_difference(hist,q,threshold,F(1),F(1))
    mid_root=_unique_positive_root_bracket(nominal_poly,bend_lower,bend_upper,bits)
    factor_low_equal = (factor_lower == factor_upper)
    if factor_low_equal:
        # Identical correction for all conformations cancels exactly.
        lo_root=hi_root=mid_root
    else:
        lo_root=_unique_positive_root_bracket(minority_poly,bend_lower,bend_upper,bits)
        hi_root=_unique_positive_root_bracket(majority_poly,bend_lower,bend_upper,bits)
    if not (factor_low_equal or lo_root[1]<=mid_root[0]<=mid_root[1]<=hi_root[0]):
        raise AssertionError('uncertainty crossover roots out of order')
    return RobustSwitchBand(seq,dim,threshold,q,factor_lower,factor_upper,
        bend_lower,bend_upper,lo_root,hi_root,mid_root,bits)


def verify_switch_band(obj: RobustSwitchBand) -> bool:
    if not isinstance(obj,RobustSwitchBand):return False
    try:
        return obj==certify_switch_band(obj.sequence,obj.dimension,obj.threshold,
            obj.contact_strength,obj.bend_domain_lower,obj.bend_domain_upper,
            obj.factor_lower,obj.factor_upper,bits=obj.bits)
    except (ValueError,TypeError,AssertionError,ZeroDivisionError,OverflowError):
        return False


def critical_factor_ratio_for_majority(p:F) -> F:
    """Sharp maximum *relative event-to-complement* odds distortion.

    For p>1/2, a guarantee of event majority survives only if
    (factor_high/factor_low) < p/(1-p). Equality permits a tie.
    """
    _strict_fraction(p,'p')
    if not F(1,2)<p<F(1):raise ValueError('requires strict uncertain majority p')
    return p/(1-p)


@dataclass(frozen=True)
class CertifiedMassIntervals:
    """Conditional input: masses must come from a separately trusted checker."""
    event_lower: F
    event_upper: F
    complement_lower: F
    complement_upper: F


@dataclass(frozen=True)
class ReweightedPopulationInterval:
    """Guaranteed enclosure assuming input mass bounds and factor bounds hold."""
    population_lower: F
    population_upper: F


def reweight_certified_mass_intervals(
        masses: CertifiedMassIntervals, factor_lower: F, factor_upper: F
        ) -> ReweightedPopulationInterval:
    """Propagate a separately *certified* E/C enclosure into energy uncertainty.

    This layer never certifies the upstream masses itself. It works even when
    complete microstate enumeration is infeasible. Degenerate denominator
    situations return conservative 0 or 1 instead of division by zero.
    """
    if not isinstance(masses,CertifiedMassIntervals):
        raise TypeError('separately certified event/complement mass bounds required')
    for name,v in vars(masses).items():_strict_fraction(v,name)
    _strict_fraction(factor_lower,'factor_lower')
    _strict_fraction(factor_upper,'factor_upper')
    le,ue,lc,uc=(masses.event_lower,masses.event_upper,
                 masses.complement_lower,masses.complement_upper)
    if not (0<=le<=ue and 0<=lc<=uc and 0<factor_lower<=factor_upper):
        raise ValueError('invalid positive mass or correction factor intervals')
    if ue+uc<=0:
        raise ValueError('upper mass bound proves zero partition function')
    denom_lo=factor_lower*le+factor_upper*uc
    denom_hi=factor_upper*ue+factor_lower*lc
    p_lo=factor_lower*le/denom_lo if denom_lo>0 else F(0)
    p_hi=factor_upper*ue/denom_hi if denom_hi>0 else F(1)
    return ReweightedPopulationInterval(p_lo,p_hi)