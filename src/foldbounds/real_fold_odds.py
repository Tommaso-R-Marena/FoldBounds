"""Exact robust fold-population constraints and transparent physical conversions.
No fitted predictor, force field, or estimated real-protein population.
"""
from fractions import Fraction
from decimal import Decimal
from math import log, sqrt

def _fraction(x):
    if not isinstance(x,(Fraction,int)) or isinstance(x,bool):
        raise TypeError('exact-rational input required')
    return Fraction(x)

def two_state_odds(p):
    p=_fraction(p)
    if not 0<p<1: raise ValueError('population must lie strictly within (0,1)')
    return p/(1-p)

def robust_population(p,faclo,fachi):
    """Sharp enclosure under independently adversarial statewise factors."""
    p=_fraction(p);l=_fraction(faclo);u=_fraction(fachi)
    if not 0<=p<=1 or l<=0 or u<l:raise ValueError('invalid population or factors')
    return (l*p/(l*p+u*(1-p)),u*p/(u*p+l*(1-p)))

def robust_interval(lo,hi,faclo,fachi):
    lo=_fraction(lo);hi=_fraction(hi)
    if not 0<=lo<=hi<=1:raise ValueError('not a probability enclosure')
    return robust_population(lo,faclo,fachi)[0],robust_population(hi,faclo,fachi)[1]

def minimum_cross_fold_odds_ratio(dominant_lower):
    """Minimal cross-mutant fold-A/fold-B odds ratio if dominance reverses.

    Valid with further unfolded/other states: competitor fractions are at most
    one minus dominant fraction. This constrains data, not sequence prediction.
    """
    p=_fraction(dominant_lower)
    if not Fraction(1,2)<p<1:raise ValueError('requires strict majority threshold')
    return (p/(1-p))**2

def min_gap_shift_kcal(dominant_lower,T_kelvin=293.15):
    """Necessary relative free-energy shift RT ln(odds ratio), kcal/mol.

    Assumes common T and equilibrium population fractions measured under
    comparable conditions. This is NOT a calculated binding/folding energy.
    """
    x=minimum_cross_fold_odds_ratio(dominant_lower)
    T=float(T_kelvin)
    if not 150<=T<=500:raise ValueError('temperature outside intended range')
    R_kcal=0.00198720425864083
    return R_kcal*T*log(float(x))

def max_admissible_independent_state_energy_error(nominal_p,target_p,T_kelvin=293.15):
    """Largest +/- delta kcal/mol per-state error permitting a robust >=target.

    Weight factors are exp(+/-delta/RT); hence ratio of adverse group errors
    exp(2delta/RT). Returns 0 if nominal <= target, not a calibrated error.
    """
    p=float(_fraction(nominal_p));target=float(_fraction(target_p));T=float(T_kelvin)
    if not 0<target<1 or not 0<p<1 or not 150<=T<=500:raise ValueError('invalid inputs')
    if p<=target:return 0.0
    relative_odds=(p/(1-p))/(target/(1-target))
    return 0.5*0.00198720425864083*T*log(relative_odds)

def report_known_threshold(T_kelvin=293.15)->dict:
    odds=minimum_cross_fold_odds_ratio(Fraction(9,10))
    p=Fraction(19,20);target=Fraction(9,10)
    error=max_admissible_independent_state_energy_error(p,target,T_kelvin)
    return {"min_cross_fold_odds_ratio_if_opposing_states_exceed_90pct":str(odds),
            "required_relative_fold_free_energy_shift_kcal_mol_at_T":min_gap_shift_kcal(Fraction(9,10),T_kelvin),
            "T_kelvin":T_kelvin,
            "nominal_95pct_to_certify_90pct_independent_per_state_energy_error_kcal_mol":error,
            "scope":"Conditional necessary lower bound from equilibrium odds. Not estimated GA95/GB95 energies."}