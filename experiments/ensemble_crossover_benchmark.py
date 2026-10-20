"""Reproduce exactly certified HP lattice conformational-crossover results."""
from __future__ import annotations
from fractions import Fraction
from time import perf_counter
import json
from foldbounds.ensemble_crossover import (
    spectrum, probability, weighted_polys, sign_changes, eval_poly,
    robust_bend_sign_certificate, unique_switch_bracket,
    certify_rectangular_population, verify_population_envelope,
    naive_decoupled_population_bounds,
)

def main():
    seq='HPHHPHHPHHHH';threshold=2
    start=perf_counter()
    h=spectrum(seq)
    data={'model':'2D anchored HP lattice, HH contact reward q, bend propensity b',
          'sequence':seq,'event':'at least two nonconsecutive HH contacts',
          'states':sum(h.values()),'joint_histogram_bins':len(h),
          'crossover_brackets':{}}
    for q in (Fraction(3,2),Fraction(2),Fraction(3)):
        elo,ehi=unique_switch_bracket(h,q,threshold,Fraction(1,4),Fraction(4),bits=40)
        e,c=weighted_polys(h,q,threshold)
        f=tuple(x-y for x,y in zip(e,c))
        data['crossover_brackets'][str(q)]={'lower':str(elo),'upper':str(ehi),
             'approximately':float((elo+ehi)/2),'descartes_sign_changes':sign_changes(f),
             'verified_straddle':[eval_poly(f,elo)<0,eval_poly(f,ehi)>=0]}
    cert=certify_rectangular_population(seq,threshold,Fraction(3,2),Fraction(2),Fraction(1,4),Fraction(4))
    assert verify_population_envelope(cert)
    nl,nu=naive_decoupled_population_bounds(seq,threshold,Fraction(3,2),Fraction(2),Fraction(1,4),Fraction(4))
    data['certified_box']={'q':['3/2','2'],'b':['1/4','4'],
        'lower':str(cert.lower_probability),'upper':str(cert.upper_probability),
        'low_float':float(cert.lower_probability),'high_float':float(cert.upper_probability),
        'positive_derivative_coefficients':cert.positive_derivative_terms,
        'zero_derivative_coefficients':cert.zero_derivative_terms,
        'naive_lower':str(nl),'naive_upper':str(nu),
        'naive_width_float':float(nu-nl),
        'sharp_width_float':float(cert.upper_probability-cert.lower_probability)}
    high=certify_rectangular_population(seq,threshold,Fraction(3,2),Fraction(2),Fraction(1174,1000),Fraction(4))
    low=certify_rectangular_population(seq,threshold,Fraction(3,2),Fraction(2),Fraction(1,4),Fraction(647,1000))
    assert low.upper_probability < Fraction(1,2) < high.lower_probability
    data['guaranteed_majority']={'below_b':'647/1000','above_b':'1174/1000',
        'all_contact_q':'[3/2,2]','minority_upper':str(low.upper_probability),
        'majority_lower':str(high.lower_probability)}
    data['runtime_seconds']=round(perf_counter()-start,4)
    print(json.dumps(data,indent=2,sort_keys=True))

if __name__=='__main__':main()