"""Deterministic exact rational benchmark of conformational-switch robustness.

All model parameters are hypothetical rational weights, not measured energies.
"""
import json
from fractions import Fraction as F
from foldbounds.ensemble_crossover import probability,spectrum
from foldbounds.ensemble_crossover_3d import spectrum_3d
from foldbounds.ensemble_uncertainty import (
    certify_error_box,verify_error_box,certify_switch_band,verify_switch_band,
    critical_factor_ratio_for_majority,
)


def f(x):return str(x)
def bracket(pair):return {'lower':f(pair[0]),'upper':f(pair[1]),'mid_decimal':float(sum(pair)/2)}

def run(dim,seq,t):
    hist=spectrum(seq) if dim==2 else spectrum_3d(seq)
    q=F(3,2);blo=F(1,4);bhi=F(4)
    p_at_1=probability(hist,q,F(1),t)
    cases=[]
    for r in (F(1),F(11,10),F(6,5)):
        lo,hi=1/r,r
        cert=certify_error_box(seq,dim,t,q,F(2),blo,bhi,lo,hi)
        assert verify_error_box(cert)
        switch=certify_switch_band(seq,dim,t,q,blo,bhi,lo,hi)
        assert verify_switch_band(switch)
        cases.append({'per_conformation_factor_range':[f(lo),f(hi)],
            'relative_odds_distortion':f(hi/lo),
            'nominal_parameter_box':[f(cert.nominal_probability_lower),f(cert.nominal_probability_upper)],
            'adversarial_population_box':[f(cert.worst_probability_lower),f(cert.worst_probability_upper)],
            'crossing_guaranteed_complement_until':bracket(switch.guaranteed_complement_until),
            'nominal_crossover':bracket(switch.nominal_crossover),
            'crossing_guaranteed_event_from':bracket(switch.guaranteed_event_from),
            'undetermined_bend_width_decimal':float(
                  sum(switch.guaranteed_event_from)/2 - sum(switch.guaranteed_complement_until)/2),
            'certificate_positive_derivative_terms':cert.positive_derivative_terms})
    return {'dimension':dim,'sequence':seq,'event_contact_threshold':t,
            'states':sum(hist.values()),'histogram_bins':len(hist),
            'nominal_probability_at_q1_5_b1':f(p_at_1),
            'critical_odds_ratio_if_prob_above_half':
                f(critical_factor_ratio_for_majority(p_at_1)) if F(1,2)<p_at_1<F(1) else None,
            'cases':cases}


def main():
    result={'scientific_status':'Exact proofs conditional on enumerated toy Hamiltonian; NOT real-protein validated',
            'uncertainty_model':'Every state weight times any factor in [1/r,r]; common normalization cancels',
            'replications':[run(2,'HPHHPHHPHHHH',2),run(3,'HPHHPHHH',1)]}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()