"""Exact nonmonotone population certificate ablation, no random sampling."""
import json
from fractions import Fraction as F
from foldbounds.ensemble_crossover import robust_bend_sign_certificate,spectrum,probability
from foldbounds.ensemble_crossover_3d import spectrum_3d
from foldbounds.bernstein_population import (
    certify_bernstein_population,verify_bernstein_population,
    certified_bernstein_with_errors)

def one(dim,seq,t):
    hist=spectrum(seq) if dim==2 else spectrum_3d(seq)
    terms,failed,positive=robust_bend_sign_certificate(hist,t,F(3,2),F(1,4))
    # Concrete counterexample for global monotonicity: b=4 vs b=8, q=1.5.
    p4=probability(hist,F(3,2),F(4),t)
    p8=probability(hist,F(3,2),F(8),t)
    assert p4>p8
    result=dict(dimension=dim,sequence=seq,event_threshold=t,
                conformations=sum(hist.values()),
                sufficient_positive_derivative_certificate_available=not bool(failed),
                genuine_decrease_at_q_1_5={'b4':str(p4),'b8':str(p8),'difference':str(p4-p8)},
                interval_bounds=[])
    previous=None
    for depth in (0,1,2,3):
        r=certify_bernstein_population(seq,dim,t,F(3,2),F(2),F(1,4),F(8),depth)
        assert verify_bernstein_population(r)
        if previous is not None:
            assert r.lower_probability>=previous.lower_probability
            assert r.upper_probability<=previous.upper_probability
        previous=r
        elo,ehi=certified_bernstein_with_errors(r,F(5,6),F(6,5))
        result['interval_bounds'].append(dict(depth=depth,subboxes=r.subboxes,
                nominal_lower=str(r.lower_probability),nominal_upper=str(r.upper_probability),
                adversarial_lower=str(elo),adversarial_upper=str(ehi),
                robust_contact_event_minority=ehi<F(1,2)))
    return result

if __name__=='__main__':
    print(json.dumps({'model':'toy lattices; exact mass histogram, factor errors [5/6,6/5]',
        'cases':[one(2,'HHHPPHHH',1),one(3,'HHPPPPPH',1)]},indent=2,sort_keys=True))