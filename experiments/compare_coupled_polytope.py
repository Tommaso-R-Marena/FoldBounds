"""Matched checked probability-envelope ablation on small HP instances."""
from fractions import Fraction
from time import perf_counter
import json

from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead

CASES=(
    ('HPHHPHHPHHHH',Fraction(3,2),2),
    ('H'*8,Fraction(2),1),
    ('H'*12,Fraction(2),2),
    ('HPHHPHHP',Fraction(2),1),
)

def comparison():
    for seq,q,t in CASES:
        cert=build_frontier(seq,q,t,2)
        kwargs=dict(lookahead=0,use_matching=True,use_transfer=True,
                    use_monotone=True,use_turn_family=True)
        a=verify_lookahead(cert,**kwargs)
        b=verify_lookahead(cert,use_coupled_polytope=True,**kwargs)
        histogram=exact_spectrum(seq)
        e=sum((v*q**k for k,v in histogram.items() if k>=t),Fraction(0))
        c=sum((v*q**k for k,v in histogram.items() if k<t),Fraction(0))
        assert a.probability_lower<=e/(e+c)<=a.probability_upper
        assert b.probability_lower<=e/(e+c)<=b.probability_upper
        assert b.probability_lower>=a.probability_lower
        assert b.probability_upper<=a.probability_upper
        assert b.event_lower>=a.event_lower and b.event_upper<=a.event_upper
        assert b.complement_lower>=a.complement_lower and b.complement_upper<=a.complement_upper
        yield dict(sequence=seq,q=str(q),threshold=t,
                  baseline_lower=str(a.probability_lower),
                  baseline_upper=str(a.probability_upper),
                  coupled_lower=str(b.probability_lower),
                  coupled_upper=str(b.probability_upper),
                  exact_probability=str(e/(e+c)))

if __name__=='__main__':
    for r in comparison():print(json.dumps(r,sort_keys=True))