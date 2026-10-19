"""End-to-end model-relative exact-rational checks for optional coupling."""
from dataclasses import replace
from fractions import Fraction
import pytest

from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead
from foldbounds.adaptive import adaptive_refine, verify_refinement_run

@pytest.mark.parametrize('seq,q,threshold',[
    ('HHHHHHHH',Fraction(2),1),
    ('HHHHHHHH',Fraction(3,2),2),
    ('HPHHPHHP',Fraction(2),2),
    ('HPHHPHHPHHHH',Fraction(3,2),2),
    ('PPPPPPPP',Fraction(2),1),
])
def test_sharp_polytope_cannot_weaken_any_bound_and_contains_exact_oracle(seq,q,threshold):
    cert=build_frontier(seq,q,threshold,2)
    options=dict(lookahead=0,use_matching=True,use_transfer=True,
                 use_monotone=True,use_turn_family=True)
    old=verify_lookahead(cert,**options)
    new=verify_lookahead(cert,use_coupled_polytope=True,**options)
    assert new.event_lower>=old.event_lower
    assert new.event_upper<=old.event_upper
    assert new.complement_lower>=old.complement_lower
    assert new.complement_upper<=old.complement_upper
    assert new.partition_lower>=old.partition_lower
    assert new.partition_upper<=old.partition_upper
    assert new.probability_lower>=old.probability_lower
    assert new.probability_upper<=old.probability_upper
    spectrum=exact_spectrum(seq)
    e=sum((n*q**k for k,n in spectrum.items() if k>=threshold),Fraction(0))
    c=sum((n*q**k for k,n in spectrum.items() if k<threshold),Fraction(0))
    assert new.event_lower<=e<=new.event_upper
    assert new.complement_lower<=c<=new.complement_upper
    assert new.partition_lower<=e+c<=new.partition_upper
    assert new.probability_lower<=e/(e+c)<=new.probability_upper


def test_adaptive_polytope_is_replayed_from_mode_not_trusted_bounds():
    cert=build_frontier('HPHHPHHP',Fraction(3,2),1,2)
    run=adaptive_refine(cert,steps=3,policy='mass',use_matching=True,
                        use_transfer=True,use_coupled_polytope=True)
    assert run.use_coupled_polytope
    assert verify_refinement_run(run)==run.bounds
    assert run.bounds.probability_lower<=run.bounds.probability_upper


@pytest.mark.parametrize('bad',(None,0,1,'True',[],Fraction(1)))
def test_type_strictness(bad):
    cert=build_frontier('HHHH',Fraction(2),1,2)
    with pytest.raises(ValueError):
        verify_lookahead(cert,use_coupled_polytope=bad)