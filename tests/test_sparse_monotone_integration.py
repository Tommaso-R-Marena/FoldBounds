"""Integration tests: checked global sparse enclosures, adaptive replay, input binding."""
from dataclasses import replace
from fractions import Fraction
import pytest
from foldbounds.lattice import ROOT, build_frontier, CertificateError
from foldbounds.lookahead import verify_lookahead
from foldbounds.adaptive import adaptive_refine,verify_refinement_run
from foldbounds.sparse_monotone import optimize_sparse_monotone,SparseMonotoneCertificate

SEQ='HPHHPHHPHHHH'
Q=Fraction(3,2)
THRESHOLD=2

def proposals(budget=25):
    return tuple(optimize_sparse_monotone(
        SEQ,ROOT,Q,THRESHOLD,leaf_budget=budget,
        max_core_length=9,objective=obj).certificate
        for obj in ('event','complement'))

def test_sparse_portfolio_strictly_improves_event_lower_without_widening_upper():
    cert=build_frontier(SEQ,Q,THRESHOLD,2)
    kw=dict(lookahead=0,use_matching=True,use_transfer=True,
            use_monotone=True,use_turn_family=True)
    old=verify_lookahead(cert,**kw)
    new=verify_lookahead(cert,**kw,sparse_certificates=proposals())
    assert new.event_lower>old.event_lower
    assert new.event_lower>=Fraction(98289,32)
    assert new.complement_lower>=old.complement_lower
    assert new.partition_upper==old.partition_upper
    assert new.event_upper==old.event_upper
    assert new.complement_upper==old.complement_upper
    assert new.probability_lower>=old.probability_lower
    assert new.probability_upper<=old.probability_upper

def test_exact_spectrum_contains_sparse_intervals():
    from foldbounds.lattice import exact_spectrum
    spectrum=exact_spectrum(SEQ)
    mass=sum((Fraction(count)*Q**score for score,count in spectrum.items()),Fraction(0))
    event=sum((Fraction(count)*Q**score for score,count in spectrum.items()
               if score>=THRESHOLD),Fraction(0))
    p=event/mass
    cert=build_frontier(SEQ,Q,THRESHOLD,2)
    b=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                       sparse_certificates=proposals(),use_coupled_polytope=True)
    assert b.event_lower<=event<=b.event_upper
    assert b.partition_lower<=mass<=b.partition_upper
    assert b.probability_lower<=p<=b.probability_upper

@pytest.mark.parametrize('invalid',[False,None,[],[1]])
def test_wrong_sparse_certificate_container_rejected(invalid):
    c=build_frontier(SEQ,Q,THRESHOLD,2)
    with pytest.raises(ValueError):
        verify_lookahead(c,lookahead=0,sparse_certificates=invalid)

def test_mismatched_model_parameters_rejected():
    cert=build_frontier(SEQ,Q,THRESHOLD,2)
    good=proposals(10)[0]
    for forged in (replace(good,q=Fraction(2)),
                   replace(good,event_at_least=0),
                   replace(good,sequence='H'*len(SEQ)),
                   replace(good,prefix=ROOT+((2,0),))):
        with pytest.raises(CertificateError):
            verify_lookahead(cert,lookahead=0,sparse_certificates=(forged,))

def test_independent_adaptive_replay_and_tamper_rejection():
    cert=build_frontier(SEQ,Q,THRESHOLD,2)
    checked=proposals(10)
    run=adaptive_refine(cert,steps=2,policy='mass',lookahead=0,
                        sparse_certificates=checked)
    assert run.bounds==verify_refinement_run(run)
    with pytest.raises(ValueError):
        verify_refinement_run(replace(run,sparse_certificates=()))
