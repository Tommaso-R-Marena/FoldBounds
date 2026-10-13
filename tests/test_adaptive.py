"""Adaptive search and proof-carrying zero-completion frontier terminals."""
from dataclasses import replace
from fractions import Fraction
import pytest

from foldbounds.lattice import (Certificate, CertificateError, Node,
                                build_frontier, exact_spectrum,
                                refine_frontier, verify_certificate, _neighbors)
from foldbounds.lookahead import verify_lookahead
from foldbounds.adaptive import adaptive_refine, verify_refinement_run


@pytest.mark.parametrize('policy', ['breadth','mass','probability'])
@pytest.mark.parametrize('sequence', ['HHHHHH', 'HHHHHHHH', 'HPHHPPHH'])
@pytest.mark.parametrize('threshold', [0,1,2])
def test_every_step_contains_exact_answer(policy,sequence,threshold):
    q=Fraction(3,2)
    spec=exact_spectrum(sequence)
    z=sum((count*q**k for k,count in spec.items()),Fraction(0))
    e=sum((count*q**k for k,count in spec.items() if k>=threshold),Fraction(0))
    initial=build_frontier(sequence,q,threshold,2)
    run=adaptive_refine(initial,steps=9,policy=policy,
                        lookahead=0,use_matching=True)
    assert 1<=len(run.history)<=10
    for step in run.history:
        assert step.partition_lower<=z<=step.partition_upper
        assert step.probability_lower<=e/z<=step.probability_upper
    for previous,next_step in zip(run.history,run.history[1:]):
        assert previous.partition_lower<=next_step.partition_lower
        assert next_step.partition_upper<=previous.partition_upper
        assert previous.probability_lower<=next_step.probability_lower
        assert next_step.probability_upper<=previous.probability_upper
    assert run.history[-1].probability_width == run.bounds.probability_upper-run.bounds.probability_lower


@pytest.mark.parametrize('policy',['breadth','mass','probability'])
def test_complete_search_converges_to_exact(policy):
    sequence='HHHHHH'
    q=Fraction(2)
    spec=exact_spectrum(sequence)
    z=sum((count*q**k for k,count in spec.items()),Fraction(0))
    event=sum((count*q**k for k,count in spec.items() if k>=1),Fraction(0))
    run=adaptive_refine(build_frontier(sequence,q,1,2),
                        steps=1000,policy=policy,use_matching=True)
    assert run.bounds.exact
    assert run.bounds.partition_lower==run.bounds.partition_upper==z
    assert run.bounds.probability_lower==run.bounds.probability_upper==event/z
    assert run.bounds.bounded_subtrees==0


def test_zero_completion_prefix_can_be_certified():
    # This nine-site spiral has no free next neighbor at its tip.
    dead=((0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2),(0,1),(1,1))
    assert not _neighbors(dead)
    sequence='H'*10
    cert=build_frontier(sequence,Fraction(2),1,9)
    assert Node(dead,'bound') in cert.frontier
    refined=refine_frontier(cert,dead)
    assert Node(dead,'dead') in refined.frontier
    old=verify_lookahead(cert,use_matching=True)
    new=verify_lookahead(refined,use_matching=True)
    baseline=verify_certificate(refined)
    assert old.partition_lower<=new.partition_lower<=new.partition_upper<=old.partition_upper
    exact=sum((count*Fraction(2)**k for k,count in exact_spectrum(sequence).items()),Fraction(0))
    assert new.partition_lower<=exact<=new.partition_upper
    assert baseline.partition_lower<=exact<=baseline.partition_upper


@pytest.mark.parametrize('kind',['dead'])
def test_fake_dead_node_is_rejected(kind):
    cert=build_frontier('HHHHHH',Fraction(2),1,2)
    tampered=Certificate(cert.sequence,cert.q,cert.event_at_least,
                         (Node(((0,0),(1,0)),kind),))
    with pytest.raises(CertificateError):
        verify_certificate(tampered)
    with pytest.raises(CertificateError):
        verify_lookahead(tampered,use_matching=True)


def test_determinism():
    cert=build_frontier('HHHHHHHH',Fraction(2),1,2)
    a=adaptive_refine(cert,steps=7,policy='probability')
    b=adaptive_refine(cert,steps=7,policy='probability')
    assert a==b


@pytest.mark.parametrize('steps', [-1,1.5,'7',True])
def test_invalid_step_budget(steps):
    cert=build_frontier('HHHH',Fraction(2),1,2)
    with pytest.raises(ValueError):
        adaptive_refine(cert,steps=steps)


def test_invalid_policy():
    cert=build_frontier('HHHH',Fraction(2),1,2)
    with pytest.raises(ValueError):
        adaptive_refine(cert,steps=2,policy='bad')


def test_adaptive_and_baseline_compatible():
    cert=build_frontier('HHPHHHHH',Fraction(2),1,2)
    run=adaptive_refine(cert,steps=16,policy='mass',lookahead=1)
    final=verify_lookahead(run.certificate,lookahead=1,use_matching=True)
    assert run.bounds.partition_lower>=final.partition_lower
    assert run.bounds.partition_upper<=final.partition_upper
    assert run.bounds.probability_lower>=final.probability_lower
    assert run.bounds.probability_upper<=final.probability_upper
    assert verify_refinement_run(run)==run.bounds



def test_replay_rejects_tampered_stored_bounds():
    cert=build_frontier('HHPHHHHH',Fraction(2),1,2)
    run=adaptive_refine(cert,steps=10,policy='probability')
    assert verify_refinement_run(run)==run.bounds
    forged=replace(run,bounds=replace(run.bounds,partition_upper=run.bounds.partition_upper+1))
    with pytest.raises(ValueError):
        verify_refinement_run(forged)


def test_initial_generator_covers_dead_ends_before_requested_depth():
    cert=build_frontier('H'*12,Fraction(2),1,11)
    assert any(n.kind=='dead' for n in cert.frontier)
    verify_certificate(cert)
    verify_lookahead(cert,use_matching=True)
