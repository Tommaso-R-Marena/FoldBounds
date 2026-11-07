"""Locally checkable state-dependent Hastings, exact receipts, frozen model."""
from dataclasses import replace
from fractions import Fraction
import math

import pytest

from foldbounds.jepa_local_hastings import (
    StabilityGuard,audit_finite_mh,exact_sample,mh_row,neural_policy,
    proposal_row,propose_edge_receipt,quantize_ranks,sample_pivot_chain,
    verify_edge_receipt,
)
from foldbounds.jepa_reversible_kinetics import (
    all_states,build_kernel,move_classes,uniform_classes,
)


def uniform_policy(path):
    return uniform_classes(len(path))


def structural_policy(path):
    """Positive, deterministic, state-dependent rational heuristic."""
    classes=move_classes(len(path))
    weight=[(3+abs(path[k][0])+int(family=="half")*len(path))
            for k,family in classes]
    z=sum(weight)
    return tuple(Fraction(x,z) for x in weight)


@pytest.mark.parametrize("sequence",["HHHHHH","HPHHPPHH"])
def test_local_mh_agrees_exactly_with_classical_metropolis_when_uniform(sequence):
    kernel=build_kernel(sequence,Fraction(2))
    for i,path in enumerate(kernel.states):
        actual=mh_row(sequence,Fraction(2),path,uniform_policy)
        expected={kernel.states[j]:p for j,p in kernel.transitions[i]}
        assert {x:p for x,p in actual.items() if p}==expected
    audit=audit_finite_mh(sequence,Fraction(2),uniform_policy)
    assert audit.states==len(kernel.states)
    assert audit.partition==kernel.partition


@pytest.mark.parametrize("sequence",["HHHHHH","HHPHHHHH","HPHHPPHH"])
def test_state_dependent_nonuniform_heuristic_preserves_boltzmann(sequence):
    audit=audit_finite_mh(sequence,Fraction(5,3),structural_policy)
    assert audit.states==len(all_states(sequence))
    assert audit.checked_edges>=audit.states
    assert audit.event_mass>=0


def test_edge_receipt_is_independently_checked_without_enumeration():
    sequence="HPHHPPHH"
    q=Fraction(2)
    source=all_states(sequence)[0]
    candidates=[y for y in proposal_row(source,structural_policy) if y!=source]
    assert candidates
    chosen=candidates[0]
    receipt=propose_edge_receipt(sequence,q,source,chosen,structural_policy)
    assert verify_edge_receipt(sequence,q,receipt,structural_policy)>0
    with pytest.raises(ValueError):
        verify_edge_receipt(sequence,q,replace(
            receipt,transition_forward=receipt.transition_forward+1,
        ),structural_policy)
    with pytest.raises(ValueError):
        propose_edge_receipt(sequence,q,source,source,structural_policy)


def test_untrusted_scores_are_deterministically_quantized_to_positive_rationals():
    scores=(0.25,0.25,-10.0,8.0)
    weights=quantize_ranks(scores)
    assert weights==quantize_ranks(scores)
    assert sum(weights)==1
    assert all(isinstance(w,Fraction) and w>0 for w in weights)
    assert weights[3]>weights[0]>weights[1]>weights[2]
    with pytest.raises(ValueError):
        quantize_ranks((float("nan"),))
    with pytest.raises(ValueError):
        quantize_ranks((True,))
    with pytest.raises(ValueError):
        quantize_ranks(())


def test_stability_guard_rejects_mutating_policy():
    original=uniform_policy
    changed=[False]
    def changing(path):
        if not changed[0]:
            return original(path)
        return structural_policy(path)
    guard=StabilityGuard(changing)
    source=all_states("HHHHHH")[0]
    guard(source)
    changed[0]=True
    with pytest.raises(ValueError,match="changed"):
        guard(source)


def test_invalid_policy_and_stochasticity_fail_closed():
    path=all_states("HHHHHH")[0]
    with pytest.raises(ValueError):
        proposal_row(path,lambda x:(Fraction(1),))
    with pytest.raises(ValueError):
        proposal_row(path,lambda x:(Fraction(0),)*len(move_classes(len(x))))
    with pytest.raises(ValueError):
        proposal_row(path,lambda x:(True,)*len(move_classes(len(x))))


def test_seeded_exact_sampler_and_longer_sequence_without_enumeration():
    seq="HHPPHHPP"*4
    policy=StabilityGuard(structural_policy)
    first=sample_pivot_chain(seq,Fraction(2),policy,steps=4,seed=8)
    second=sample_pivot_chain(seq,Fraction(2),policy,steps=4,seed=8)
    assert first==second
    assert len(first)==5
    assert all(len(path)==len(seq) and len(set(path))==len(path) for path in first)


def test_frozen_jepa_actor_can_propose_exact_rational_actions():
    pytest.importorskip("torch")
    from foldbounds.jepa_kinetic_teacher import collect_pivot_examples,fit_pivot_jepa
    training=build_kernel("HHHHHH",Fraction(2))
    examples=collect_pivot_examples(training,1,max_states=10,seed=4)
    model,loss=fit_pivot_jepa(examples,epochs=3,seed=4)
    assert all(math.isfinite(x) for x in loss)
    policy=StabilityGuard(neural_policy(model,"HPHHPPHH",Fraction(2),1))
    # Exact finite oracle audits the entire learned state-dependent kernel.
    report=audit_finite_mh("HPHHPPHH",Fraction(2),policy,max_states=1000)
    assert report.states==543
    sample=sample_pivot_chain("HPHHPPHH",Fraction(2),policy,steps=3,seed=7)
    assert len(sample)==4
