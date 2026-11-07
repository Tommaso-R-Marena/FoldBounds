"""JEPA exact-transition teacher and independently verified Boltzmann gate."""
from fractions import Fraction
import math

import pytest

from foldbounds.jepa_kinetic_teacher import (
    collect_pivot_examples, conditional_pivot_outcomes,
    conformation_features, fit_pivot_jepa, gate_model_proposal,
    pivot_class_features, rank_classes_from_model,
)
from foldbounds.jepa_reversible_kinetics import (
    build_kernel, class_flux_contributions, event_flux, optimize_exact_flux,
    verify_kernel,
)


def test_each_conditional_pivot_action_is_normalized():
    kernel=build_kernel("HPHHPPHH",Fraction(2))
    for state_index in range(0,len(kernel.states),37):
        for action in kernel.classes:
            distribution=conditional_pivot_outcomes(kernel,state_index,action)
            assert sum((p for j,p in distribution),Fraction())==1
            assert all(p>0 for _,p in distribution)
    assert len(conformation_features(kernel.sequence,kernel.q,1,kernel.states[0]))==10
    assert len(pivot_class_features(kernel.sequence,kernel.classes[0]))==9


def test_exact_expected_labels_not_empirical_or_forged():
    kernel=build_kernel("HHHHHH",Fraction(2))
    data=collect_pivot_examples(kernel,1,max_states=7,seed=0)
    assert len(data)==7*len(kernel.classes)
    for item in data:
        assert sum((p for _,p in item.successors),Fraction())==1
        assert Fraction()<=item.event_escape<=Fraction(1)
        assert all(math.isfinite(x) for x in item.state+item.action)
    assert data==collect_pivot_examples(kernel,1,max_states=7,seed=0)


def test_bad_model_vector_refused_before_policy_acceptance():
    kernel=build_kernel("HHHHHH",Fraction(2))
    with pytest.raises(ValueError):
        gate_model_proposal(kernel,1,(float("nan"),)*len(kernel.classes))
    with pytest.raises(ValueError):
        gate_model_proposal(kernel,1,(0.0,))
    with pytest.raises(ValueError):
        gate_model_proposal(kernel,1,(True,)*len(kernel.classes))


def test_weak_model_is_rejected_to_uniform_exactly():
    kernel=build_kernel("HPHPHPHH",Fraction(2))
    contributions=class_flux_contributions(kernel,1)
    worst=min(range(len(contributions)),key=lambda i:contributions[i])
    scores=tuple(1.0 if i==worst else 0.0 for i in range(len(contributions)))
    result=gate_model_proposal(kernel,1,scores)
    assert not result.accepted
    assert result.achieved_flux==result.baseline_flux
    assert result.raw_candidate_flux<result.baseline_flux
    verify_kernel(build_kernel(kernel.sequence,kernel.q,class_probs=result.actual_probs))


def test_oracle_policy_accepted_and_matches_exact_linear_optimum():
    kernel=build_kernel("HPHPHPHH",Fraction(2))
    oracle=optimize_exact_flux(kernel,1)
    scores=tuple(float(v) for v in oracle.class_flux)
    result=gate_model_proposal(kernel,1,scores)
    assert result.accepted
    assert result.achieved_flux==result.oracle_flux==oracle.optimized_flux


def test_jepa_learns_stochastic_successor_and_gate_stays_authoritative():
    pytest.importorskip("torch")
    train=build_kernel("HHHHHH",Fraction(2))
    data=collect_pivot_examples(train,1,max_states=15,seed=2)
    model,losses=fit_pivot_jepa(data,epochs=4,seed=2)
    assert len(losses)==4 and all(math.isfinite(value) for value in losses)
    test=build_kernel("HPHHPPHH",Fraction(2))
    predictions=rank_classes_from_model(model,test,1)
    assert len(predictions)==len(test.classes)
    assert all(math.isfinite(value) for value in predictions)
    gated=gate_model_proposal(test,1,predictions)
    assert gated.achieved_flux>=gated.baseline_flux
    assert gated.achieved_flux<=gated.oracle_flux


def test_small_exact_world_has_no_forged_kinetic_guarantees():
    kernel=build_kernel("PPPPPP",Fraction(2))
    result=gate_model_proposal(kernel,1,tuple(0.0 for _ in kernel.classes))
    assert result.baseline_flux==result.achieved_flux==result.oracle_flux==0
    assert event_flux(kernel,1)[0]==0
