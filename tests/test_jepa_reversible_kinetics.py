"""Independent exact stationary-state, detailed-balance, and flux tests."""
from dataclasses import replace
from fractions import Fraction

import pytest

from foldbounds.jepa_reversible_kinetics import (
    all_states, build_kernel, class_flux_contributions, event_flux,
    move_classes, optimize_exact_flux, pivot, uniform_classes, verify_kernel,
)


@pytest.mark.parametrize("sequence", ["HHHHH", "HHHHHH", "HPHPHHH"])
@pytest.mark.parametrize("q", [Fraction(1), Fraction(2), Fraction(5, 3)])
def test_all_markov_rows_are_stochastic_and_exactly_reversible(sequence, q):
    kernel = build_kernel(sequence, q)
    verify_kernel(kernel)
    assert kernel.partition > 0
    assert len(kernel.states) > 1
    assert len(kernel.classes) == 2*(len(sequence)-2)
    assert sum(kernel.class_probs) == 1
    assert all(sum(p for j,p in row) == 1 for row in kernel.transitions)


@pytest.mark.parametrize("sequence", ["HHHHHH", "HPHHPPHH", "HPHPHPHH"])
@pytest.mark.parametrize("threshold", [1, 2])
def test_exact_flux_optimizer_cannot_decrease_uniform_flux(sequence,threshold):
    kernel = build_kernel(sequence,Fraction(3,2))
    result = optimize_exact_flux(kernel,threshold)
    selected = build_kernel(sequence,Fraction(3,2),
                            class_probs=result.optimized_weights)
    verify_kernel(selected)
    assert event_flux(selected,threshold)[1] == result.optimized_flux
    assert event_flux(kernel,threshold)[1] == result.baseline_flux
    assert result.optimized_flux >= result.baseline_flux
    assert all(w >= Fraction(1,12*(len(sequence)-2))
               for w in result.optimized_weights)
    assert sum(result.optimized_weights) == 1
    assert all(isinstance(v,Fraction) for v in result.class_flux)


def test_strictly_exact_stationary_probabilities_match_full_enumeration():
    from foldbounds.lattice import exact_spectrum
    sequence = "HHHHHHHH"
    q = Fraction(5,3)
    kernel = build_kernel(sequence,q)
    spectrum = exact_spectrum(sequence)
    partition = sum((count*q**score for score,count in spectrum.items()),Fraction())
    event_mass = sum((count*q**score for score,count in spectrum.items()
                      if score >= 1),Fraction())
    assert partition == kernel.partition
    assert event_mass/partition == event_flux(kernel,1)[0]


def test_reversed_crossing_flux_equals_forward_flux_at_equilibrium():
    kernel = build_kernel("HPHHPPHH",Fraction(2))
    z = kernel.partition
    positive = sum((
        kernel.weights[i]*probability
        for i,row in enumerate(kernel.transitions)
        if kernel.scores[i] >= 1
        for j,probability in row
        if kernel.scores[j] < 1
    ), Fraction())/z
    negative = sum((
        kernel.weights[i]*probability
        for i,row in enumerate(kernel.transitions)
        if kernel.scores[i] < 1
        for j,probability in row
        if kernel.scores[j] >= 1
    ), Fraction())/z
    assert positive == negative == event_flux(kernel,1)[1]


def test_bad_data_and_tampering_refused():
    kernel = build_kernel("HHHHHH",Fraction(2))
    with pytest.raises(ValueError):
        build_kernel("HHHHHH",Fraction(2),class_probs=(Fraction(1),))
    with pytest.raises(ValueError):
        build_kernel("HHHHHH",Fraction(2),
                     class_probs=(Fraction(0),)*len(kernel.classes))
    with pytest.raises(ValueError):
        build_kernel("HHHHHH",Fraction(2),
                     class_probs=(True,)*len(kernel.classes))
    with pytest.raises(ValueError):
        verify_kernel(replace(kernel,weights=(Fraction(999),)+kernel.weights[1:]))
    with pytest.raises(ValueError):
        verify_kernel(replace(kernel,transitions=((),)+kernel.transitions[1:]))


def test_finite_universe_limit_fails_closed_and_proposals_are_invertible():
    with pytest.raises(ValueError,match="max_states"):
        all_states("H"*11,max_states=10)
    path = all_states("HHHHHH")[0]
    for i in range(1,len(path)-1):
        forward = pivot(path,i,90)
        if forward is not None:
            assert pivot(forward,i,-90) == path
        half = pivot(path,i,180)
        if half is not None:
            assert pivot(half,i,180) == path


def test_zero_flux_is_honest_without_infinite_gain():
    kernel = build_kernel("PPPPPP",Fraction(2))
    result = optimize_exact_flux(kernel,1)
    assert result.baseline_flux == result.optimized_flux == 0
    assert result.relative_gain is None


def test_nonuniform_input_kernel_rechecks_class_flux_identity():
    n=8
    unif=uniform_classes(n)
    probs=[w/2 for w in unif]
    probs[0]+=Fraction(1,2)
    kernel=build_kernel("HHPPHHHH",Fraction(2),class_probs=tuple(probs))
    verify_kernel(kernel)
    selected=optimize_exact_flux(kernel,1)
    assert selected.baseline_flux >= 0
    assert selected.optimized_flux >= selected.baseline_flux
    assert len(class_flux_contributions(kernel,1))==len(move_classes(n))
