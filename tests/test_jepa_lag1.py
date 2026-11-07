"""Independent covariance-oracle tests for exact stationary event lag one."""
from fractions import Fraction

import pytest

from foldbounds.jepa_flux_maxflow import (
    construct_reversible_kernel,propose_maxflow_certificate,
    stationary_event_lag1,stationary_event_lag1_reversible,
)
from foldbounds.jepa_reversible_kinetics import build_kernel,event_flux


@pytest.mark.parametrize("sequence",[
    "HHHHHH","HPHHPPHH","HHPPHHHH","HPHPHPHH",
])
def test_exact_lag_one_matches_direct_stationary_joint_covariance(sequence):
    kernel=build_kernel(sequence,Fraction(2))
    p,flux,_=event_flux(kernel,1)
    assert 0<p<1
    direct_joint=sum((
        kernel.weights[i]*probability
        for i,row in enumerate(kernel.transitions)
        if kernel.scores[i]>=1
        for j,probability in row
        if kernel.scores[j]>=1
    ),Fraction())/kernel.partition
    rho_from_joint=(direct_joint-p*p)/(p*(1-p))
    assert rho_from_joint==stationary_event_lag1(kernel,1)
    assert rho_from_joint==1-flux/(p*(1-p))


@pytest.mark.parametrize("sequence",[
    "HHHHHH","HPHHPPHH","HPHPHPHH",
])
def test_checked_mixture_decreases_event_lag_one_exactly(sequence):
    kernel=build_kernel(sequence,Fraction(2))
    certificate=propose_maxflow_certificate(kernel,1)
    baseline=stationary_event_lag1(kernel,1)
    mixed=construct_reversible_kernel(kernel,certificate,blend=Fraction(1,2))
    result=stationary_event_lag1_reversible(mixed,1)
    assert baseline is not None and result is not None
    assert result<=baseline
    # Direct independent covariance check on the mixed chain.
    p,_,_=event_flux(kernel,1)
    joint=sum((
        kernel.weights[i]*pr
        for i,row in enumerate(mixed.transitions)
        if kernel.scores[i]>=1
        for j,pr in row
        if kernel.scores[j]>=1
    ),Fraction())/kernel.partition
    assert result==(joint-p*p)/(p*(1-p))


@pytest.mark.parametrize("threshold",[-100,100])
def test_undefined_correlation_for_constant_event(threshold):
    kernel=build_kernel("HHHHHH",Fraction(2))
    assert stationary_event_lag1(kernel,threshold) is None
    mixture=construct_reversible_kernel(kernel,
        propose_maxflow_certificate(kernel,threshold))
    assert stationary_event_lag1_reversible(mixture,threshold) is None
