"""Exact primal/dual cut certificates and reversible witness regression."""
from dataclasses import replace
from fractions import Fraction

import pytest

from foldbounds.jepa_reversible_kinetics import build_kernel,event_flux
from foldbounds.jepa_flux_maxflow import (
    propose_maxflow_certificate,verify_maxflow_certificate,
    construct_reversible_kernel,verify_reversible_kernel,event_flux_reversible,
)


@pytest.mark.parametrize("sequence",[
    "HHHHHH","HPHHPPHH","HHPPHHHH","HPHPHPHH",
])
@pytest.mark.parametrize("q",[Fraction(1),Fraction(2),Fraction(5,3)])
def test_exact_maxflow_cuts_are_dual_optimal_and_preserve_equilibrium(sequence,q):
    kernel=build_kernel(sequence,q)
    cert=propose_maxflow_certificate(kernel,1)
    verify_maxflow_certificate(kernel,cert)
    mixed=construct_reversible_kernel(kernel,cert,blend=Fraction(1,2))
    verify_reversible_kernel(mixed)
    baseline=event_flux(kernel,1)[1]
    optimum=cert.total_mass/kernel.partition
    achieved=event_flux_reversible(mixed,1)
    assert optimum>=baseline
    assert achieved==(optimum+baseline)/2
    assert achieved>=baseline
    assert cert.total_mass==cert.cut_mass


@pytest.mark.parametrize("threshold",[0,1,2,100])
def test_event_or_complement_empty_still_certifies_zero_flux(threshold):
    kernel=build_kernel("PPPPPP",Fraction(2))
    cert=propose_maxflow_certificate(kernel,threshold)
    verify_maxflow_certificate(kernel,cert)
    mixed=construct_reversible_kernel(kernel,cert)
    assert event_flux_reversible(mixed,threshold)==0
    assert cert.total_mass==0


def test_malicious_flow_and_cut_certificates_fail_closed():
    kernel=build_kernel("HPHHPPHH",Fraction(2))
    cert=propose_maxflow_certificate(kernel,1)
    assert cert.flows
    i,j,m=cert.flows[0]
    with pytest.raises(ValueError):
        verify_maxflow_certificate(kernel,replace(
            cert,flows=((i,j,m+1),)+cert.flows[1:],
        ))
    with pytest.raises(ValueError):
        verify_maxflow_certificate(kernel,replace(
            cert,cut_mass=cert.cut_mass+1,
        ))
    with pytest.raises(ValueError):
        verify_maxflow_certificate(kernel,replace(
            cert,source_side_states=(),
        ))
    with pytest.raises(ValueError):
        verify_maxflow_certificate(kernel,replace(
            cert,flows=((j,i,m),)+cert.flows[1:],
        ))
    with pytest.raises(ValueError):
        verify_maxflow_certificate(kernel,replace(
            cert,flows=((i,j,m),)+cert.flows,
        ))


def test_falsified_reversible_witness_and_invalid_mixing_rejected():
    kernel=build_kernel("HPHHPPHH",Fraction(2))
    cert=propose_maxflow_certificate(kernel,1)
    mixed=construct_reversible_kernel(kernel,cert)
    with pytest.raises(ValueError):
        verify_reversible_kernel(replace(
            mixed,transitions=((),)+mixed.transitions[1:],
        ))
    with pytest.raises(ValueError):
        construct_reversible_kernel(kernel,cert,blend=Fraction(1))
    with pytest.raises(ValueError):
        construct_reversible_kernel(kernel,cert,blend=Fraction(0))


def test_support_retained_by_strictly_positive_baseline_mixture():
    kernel=build_kernel("HPHPHPHH",Fraction(2))
    cert=propose_maxflow_certificate(kernel,1)
    mixed=construct_reversible_kernel(kernel,cert,blend=Fraction(3,4))
    for i,row in enumerate(kernel.transitions):
        new=dict(mixed.transitions[i])
        for j,p in row:
            if i!=j and p>0:
                assert new[j]>=Fraction(1,4)*p>0
    assert event_flux_reversible(mixed,1)>=event_flux(kernel,1)[1]
