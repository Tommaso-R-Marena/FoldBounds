"""Independent exact polynomial and explicit ensemble oracles for Bernstein bounds."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
from random import Random
import pytest
from foldbounds.ensemble_crossover import probability,spectrum
from foldbounds.ensemble_crossover_3d import spectrum_3d
from foldbounds.ensemble_uncertainty import upper_probability,lower_probability
from foldbounds.bernstein_population import (
    _monomial_to_bernstein,_bernstein_ratio_bounds,
    certify_bernstein_population,verify_bernstein_population,
    certified_bernstein_with_errors,
)

def bernstein(n,k,t):
    from math import comb
    return F(comb(n,k))*t**k*(1-t)**(n-k)

@pytest.mark.parametrize('exp,degree',[(0,0),(0,5),(1,5),(3,5),(5,5)])
def test_exact_bernstein_basis_matches_independent_power_evaluations(exp,degree):
    lo=F(1,4);hi=F(7,2)
    coeffs=_monomial_to_bernstein(exp,degree,lo,hi)
    for p in (F(0),F(1,12),F(1,4),F(1,2),F(3,4),F(1)):
        val=sum((c*bernstein(degree,k,p) for k,c in enumerate(coeffs)),F(0))
        assert val==(lo+(hi-lo)*p)**exp

@pytest.mark.parametrize('dim,seq,t',[ (2,'HHHPPHHH',1),(3,'HHPPPPPH',1),
    (2,'HPHHPHHPHHHH',2),(3,'HPHHPHHH',1)])
def test_certified_nonmonotone_full_rectangle_encloses_exact_grid(dim,seq,t):
    h=spectrum(seq) if dim==2 else spectrum_3d(seq)
    prior=None
    for depth in range(4):
        cert=certify_bernstein_population(seq,dim,t,F(3,2),F(2),F(1,4),F(8),depth)
        assert verify_bernstein_population(cert)
        if prior is not None:
            assert prior.lower_probability<=cert.lower_probability
            assert cert.upper_probability<=prior.upper_probability
        prior=cert
        for i in range(13):
            q=F(3,2)+F(i,24)
            for j in range(13):
                b=F(1,4)+F(31*j,48)
                p=probability(h,q,b,t)
                assert cert.lower_probability<=p<=cert.upper_probability


def test_bending_counterexamples_are_actual_not_failed_sufficient_witnesses_only():
    # Exact fractions prove two population decreases. Thus universal
    # bend-monotonicity cannot be claimed for arbitrary sequences.
    h2=spectrum('HHHPPHHH')
    assert probability(h2,F(3,2),F(4),1)>probability(h2,F(3,2),F(8),1)
    h3=spectrum_3d('HHPPPPPH')
    assert probability(h3,F(3,2),F(4),1)>probability(h3,F(3,2),F(8),1)
    from foldbounds.ensemble_crossover import robust_bend_sign_certificate
    for h in (h2,h3):
        _,negative,_=robust_bend_sign_certificate(h,1,F(3,2),F(1,4))
        assert negative


def test_bernstein_encloses_parameter_induced_nonmonotone_response_and_errors():
    seq='HHPPPPPH'
    cert=certify_bernstein_population(seq,3,1,F(3,2),F(2),F(1,4),F(8),depth=3)
    lo,hi=certified_bernstein_with_errors(cert,F(5,6),F(6,5))
    assert lo<=cert.lower_probability
    assert hi>=cert.upper_probability
    assert hi<F(1,2)  # Definitely no majority contact event across continuum + error.
    hist=spectrum_3d(seq)
    rng=Random(98744)
    for _ in range(200):
        q=F(3,2)+F(rng.randrange(101),200)
        b=F(1,4)+F(31*rng.randrange(101),400)
        nominal=probability(hist,q,b,1)
        assert lo<=lower_probability(nominal,F(5,6),F(6,5))
        assert upper_probability(nominal,F(5,6),F(6,5))<=hi


def test_fail_closed_and_tamper_rejection():
    cert=certify_bernstein_population('HPHHPHHH',3,1,F(3,2),F(2),F(1,4),F(4),2)
    assert not verify_bernstein_population(replace(cert,lower_probability=F(0)))
    assert not verify_bernstein_population(replace(cert,microstates=0))
    assert not verify_bernstein_population(replace(cert,dimension=7))
    assert not verify_bernstein_population(object())
    with pytest.raises(ValueError):certified_bernstein_with_errors(replace(cert,subboxes=0),F(1),F(2))
    with pytest.raises(TypeError):certify_bernstein_population('HPHHPHHH',3,1,1.5,F(2),F(1,4),F(4))
    with pytest.raises(ValueError):certify_bernstein_population('HPHHPHHH',3,1,F(3,2),F(2),F(1,4),F(4),depth=6)


def test_exact_reductions_for_degenerate_parameter_rectangles():
    for dim,seq,t in ((2,'HHHPPHHH',1),(3,'HHPPPPPH',1)):
        q,b=F(3,2),F(3,4)
        cert=certify_bernstein_population(seq,dim,t,q,q,b,b,depth=0)
        h=spectrum(seq) if dim==2 else spectrum_3d(seq)
        exact=probability(h,q,b,t)
        assert cert.lower_probability==cert.upper_probability==exact


def test_nonmonotone_certificates_across_32_stratified_eight_mer_2d_sequences():
    # The separate sequence-falsification experiment scans all 256; this
    # Bernstein CI regression uses a deterministic 32-sequence subset.
    for i,word in enumerate(product('HP',repeat=8)):
        if i % 8:continue
        seq=''.join(word)
        a=certify_bernstein_population(seq,2,1,F(3,2),F(2),F(1,4),F(8),depth=0)
        b=certify_bernstein_population(seq,2,1,F(3,2),F(2),F(1,4),F(8),depth=1)
        assert a.lower_probability<=b.lower_probability<=b.upper_probability<=a.upper_probability