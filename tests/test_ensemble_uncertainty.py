"""Independent worst-case checks for robust fold-switch population certificates."""
from dataclasses import replace
from fractions import Fraction as F
from random import Random
import pytest
from foldbounds.ensemble_crossover import geometry, spectrum, probability
from foldbounds.ensemble_crossover_3d import spectrum_3d
from foldbounds.ensemble_uncertainty import (
    lower_probability, upper_probability, certify_error_box,
    verify_error_box, certify_switch_band, verify_switch_band,
    critical_factor_ratio_for_majority,
)

SEQ2='HPHHPHHPHHHH';SEQ3='HPHHPHHH'

@pytest.mark.parametrize('dim,seq,t',[ (2,SEQ2,2),(3,SEQ3,1)])
def test_zero_error_recovers_exact_old_parameter_box(dim,seq,t):
    r=certify_error_box(seq,dim,t,F(3,2),F(2),F(1,4),F(4),F(1),F(1))
    assert r.worst_probability_lower==r.nominal_probability_lower
    assert r.worst_probability_upper==r.nominal_probability_upper
    assert verify_error_box(r)

@pytest.mark.parametrize('dim,seq,t,expected_lo,expected_hi',[ 
    (2,SEQ2,2,F(84529143,508499807),F(1108499200,1368704059)),
    (3,SEQ3,1,F(149019,566171),F(27968768,37452129))
])
def test_robust_box_recomputes_exact_corner_probabilities(dim,seq,t,expected_lo,expected_hi):
    lo,hi=F(5,6),F(6,5)
    r=certify_error_box(seq,dim,t,F(3,2),F(2),F(1,4),F(4),lo,hi)
    assert r.nominal_probability_lower==expected_lo
    assert r.nominal_probability_upper==expected_hi
    assert r.worst_probability_lower==lower_probability(expected_lo,lo,hi)
    assert r.worst_probability_upper==upper_probability(expected_hi,lo,hi)
    assert r.status=='undetermined'

@pytest.mark.parametrize('dim,seq,t',[ (2,SEQ2,2),(3,SEQ3,1)])
def test_unique_adversarial_switch_band_has_ordered_exact_sign_changes(dim,seq,t):
    r=certify_switch_band(seq,dim,t,F(3,2),F(1,4),F(4),F(5,6),F(6,5))
    assert verify_switch_band(r)
    lo,mid,hi=(sum(pair)/2 for pair in (r.guaranteed_complement_until,
                 r.nominal_crossover,r.guaranteed_event_from))
    assert lo<mid<hi
    assert all(x[1]-x[0]<F(1,10**12) for x in
               (r.guaranteed_complement_until,r.nominal_crossover,r.guaranteed_event_from))
    h=spectrum(seq) if dim==2 else spectrum_3d(seq)
    b1=r.guaranteed_complement_until[0]-F(1,10000)
    b2=r.guaranteed_event_from[1]+F(1,10000)
    p1=probability(h,F(3,2),b1,t)
    p2=probability(h,F(3,2),b2,t)
    assert upper_probability(p1,F(5,6),F(6,5))<F(1,2)
    assert lower_probability(p2,F(5,6),F(6,5))>F(1,2)

@pytest.mark.parametrize('dim,seq,t',[ (2,SEQ2,2),(3,SEQ3,1)])
def test_nominal_zero_width_switch_band(dim,seq,t):
    r=certify_switch_band(seq,dim,t,F(3,2),F(1,4),F(4),F(1),F(1))
    assert r.guaranteed_complement_until==r.guaranteed_event_from==r.nominal_crossover


def test_adversarial_extrema_are_attained_by_allowed_per_state_corrections():
    # Independent from closed-form ratio formula: directly enumerate all
    # 30,073 walks and reweight each. Extreme assignments are actually legal.
    seq=SEQ2; q=F(3,2); b=F(1); t=2
    lo,hi=F(5,6),F(6,5)
    EL=CL=EU=CU=F(0)
    for geom in geometry(len(seq)):
        k=sum(seq[i]=='H' and seq[j]=='H' for i,j in geom.contacts)
        weight=q**k*b**geom.bends
        if k>=t:
            EL+=weight*lo;EU+=weight*hi
        else:
            CL+=weight*hi;CU+=weight*lo
    p=probability(spectrum(seq),q,b,t)
    assert EL/(EL+CL)==lower_probability(p,lo,hi)
    assert EU/(EU+CU)==upper_probability(p,lo,hi)
    assert 0<EL/(EL+CL)<p<EU/(EU+CU)<1


def test_sampled_uncertainties_within_certified_bounds_2d_and_3d():
    rng=Random(4179)
    for dim,seq,t in ((2,SEQ2,2),(3,SEQ3,1)):
        h=spectrum(seq) if dim==2 else spectrum_3d(seq)
        r=certify_error_box(seq,dim,t,F(3,2),F(2),F(1,4),F(4),F(4,5),F(5,4))
        keys=list(h.items())
        # Group-dependent model errors are sampled adversarially without
        # independence assumptions. A single microstate has its own factor;
        # groupwise factors are a subset of allowed statewise assignments.
        for i in range(60):
            q=F(3,2)+F(rng.randrange(101),200)
            b=F(1,4)+F(15*rng.randrange(101),400)
            e=c=F(0)
            for (k,bends),count in keys:
                f=F(4,5)+F(9*rng.randrange(101),2000)
                w=count*q**k*b**bends*f
                if k>=t:e+=w
                else:c+=w
            p=e/(e+c)
            assert r.worst_probability_lower<=p<=r.worst_probability_upper


def test_tamper_rejection_and_failclosed_unsupported_model():
    r=certify_error_box(SEQ2,2,2,F(3,2),F(2),F(1,4),F(4),F(5,6),F(6,5))
    assert not verify_error_box(replace(r,worst_probability_lower=F(0)))
    assert not verify_error_box(replace(r,derivative_terms=r.derivative_terms+1))
    assert not verify_error_box(replace(r,dimension=4))
    assert not verify_error_box(object())
    switch=certify_switch_band(SEQ2,2,2,F(3,2),F(1,4),F(4),F(5,6),F(6,5))
    assert not verify_switch_band(replace(switch,guaranteed_event_from=(F(0),F(1))))
    assert not verify_switch_band(object())
    with pytest.raises(ValueError):certify_error_box('HHHHPHHP',2,2,F(3,2),F(2),F(1,4),F(4),F(1),F(1))
    with pytest.raises(TypeError):certify_error_box(SEQ2,2,2,1.5,F(2),F(1,4),F(4),F(1),F(1))
    with pytest.raises(ValueError):certify_switch_band(SEQ3,3,1,F(3,2),F(1,4),F(4),F(1,2),F(2))


def test_critical_error_odds_and_global_factor_cancellation():
    assert critical_factor_ratio_for_majority(F(9,10))==9
    assert critical_factor_ratio_for_majority(F(3,5))==F(3,2)
    for p in (F(1,10),F(1,2),F(9,10)):
        for f in (F(1,3),F(3),F(111,100)):
            assert lower_probability(p,f,f)==p
            assert upper_probability(p,f,f)==p
    with pytest.raises(ValueError):critical_factor_ratio_for_majority(F(1,2))

@pytest.mark.parametrize('ratio',[F(11,10),F(6,5),F(3,2)])
def test_wider_energy_uncertainty_never_narrows_switch_band_2d(ratio):
    r=certify_switch_band(SEQ2,2,2,F(3,2),F(1,4),F(4),F(1)/ratio,ratio)
    assert verify_switch_band(r)
    narrow=certify_switch_band(SEQ2,2,2,F(3,2),F(1,4),F(4),F(10,11),F(11,10))
    assert r.guaranteed_complement_until[0]<=narrow.guaranteed_complement_until[0]
    assert r.guaranteed_event_from[1]>=narrow.guaranteed_event_from[1]

from foldbounds.ensemble_uncertainty import (
    CertifiedMassIntervals, reweight_certified_mass_intervals)

@pytest.mark.parametrize('lo,hi,lc,uc',[ 
    (F(0),F(10),F(0),F(20)),
    (F(2),F(12),F(3),F(17)),
    (F(0),F(1),F(10),F(50)),
    (F(0),F(0),F(1),F(2)),
    (F(1),F(2),F(0),F(0)),
])
def test_interval_mass_reweighting_contains_all_bounded_exact_masses(lo,hi,lc,uc):
    m=CertifiedMassIntervals(lo,hi,lc,uc)
    b=reweight_certified_mass_intervals(m,F(5,6),F(6,5))
    for e in (lo,hi,(lo+hi)/2):
        for c in (lc,uc,(lc+uc)/2):
            if e+c==0:continue
            for fe in (F(5,6),F(6,5)):
                for fc in (F(5,6),F(6,5)):
                    true_p=e*fe/(e*fe+c*fc)
                    assert b.population_lower<=true_p<=b.population_upper


def test_mass_reweighting_rational_sharpness_and_invalid_input_rejection():
    m=CertifiedMassIntervals(F(5),F(9),F(7),F(13))
    x=reweight_certified_mass_intervals(m,F(5,6),F(6,5))
    assert x.population_lower==F(5)*F(5,6)/(F(5)*F(5,6)+F(13)*F(6,5))
    assert x.population_upper==F(9)*F(6,5)/(F(9)*F(6,5)+F(7)*F(5,6))
    for bad in (CertifiedMassIntervals(F(-1),F(10),F(1),F(2)),
                CertifiedMassIntervals(F(2),F(1),F(1),F(2)),
                CertifiedMassIntervals(F(0),F(0),F(0),F(0))):
        with pytest.raises(ValueError):reweight_certified_mass_intervals(bad,F(1),F(2))


def test_variable_bit_length_witness_verifies_at_32_and_64():
    a=certify_switch_band(SEQ2,2,2,F(3,2),F(1,4),F(4),F(5,6),F(6,5),bits=32)
    b=certify_switch_band(SEQ2,2,2,F(3,2),F(1,4),F(4),F(5,6),F(6,5),bits=64)
    assert verify_switch_band(a) and verify_switch_band(b)
    assert b.nominal_crossover[1]-b.nominal_crossover[0] < a.nominal_crossover[1]-a.nominal_crossover[0]