"""Independent checks of reference-to-mutant exact reweighting & thermodynamic cycles."""
from fractions import Fraction as F
from dataclasses import replace
from itertools import combinations
import math
import pytest
from foldbounds.counterfactual_reweighting import (
    reweight_cycle,critical_energy_error_kcal,population_with_factor_error)
from foldbounds.mutation_switch import nominal_probability,certify_mutation_switch
from foldbounds.mutation_switch_10 import population10,certify_ten_mutation

@pytest.mark.parametrize('seq,pairs',[
    ('PHPHPHPH',(1,3)),('PPHPPHPP',(1,4)),
    ('HHPHPHPH',(2,6)),('PHPHPHPHPH',(1,3)),
    ('HHPPHPHPPH',(2,5))])
def test_exact_reference_reweighting_agrees_with_all_four_direct_partitions(seq,pairs):
    n=len(seq); q=F(8);b=F(1)
    r=reweight_cycle(seq,*pairs,q,b)
    basis=[seq]
    for sites in [(pairs[0],),(pairs[1],),pairs]:
        a=list(seq)
        for i in sites:a[i-1]='H' if a[i-1]=='P' else 'P'
        basis.append(''.join(a))
    direct=(lambda sequence:nominal_probability(sum((c=='H')<<i for i,c in enumerate(sequence)),q,b)) if n==8 else (lambda sequence:population10(sequence,q,b))
    values=[direct(x) for x in basis]
    assert values==[r.baseline_event_probability,r.single_i_event_probability,r.single_j_event_probability,r.double_event_probability]
    if (pairs[0]-pairs[1])%2==0:
        assert not r.direct_contact_possible
        assert r.odds_epistasis_factor==r.event_interaction_factor/r.complement_interaction_factor
    assert r.enumerated_microstates==(13565 if n==8 else 308981)

@pytest.mark.parametrize('n,seq',[(8,'PHPHPHPH'),(10,'PHPHPHPHPH')])
def test_no_direct_interaction_but_ensemble_covariance_generates_epistasis(n,seq):
    x=reweight_cycle(seq,1,3,F(8),F(1))
    assert not x.direct_contact_possible
    assert x.odds_epistasis_factor != 1
    assert x.event_conditional_covariance != 0
    assert x.complement_conditional_covariance != 0
    assert x.odds_epistasis_factor>0

@pytest.mark.parametrize('seq,sitepair',[('PHPHPHPH',(1,4)),('PHPHPHPH',(2,5)),('PHPHPHPHPH',(1,4))])
def test_adjacent_parity_may_have_direct_contact_and_does_not_assume_factorization(seq,sitepair):
    x=reweight_cycle(seq,*sitepair,F(3,2),F(1))
    assert x.direct_contact_possible

@pytest.mark.parametrize('bad',[(0,3),(1,11),(2,2),(3,2),(True,3),(1,False)])
def test_malformed_mutation_sites_fail_closed(bad):
    with pytest.raises(ValueError):reweight_cycle('PHPHPHPH',*bad,F(8),F(1))

@pytest.mark.parametrize('p,lo,hi',[(F(1,3),F(1,2),F(2)),(F(1,2),F(1),F(1)),(F(0),F(1,2),F(3)),(F(1),F(2,3),F(3,2))])
def test_adversarial_reweight_envelope_exact_and_attainable(p,lo,hi):
    lower,upper=population_with_factor_error(p,lo,hi)
    assert lower<=p<=upper
    assert lower==p*lo/(p*lo+(1-p)*hi)
    assert upper==p*hi/(p*hi+(1-p)*lo)


def test_room_temperature_energy_error_no_go_for_strong_guarantees():
    c8=certify_mutation_switch('PHPHPHPH','HHPHPHPH',F(8),F(12),F(3,4),F(5,4))
    c10=certify_ten_mutation()
    t8=critical_energy_error_kcal(c8.parent.event_upper,c8.mutant.event_lower)
    t10=critical_energy_error_kcal(c10.original_nominal_upper,c10.mutant_nominal_lower)
    assert .228<t8<.231 and .111<t10<.113
    assert t8>t10
    # Weight factors 5/6 .. 6/5 correspond to just ~0.106 kcal/mol
    energy=0.00198720425864083*293.15*math.log(1.2)
    assert .105<energy<.107 and energy<t10
    # Already a factor-of-2 uncertainty destroys the guaranteed majority
    # despite same finite-model nominal populations.
    for old,mut in [(c8.parent.event_upper,c8.mutant.event_lower),
                    (c10.original_nominal_upper,c10.mutant_nominal_lower)]:
        lower_old,upper_old=population_with_factor_error(old,F(1,2),F(2))
        lower_mut,upper_mut=population_with_factor_error(mut,F(1,2),F(2))
        assert not (upper_old<F(1,2)<lower_mut)


def test_ga_gb_population_thermodynamic_minimum_from_published_90_percent():
    # p>0.9 vs p<0.1 imply delta G difference exceeds 2RT ln 9.
    RT=0.00198720425864083*(273.15+20)
    assert 2.55<2*RT*math.log(9)<2.57

def test_full_768_sublattice_cycle_falsification():
    # All possible P->H two-site mutations at parity-matched indices of 8-bead chains.
    # Their mutated sites can never form a nonconsecutive cubic-lattice contact.
    count=0;nonzero=0;pos=neg=0
    for i,j in combinations(range(1,9),2):
        if (i-j)%2:continue
        sites=[s for s in range(8) if s not in (i-1,j-1)]
        for othermask in range(1<<6):
            bits=sum((1<<s) for idx,s in enumerate(sites) if othermask&(1<<idx))
            seq=''.join('H' if bits&(1<<k) else 'P' for k in range(8))
            r=reweight_cycle(seq,i,j,F(8),F(1))
            assert not r.direct_contact_possible
            if r.odds_epistasis_factor!=1:
                nonzero+=1
                pos+=r.odds_epistasis_factor>1
                neg+=r.odds_epistasis_factor<1
            count+=1
    assert (count,nonzero,pos,neg)==(768,504,140,364)