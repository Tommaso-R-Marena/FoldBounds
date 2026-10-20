"""Independent, falsification-first tests for exact toy lattice mutation transitions."""
from fractions import Fraction as F
from collections import Counter
from dataclasses import replace
import pytest
from foldbounds.mutation_switch import (
    _sequence,_mask,compact_spectrum,nominal_probability,certify_population_box,
    verify_population_box,certify_mutation_switch,verify_mutation_switch,
    enumerate_mutation_landscape)
from foldbounds.ensemble_uncertainty import lower_probability,upper_probability

PARENT='PHPHPHPH'
MUTANT='HHPHPHPH'


def _independent_dfs_spectra():
    """Separate enumerator, independent of cached geometry_3d and production scoring."""
    S=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
    path=[(0,0,0),(1,0,0)]; visited=set(path)
    out={170:Counter(),171:Counter()}
    def rec(last,turns):
        if len(path)==8:
            pairs=[]
            for i in range(8):
                for j in range(i+2,8):
                    if sum(abs(a-b) for a,b in zip(path[i],path[j]))==1:
                        pairs.append((i,j))
            for mask in out:
                h=sum(bool((mask>>i)&1) and bool((mask>>j)&1) for i,j in pairs)
                out[mask][(h,turns,int(len(pairs)>=2))]+=1
            return
        for dx,dy,dz in S:
            nx=(path[-1][0]+dx,path[-1][1]+dy,path[-1][2]+dz)
            if nx not in visited:
                visited.add(nx);path.append(nx)
                rec((dx,dy,dz),turns+((dx,dy,dz)!=last))
                path.pop();visited.remove(nx)
    rec((1,0,0),0)
    return out


def test_independent_full_3d_histograms():
    independent=_independent_dfs_spectra()
    for m in (170,171):
        assert dict(compact_spectrum(m))==dict(independent[m])
        assert sum(independent[m].values())==13565
        assert sum(v for (k,r,event),v in independent[m].items() if event)==2760


def test_single_mutation_is_geometric_not_target_label_change():
    assert _mask(PARENT)==170 and _mask(MUTANT)==171
    assert all(k==0 for (k,r,event),v in compact_spectrum(170))
    assert len([1 for (k,r,event),v in compact_spectrum(171) if k>0])>0


def test_exact_continuum_switch_worst_case_independent_factors():
    c=certify_mutation_switch(PARENT,MUTANT,F(8),F(12),F(3,4),F(5,4),F(5,6),F(6,5))
    assert verify_mutation_switch(c)
    assert c.mutated_residue==1
    assert c.parent.event_lower==nominal_probability(170,F(8),F(3,4))
    assert c.parent.event_upper==nominal_probability(170,F(12),F(5,4))
    assert c.mutant.event_lower==nominal_probability(171,F(8),F(3,4))
    assert c.mutant.event_upper==nominal_probability(171,F(12),F(5,4))
    assert c.parent_robust_upper < F(3,10) < F(1,2)
    assert c.mutant_robust_lower > F(3,5) > F(1,2)
    assert c.parent_robust_upper<F(1,2)<c.mutant_robust_lower


def test_bivariate_bernstein_includes_independent_dense_rational_grid():
    for seq in (PARENT,MUTANT):
        c=certify_population_box(seq,F(8),F(12),F(3,4),F(5,4))
        assert verify_population_box(c)
        for qi in range(9):
            for bi in range(9):
                q=F(8)+F(qi,2);b=F(3,4)+F(bi,16)
                p=nominal_probability(_mask(seq),q,b)
                assert c.event_lower<=p<=c.event_upper
        # Rational Bernstein coefficient convex-hull is tight in this case.
        assert c.event_lower==nominal_probability(_mask(seq),F(8),F(3,4))
        assert c.event_upper==nominal_probability(_mask(seq),F(12),F(5,4))


def test_explicit_statewise_adversarial_reweighting_attains_extrema():
    for m in (170,171):
        for q,b in ((F(8),F(3,4)),(F(12),F(5,4))):
            h=dict(compact_spectrum(m))
            e=sum((v*q**k*b**r for (k,r,compact),v in h.items() if compact),F())
            c=sum((v*q**k*b**r for (k,r,compact),v in h.items() if not compact),F())
            p=e/(e+c)
            lo=(F(5,6)*e)/(F(5,6)*e+F(6,5)*c)
            hi=(F(6,5)*e)/(F(6,5)*e+F(5,6)*c)
            assert lo==lower_probability(p,F(5,6),F(6,5))
            assert hi==upper_probability(p,F(5,6),F(6,5))


def test_all_1024_single_mutation_edges_classified_without_selection():
    # Complete binary sequence space and all 8*2^7 one-site mutation edges.
    assert len(tuple((m,m^(1<<i)) for m in range(256) for i in range(8) if m < (m^(1<<i))))==1024
    assert len(enumerate_mutation_landscape(6))==40
    assert len(enumerate_mutation_landscape(8))==92
    assert len(enumerate_mutation_landscape(12))==134
    for a,z,i in enumerate_mutation_landscape(8):
        assert (a^z)==1<<(i-1)
        assert upper_probability(nominal_probability(a,F(8),F(1)),F(5,6),F(6,5))<F(1,2)
        assert lower_probability(nominal_probability(z,F(8),F(1)),F(5,6),F(6,5))>F(1,2)


@pytest.mark.parametrize('mode',[0,1,2,3,4,5])
def test_forged_certificates_are_rejected(mode):
    obj=certify_mutation_switch(PARENT,MUTANT,F(8),F(12),F(3,4),F(5,4))
    mutations=(replace(obj,parent_robust_upper=F(0)),
               replace(obj,mutant_robust_lower=F(1)),
               replace(obj,mutated_residue=2),
               replace(obj,orientation='high-to-low'),
               replace(obj,parent=replace(obj.parent,event_upper=F(0))),
               replace(obj,mutant=replace(obj.mutant,sequence=PARENT)))
    assert not verify_mutation_switch(mutations[mode])


@pytest.mark.parametrize('seq', ['', 'HPHP', 'Z'*8, 'h'*8, 'H'*9, None])
def test_reject_invalid_sequences(seq):
    with pytest.raises((ValueError,TypeError)):
        certify_population_box(seq,F(8),F(12),F(3,4),F(5,4))


def test_reject_non_single_mutation_or_unproved_switch():
    with pytest.raises(ValueError):
        certify_mutation_switch(PARENT,PARENT,F(8),F(12),F(3,4),F(5,4))
    with pytest.raises(ValueError):
        certify_mutation_switch(PARENT,'HHPHHPPH',F(8),F(12),F(3,4),F(5,4))
    with pytest.raises(ValueError):
        certify_mutation_switch(PARENT,MUTANT,F(2),F(3),F(3,4),F(5,4))


def test_no_dishonest_certification_under_large_unknown_energy_error():
    with pytest.raises(ValueError):
        certify_mutation_switch(PARENT,MUTANT,F(8),F(12),F(3,4),F(5,4),F(1,2),F(2))


def test_refinement_never_widens_population_enclosure():
    for seq in (PARENT,MUTANT):
        a=certify_population_box(seq,F(8),F(12),F(3,4),F(5,4),0)
        b=certify_population_box(seq,F(8),F(12),F(3,4),F(5,4),2)
        assert a.event_lower<=b.event_lower<=b.event_upper<=a.event_upper


def test_entire_uniform_parameter_box_landscape_recomputed():
    from foldbounds.mutation_switch import certify_uniform_mutation_landscape,verify_uniform_mutation_landscape
    c=certify_uniform_mutation_landscape(F(8),F(12),F(3,4),F(5,4))
    assert verify_uniform_mutation_landscape(c)
    assert c.checked_edges==1024
    assert c.microstates_per_sequence==13565
    assert len(c.robust_switch_edges)==36
    assert (170,171,1) in c.robust_switch_edges
    # Every P->H switch has exactly one mutated site and retains a fixed geometry event.
    for lo,hi,site in c.robust_switch_edges:
        assert lo^hi == 1<<(site-1)
        assert not (lo & (1<<(site-1))) and bool(hi & (1<<(site-1)))


def test_landscape_forgery_rejected():
    from foldbounds.mutation_switch import certify_uniform_mutation_landscape,verify_uniform_mutation_landscape
    c=certify_uniform_mutation_landscape(F(8),F(12),F(3,4),F(5,4))
    assert not verify_uniform_mutation_landscape(replace(c,checked_edges=1023))
    assert not verify_uniform_mutation_landscape(replace(c,robust_switch_edges=((170,171,1),)))
    assert not verify_uniform_mutation_landscape(replace(c,factor_upper=F(2)))


def test_conformational_epistasis_without_possible_direct_mutant_contact():
    # First, mutant sites 4 and 6 differ by two covalent bonds, so in a cubic
    # bipartite lattice they CANNOT be adjacent. Nonetheless equilibrium
    # geometry-event odds exhibit strong population-mediated epistasis.
    m=1; i=3; j=5
    a=m|(1<<i);b=m|(1<<j);ab=a|(1<<j)
    p=[nominal_probability(x,F(8),F(1)) for x in (m,a,b,ab)]
    assert p[0] == F(552,2713)
    assert max(p[:3])<F(1,2)<p[3]
    # Exact odds epistasis cannot result from an i-j contact energy term;
    # mutations at even sequence separation cannot directly lattice-contact.
    odds=lambda z:z/(1-z)
    assert odds(p[3])*odds(p[0])!=odds(p[1])*odds(p[2])
    assert all(upper_probability(x,F(10,11),F(11,10))<F(1,2) for x in p[:3])
    assert lower_probability(p[3],F(10,11),F(11,10))>F(1,2)


def test_full_experiment_determinism_and_reversal_invariance():
    from foldbounds.mutation_switch import geometric_fingerprints
    assert len(geometric_fingerprints())==171
    for m in range(256):
        rev=int(f'{m:08b}'[::-1],2)
        assert compact_spectrum(m)==compact_spectrum(rev)
    from pathlib import Path
    from runpy import run_path
    run=run_path(str(Path(__file__).resolve().parents[1]/"experiments"/"certified_mutation_screen.py"))["run"]
    data=run()
    assert data['robust_switch_reversal_orbits']==18
    assert data['double_mutation_epistasis']['nonzero_exact_odds_epistasis_despite_no_direct_contact']==504
    assert data['double_mutation_epistasis']['example_parent_single_single_double'][0]=='HPPPPPPP'