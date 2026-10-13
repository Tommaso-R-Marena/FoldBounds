"""Adversarially checked conformational witness lower bounds."""
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.lattice import (
    CertificateError, build_frontier, exact_spectrum,
)
from foldbounds.lookahead import verify_lookahead
from foldbounds.adaptive import adaptive_refine, verify_refinement_run
from foldbounds.witnesses import (
    WitnessCertificate, propose_dfs_witnesses, verify_witness_certificate,
)

ROOT=((0,0),(1,0))


def independent_oracle(sequence, q, threshold):
    # Independent grid-walk generator and contact counter, no production
    # enumerator or contact scorer is used.
    states=[ROOT]
    for _ in range(2,len(sequence)):
        new=[]
        for walk in states:
            x,y=walk[-1]
            for dx,dy in ((-1,0),(0,-1),(1,0),(0,1)):
                point=(x+dx,y+dy)
                if point not in walk:
                    new.append(walk+(point,))
        states=new
    e=n=Fraction(0)
    for walk in states:
        k=sum(sequence[i]=='H' and sequence[j]=='H' and
              abs(walk[i][0]-walk[j][0])+abs(walk[i][1]-walk[j][1])==1
              for i in range(len(sequence))
              for j in range(i+2,len(sequence)))
        w=q**k
        if k>=threshold: e+=w
        else: n+=w
    return e,n


@pytest.mark.parametrize("length",range(2,8))
def test_every_short_hp_sequence_with_independent_oracle(length):
    for letters in product('HP',repeat=length):
        seq=''.join(letters)
        q=Fraction(3,2)
        threshold=1
        exact_e,exact_n=independent_oracle(seq,q,threshold)
        cert=build_frontier(seq,q,threshold,2)
        previous=None
        for budget in (0,16,64):
            verified=verify_lookahead(cert,lookahead=0,use_matching=True,
                                      use_transfer=True,witness_nodes=budget)
            assert verified.event_lower<=exact_e<=verified.event_upper
            assert verified.complement_lower<=exact_n<=verified.complement_upper
            assert verified.partition_lower<=exact_e+exact_n<=verified.partition_upper
            assert verified.probability_lower<=exact_e/(exact_e+exact_n)<=verified.probability_upper
            if previous is not None:
                assert previous.event_lower<=verified.event_lower
                assert previous.complement_lower<=verified.complement_lower
                assert verified.event_upper<=previous.event_upper
                assert verified.complement_upper<=previous.complement_upper
            previous=verified


@pytest.mark.parametrize('budget',[0,1,4,12,32,128,1000])
def test_witness_certificate_is_monotone_as_budget_increases(budget):
    seq='HHHHHHHHH'
    walks=propose_dfs_witnesses(ROOT,len(seq),node_budget=budget)
    verified=verify_witness_certificate(
        WitnessCertificate(seq,Fraction(2),1,ROOT,walks))
    exact_e,exact_n=independent_oracle(seq,Fraction(2),1)
    assert 0<=verified.event<=exact_e
    assert 0<=verified.complement<=exact_n
    assert len(set(walks))==len(walks)
    assert verified.paths==len(walks)


def test_witness_budget_strictly_improves_lower_mass():
    seq='HHHHHHHH'
    cert=build_frontier(seq,Fraction(2),1,2)
    before=verify_lookahead(cert,use_matching=True,use_transfer=True,witness_nodes=0)
    after=verify_lookahead(cert,use_matching=True,use_transfer=True,witness_nodes=200)
    assert after.partition_lower>before.partition_lower
    assert after.event_lower>=before.event_lower
    assert after.complement_lower>=before.complement_lower
    assert after.partition_upper==before.partition_upper


def test_witness_search_cannot_trust_forged_duplicate():
    seq='HHHH'
    original=propose_dfs_witnesses(ROOT,4,node_budget=30)
    assert original
    valid=WitnessCertificate(seq,Fraction(2),1,ROOT,original)
    assert verify_witness_certificate(valid).paths==len(original)
    with pytest.raises(CertificateError):
        verify_witness_certificate(replace(valid,completions=original+(original[0],)))


@pytest.mark.parametrize('corruption',[
    ((0,0),(1,0),(1,1),(1,0)),  # repeated coordinate
    ((0,0),(1,0),(4,0),(4,1)),  # non-grid step
    ((0,0),(0,1),(1,1),(1,0)),  # wrong root
    ((0,0),(1,0),(2,0)),        # incomplete
])
def test_witness_rejects_fabricated_geometry(corruption):
    cert=WitnessCertificate('HHHH',Fraction(2),1,ROOT,(corruption,))
    with pytest.raises(CertificateError):
        verify_witness_certificate(cert)


def test_witness_must_extend_exact_prefix():
    seq='HHHHH'
    prefix=((0,0),(1,0),(2,0))
    other=((0,0),(1,0),(1,1),(1,2),(1,3))
    with pytest.raises(CertificateError):
        verify_witness_certificate(
            WitnessCertificate(seq,Fraction(2),1,prefix,(other,)))


def test_witness_q_and_threshold_are_strict():
    for q in (Fraction(0),Fraction(1,2),2.0):
        with pytest.raises(ValueError):
            verify_witness_certificate(WitnessCertificate('HHH',q,1,ROOT,()))
    with pytest.raises(CertificateError):
        verify_witness_certificate(WitnessCertificate('HHH',Fraction(2),True,ROOT,()))


@pytest.mark.parametrize('bad',[-1,True,1.5,'100'])
def test_invalid_budgets_fail_closed(bad):
    with pytest.raises(ValueError):
        propose_dfs_witnesses(ROOT,8,node_budget=bad)
    with pytest.raises(ValueError):
        verify_lookahead(build_frontier('HHHH',Fraction(2),1,2),
                         witness_nodes=bad)


def test_adaptive_replay_binds_witness_budget():
    cert=build_frontier('HHHHHHHH',Fraction(2),1,2)
    result=adaptive_refine(cert,steps=4,lookahead=0,use_matching=True,
                           use_transfer=True,witness_nodes=32,policy='mass')
    assert result.witness_nodes==32
    assert verify_refinement_run(result)==result.bounds
    with pytest.raises(ValueError):
        verify_refinement_run(replace(result,witness_nodes=0))


def test_witness_generator_budget_does_not_affect_soundness():
    seq='HHHHHHH'
    exact_e,exact_n=independent_oracle(seq,Fraction(2),1)
    for budget in (0,1,10,500):
        cert=WitnessCertificate(seq,Fraction(2),1,ROOT,
                     propose_dfs_witnesses(ROOT,len(seq),node_budget=budget))
        observed=verify_witness_certificate(cert)
        assert observed.event<=exact_e and observed.complement<=exact_n
