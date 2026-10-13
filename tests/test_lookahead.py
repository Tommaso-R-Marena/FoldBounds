"""Exhaustive independent oracles and adversarial certificates."""
from dataclasses import replace
from fractions import Fraction
from itertools import product
import pytest
from foldbounds.lattice import (CertificateError,Node,ROOT,build_frontier,
                                exact_spectrum,verify_certificate,refine_frontier)
from foldbounds.lookahead import verify_lookahead,score_cap,walk_count_upper
from foldbounds.lattice import _possible_hh_pairs,_neighbors


def oracle(seq,q,event):
    density=exact_spectrum(seq)
    z=sum((count*q**k for k,count in density.items()),Fraction(0))
    a=sum((count*q**k for k,count in density.items() if k>=event),Fraction(0))
    return z,a/z


@pytest.mark.parametrize('n',range(2,9))
@pytest.mark.parametrize('q',[Fraction(1),Fraction(3,2),Fraction(3)])
@pytest.mark.parametrize('threshold',[0,1,2])
def test_all_hp_sequences_up_to_eight(n,q,threshold):
    for seq_tuple in product('HP',repeat=n):
        seq=''.join(seq_tuple)
        z,p=oracle(seq,q,threshold)
        cert=build_frontier(seq,q,threshold,2)
        a=verify_lookahead(cert,lookahead=0)
        b=verify_lookahead(cert,lookahead=1)
        assert a.partition_lower <= z <= a.partition_upper
        assert a.probability_lower <= p <= a.probability_upper
        assert a.event_lower <= b.event_lower <= b.event_upper <= a.event_upper
        assert a.complement_lower <= b.complement_lower <= b.complement_upper <= a.complement_upper
        assert b.partition_lower <= z <= b.partition_upper
        assert b.probability_lower <= p <= b.probability_upper


@pytest.mark.parametrize('sequence',['HHHHHH','HPPHPHH','HHHHHHHH','PPPPPPPP'])
@pytest.mark.parametrize('q',[Fraction(1),Fraction(5,2)])
@pytest.mark.parametrize('event',[0,1,2,3])
def test_hierarchy_narrows_to_exact(sequence,q,event):
    cert=build_frontier(sequence,q,event,2)
    z,p=oracle(sequence,q,event)
    levels=[verify_lookahead(cert,lookahead=i) for i in range(len(sequence)-1)]
    for old,new in zip(levels,levels[1:]):
        assert old.event_lower <= new.event_lower <= new.event_upper <= old.event_upper
        assert old.complement_lower <= new.complement_lower <= new.complement_upper <= old.complement_upper
        assert old.partition_lower <= new.partition_lower <= z <= new.partition_upper <= old.partition_upper
        assert old.probability_lower <= new.probability_lower <= p <= new.probability_upper <= old.probability_upper
    assert levels[-1].exact
    assert levels[-1].partition_lower == z
    assert levels[-1].probability_lower == p == levels[-1].probability_upper


def test_bounded_contact_cap():
    seq='HHHHHHHHHH'; pairs=_possible_hh_pairs(seq)
    for path in _neighbors(ROOT):
        lo,hi=score_cap(seq,path,pairs)
        assert 0 <= lo <= hi <= len(pairs)


def test_refinement_preserves_enclosures():
    cert=build_frontier('HHHHHH',Fraction(2),1,3)
    new=refine_frontier(cert,cert.frontier[0].path)
    z,p=oracle(cert.sequence,cert.q,cert.event_at_least)
    for c in (cert,new):
        result=verify_lookahead(c,lookahead=1)
        assert result.partition_lower <= z <= result.partition_upper
        assert result.probability_lower <= p <= result.probability_upper


def test_tampered_certificates_rejected():
    c=build_frontier('HHHHHHH',Fraction(2),1,3)
    for corrupt in (replace(c,frontier=c.frontier[:-1]),
                    replace(c,frontier=c.frontier+c.frontier[:1]),
                    replace(c,frontier=c.frontier+(Node(ROOT,'bound'),)),
                    replace(c,frontier=(Node(((0,0),(1,0),(0,0)),'bound'),))):
        with pytest.raises(CertificateError):
            verify_lookahead(corrupt,lookahead=2)


def test_trivial_events_certified_immediately():
    c=build_frontier('HHHHHH',Fraction(2),-1,2)
    assert verify_lookahead(c).probability_lower == verify_lookahead(c).probability_upper == 1
    c=replace(c,event_at_least=99)
    assert verify_lookahead(c).probability_lower == verify_lookahead(c).probability_upper == 0


def test_nonzero_witnesses_and_stronger_upper():
    c=build_frontier('HHHHHHHH',Fraction(2),1,2)
    old=verify_certificate(c)
    new=verify_lookahead(c,lookahead=3)
    z,p=oracle(c.sequence,c.q,c.event_at_least)
    assert new.partition_lower > old.partition_lower
    assert new.partition_upper < old.partition_upper
    assert new.partition_lower <= z <= new.partition_upper
    assert new.probability_lower <= p <= new.probability_upper


@pytest.mark.parametrize('remaining',range(1,13))
def test_universal_walk_upper_is_below_simple_branch_count(remaining):
    value=walk_count_upper(remaining)
    assert value <= 3 ** remaining
    if remaining >= 4:
        assert value < 3 ** remaining


def test_universal_walk_bounds_real_prefix_subtrees():
    # An independent DFS counts true completions for many trapped prefixes.
    for n in range(2,10):
        layers={ROOT:1}
        prefix_counts={}
        def enumerate_paths(path):
            if len(path)==n:
                return 1
            count=sum(enumerate_paths(child) for child in _neighbors(path))
            prefix_counts[path]=count
            return count
        enumerate_paths(ROOT)
        for path,count in prefix_counts.items():
            assert count <= walk_count_upper(n-len(path))
