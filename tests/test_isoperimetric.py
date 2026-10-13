"""Independent lattice-perimeter identities and bound preservation tests."""
from itertools import product
from fractions import Fraction

import pytest

from foldbounds.isoperimetric import (
    minimum_final_box_sum, prefix_contact_cap, universal_contact_cap,
)
from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead
from foldbounds.matching import matching_score_cap
from foldbounds.lattice import _possible_hh_pairs


def independent_walks(n):
    walks=[((0,0),(1,0))]
    for _ in range(2,n):
        successor=[]
        for path in walks:
            x,y=path[-1]
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                p=(x+dx,y+dy)
                if p not in path:
                    successor.append(path+(p,))
        walks=successor
    return walks


def independent_all_contacts(path):
    occupied=set(path)
    internal_edges=sum((x+1,y) in occupied for x,y in path)
    internal_edges+=sum((x,y+1) in occupied for x,y in path)
    return internal_edges-(len(path)-1)


def independent_edge_perimeter(path):
    cells=set(path)
    return sum((x+dx,y+dy) not in cells for x,y in path
               for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)))


@pytest.mark.parametrize('n',range(2,11))
def test_all_anchored_walks_exact_perimeter_identity_and_prefix_caps(n):
    universal=universal_contact_cap(n)
    for path in independent_walks(n):
        contacts=independent_all_contacts(path)
        perimeter=independent_edge_perimeter(path)
        assert perimeter % 2 == 0
        assert contacts == n+1-perimeter//2
        assert contacts <= universal
        for m in range(2,n+1):
            assert contacts <= prefix_contact_cap(path[:m],n) <= universal


@pytest.mark.parametrize('n',range(2,15))
def test_global_cap_and_matching_no_weaker(n):
    path=tuple((k,0) for k in range(2))
    # A straight anchored prefix of width n-1 permits at most one contact.
    long=tuple((k,0) for k in range(n-1))
    assert prefix_contact_cap(long,n) <= universal_contact_cap(n)
    assert prefix_contact_cap(path,n) == universal_contact_cap(n)
    assert minimum_final_box_sum(n,n-1,1) >= n


@pytest.mark.parametrize('n,cap',[(2,0),(3,0),(4,1),(5,1),(6,2),
                                  (7,2),(8,3),(9,4),(10,4)])
def test_expected_discrete_isoperimetry(n,cap):
    assert universal_contact_cap(n)==cap


@pytest.mark.parametrize('sequence', ['HHHHHHHH','HHPHHHPH','HPPHPHHP','PPPPPPPP','HHHHHHHHHH'])
@pytest.mark.parametrize('threshold',[0,1,2,3])
@pytest.mark.parametrize('q',[Fraction(1),Fraction(3,2),Fraction(2)])
def test_perimeter_bound_tightens_matching_and_preserves_partition(sequence,threshold,q):
    cert=build_frontier(sequence,q,threshold,2)
    without=verify_lookahead(cert,lookahead=1,use_matching=True,use_perimeter=False)
    with_perimeter=verify_lookahead(cert,lookahead=1,use_matching=True,use_perimeter=True)
    spectrum=exact_spectrum(sequence)
    z=sum((count*q**k for k,count in spectrum.items()),Fraction(0))
    e=sum((count*q**k for k,count in spectrum.items() if k>=threshold),Fraction(0))
    assert without.partition_lower <= with_perimeter.partition_lower <= z
    assert z <= with_perimeter.partition_upper <= without.partition_upper
    assert without.probability_lower <= with_perimeter.probability_lower <= e/z
    assert e/z <= with_perimeter.probability_upper <= without.probability_upper


def test_elongated_prefix_reduces_contact_upper_bound():
    path=((0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(6,0))
    n=8
    pairs=_possible_hh_pairs('H'*n)
    no_perimeter=matching_score_cap('H'*n,path,pairs,use_perimeter=False)
    with_perimeter=matching_score_cap('H'*n,path,pairs,use_perimeter=True)
    assert with_perimeter[1]==0
    assert with_perimeter[1]<=no_perimeter[1]


@pytest.mark.parametrize('bad',[0,1,2.0,True,-1])
def test_invalid_global_size(bad):
    with pytest.raises(ValueError):
        universal_contact_cap(bad)


@pytest.mark.parametrize('bad',[0,1.5,'yes',None])
def test_invalid_perimeter_selector(bad):
    cert=build_frontier('HHHHHH',Fraction(2),1,2)
    with pytest.raises(ValueError):
        verify_lookahead(cert,use_matching=True,use_perimeter=bad)
