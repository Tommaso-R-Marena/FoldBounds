"""Exhaustive oracles for the optional graph-theoretic contact bound."""
from fractions import Fraction
from itertools import product
import pytest

from foldbounds.lattice import build_frontier, exact_spectrum, _possible_hh_pairs
from foldbounds.lookahead import score_cap, verify_lookahead
from foldbounds.matching import matching_score_cap


def independent_complete_paths(n):
    paths=[((0,0),(1,0))]
    for _ in range(2,n):
        nxt=[]
        for path in paths:
            x,y=path[-1]
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                site=(x+dx,y+dy)
                if site not in path:
                    nxt.append(path+(site,))
        paths=nxt
    return paths


def independent_contacts(seq,path):
    score=0
    for i in range(len(seq)):
        for j in range(i+2,len(seq)):
            if seq[i]!='H' or seq[j]!='H':
                continue
            if abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1:
                score+=1
    return score


@pytest.mark.parametrize('n', range(2,9))
def test_every_short_sequence_and_reachable_prefix(n):
    paths=independent_complete_paths(n)
    for residues in product('HP',repeat=n):
        seq=''.join(residues)
        edges=_possible_hh_pairs(seq)
        best={}
        for path in paths:
            contacts=independent_contacts(seq,path)
            for m in range(2,n):
                prefix=path[:m]
                best[prefix]=max(best.get(prefix,0),contacts)
        for path,max_possible in best.items():
            minimum,upper=matching_score_cap(seq,path,edges)
            old_minimum,old_upper=score_cap(seq,path,edges)
            assert minimum == old_minimum
            assert minimum <= max_possible <= upper <= old_upper


@pytest.mark.parametrize('seq', ['HHHH', 'HHHHHH', 'HHHHHHHH', 'HPHHPHHP', 'HHPHHHPH'])
@pytest.mark.parametrize('q', [Fraction(1), Fraction(3,2), Fraction(2)])
@pytest.mark.parametrize('threshold', [0,1,2,3])
@pytest.mark.parametrize('horizon', [0,2])
def test_matching_encloses_exact_and_narrows_old(seq,q,threshold,horizon):
    spec=exact_spectrum(seq)
    z=sum((count*q**k for k,count in spec.items()),Fraction(0))
    event=sum((count*q**k for k,count in spec.items() if k>=threshold),Fraction(0))
    certificate=build_frontier(seq,q,threshold,2)
    basic=verify_lookahead(certificate,lookahead=horizon)
    matched=verify_lookahead(certificate,lookahead=horizon,use_matching=True)
    assert basic.partition_lower <= matched.partition_lower <= z
    assert z <= matched.partition_upper <= basic.partition_upper
    assert basic.event_lower <= matched.event_lower <= matched.event_upper <= basic.event_upper
    assert basic.complement_lower <= matched.complement_lower <= matched.complement_upper <= basic.complement_upper
    assert basic.probability_lower <= matched.probability_lower <= event/z
    assert event/z <= matched.probability_upper <= basic.probability_upper


def test_matching_strictly_improves_a_constructed_case():
    seq='HHHHHHHH'
    path=((0,0),(1,0),(2,0),(3,0),(4,0))
    edges=_possible_hh_pairs(seq)
    assert matching_score_cap(seq,path,edges)[1] < score_cap(seq,path,edges)[1]


@pytest.mark.parametrize('bad', [0,1,'yes',None])
def test_matching_option_is_explicit_boolean(bad):
    if bad is True or bad is False:
        return
    with pytest.raises(ValueError):
        verify_lookahead(build_frontier('HHHH',Fraction(2),1,2),use_matching=bad)
