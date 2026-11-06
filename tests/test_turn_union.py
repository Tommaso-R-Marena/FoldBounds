"""Independent exact-oracle union of monotone and <=2 turn suffix paths."""
from fractions import Fraction
import pytest
from foldbounds.lattice import ROOT,_neighbors
from foldbounds.turn_union import monotone_turn_union
from foldbounds.monotone import monotone_family_mass

CONES=(frozenset(((1,0),(0,1))),frozenset(((1,0),(0,-1))),
       frozenset(((-1,0),(0,1))),frozenset(((-1,0),(0,-1))))

def exact_walks(n):
    paths=[ROOT]
    for _ in range(2,n):
        paths=[child for path in paths for child in _neighbors(path)]
    return paths

def membership(full,prefix):
    if full[:len(prefix)]!=prefix:return False
    moves=tuple((full[j][0]-full[j-1][0],full[j][1]-full[j-1][1])
                for j in range(len(prefix),len(full)))
    return (any(set(moves)<=cone for cone in CONES)
            or sum(a!=b for a,b in zip(moves,moves[1:]))<=2)

def oracle(seq,prefix,q,t):
    e=c=Fraction(0);count=0
    for full in exact_walks(len(seq)):
        if not membership(full,prefix):continue
        score=sum(seq[i]=='H' and seq[j]=='H'
                 and abs(full[i][0]-full[j][0])+abs(full[i][1]-full[j][1])==1
                 for i in range(len(seq)) for j in range(i+2,len(seq)))
        weight=q**score
        if score>=t:e+=weight
        else:c+=weight
        count+=1
    return e,c,count

@pytest.mark.parametrize('seq',['HHHHHHHH','HPHHPHHP','PPPPPPPP','HHPHHHHPH'])
@pytest.mark.parametrize('depth',[2,4,6])
@pytest.mark.parametrize('threshold',[0,1,2,3])
def test_turn_union_matches_independent_full_walk_oracle(seq,depth,threshold):
    prefixes=sorted({p[:min(depth,len(seq))] for p in exact_walks(len(seq))})
    for prefix in prefixes[:4]:
        q=Fraction(3,2)
        result=monotone_turn_union(seq,prefix,q,threshold)
        assert (result.event,result.complement,result.distinct_paths)==oracle(seq,prefix,q,threshold)
        base=monotone_family_mass(seq,prefix,q,threshold)
        assert result.event>=base.event
        assert result.complement>=base.complement

@pytest.mark.parametrize('bad',[True,False,3,-1,'2',1.5])
def test_invalid_turn_budget_fails_closed(bad):
    with pytest.raises(ValueError):
        monotone_turn_union('HHHHHH',ROOT,Fraction(2),1,bad)

def test_mixed_sequence_root_nonmonotone_contributions():
    seq='HPHHPHHPHHHH';q=Fraction(3,2)
    result=monotone_turn_union(seq,ROOT,q,2)
    assert result.event==Fraction(459,8)
    assert (result.event,result.complement,result.distinct_paths)==oracle(seq,ROOT,q,2)
