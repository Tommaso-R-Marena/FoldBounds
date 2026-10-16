"""Independent exhaustive oracle for finite-memory transfer upper bounds."""
from collections import defaultdict
from dataclasses import replace
from fractions import Fraction
from itertools import product
import pytest
from foldbounds.adaptive import adaptive_refine,verify_refinement_run
from foldbounds.lattice import ROOT,build_frontier
from foldbounds.lookahead import verify_lookahead
from foldbounds.transfer_bound import (
    nonbacktracking_mass_upper,relaxed_cross_spectrum,recent_memory_cross_spectrum)

MOVES=((1,0),(-1,0),(0,1),(0,-1))
def oracle_walks(n):
    paths=[ROOT]
    for _ in range(2,n):
        nxt=[]
        for path in paths:
            x,y=path[-1]
            for dx,dy in MOVES:
                point=(x+dx,y+dy)
                if point not in path:nxt.append(path+(point,))
        paths=nxt
    return paths

def oracle_contacts(seq,path):
    return sum(seq[i]=='H' and seq[j]=='H' and
               abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
               for i in range(len(seq)) for j in range(i+2,len(seq)))

@pytest.mark.parametrize('n',range(2,8))
def test_exhaustive_oracle_mass_upper_for_every_sampled_prefix(n):
    walks=oracle_walks(n)
    seqs=([''.join(x) for x in product('HP',repeat=n)] if n<=5 else
          ['H'*n,'P'*n,('HP'*n)[:n],('PH'*n)[:n],('HHPHHPH'*n)[:n]])
    for seq in seqs:
        by_prefix=defaultdict(list)
        for walk in walks:
            score=oracle_contacts(seq,walk)
            for k in range(2,n): by_prefix[walk[:k]].append(score)
        keys=sorted(by_prefix)
        if len(keys)>12:keys=keys[:6]+keys[len(keys)//2:len(keys)//2+3]+keys[-3:]
        for path in keys:
            scores=by_prefix[path]
            for memory in (2,4,6):
                for q in (Fraction(1),Fraction(3,2),Fraction(2)):
                    for threshold in (-1,0,1,2,4):
                        b=nonbacktracking_mass_upper(seq,path,q,threshold,
                                                     suffix_memory=memory)
                        z=sum((q**k for k in scores),Fraction(0))
                        ev=sum((q**k for k in scores if k>=threshold),Fraction(0))
                        assert z<=b.partition
                        assert ev<=b.event
                        assert z-ev<=b.complement

@pytest.mark.parametrize('seq',['HHHHHH','HPHPHHHH','HHHHHHHHHH','PPHHPHPHPH'])
@pytest.mark.parametrize('memory',[2,3,4,5,6,8])
def test_nested_transfer_bounds(seq,memory):
    b=nonbacktracking_mass_upper(seq,ROOT,Fraction(3,2),2,
                                 suffix_memory=memory)
    old=nonbacktracking_mass_upper(seq,ROOT,Fraction(3,2),2)
    assert b.partition<=old.partition
    assert b.event<=old.event
    assert b.complement<=old.complement
    assert b.relaxed_walk_count<=old.relaxed_walk_count
    if memory>=3:
        prev=nonbacktracking_mass_upper(seq,ROOT,Fraction(3,2),2,
                                        suffix_memory=memory-1)
        assert b.partition<=prev.partition
        assert b.event<=prev.event
        assert b.complement<=prev.complement

def test_original_memory_two_is_exactly_preserved():
    for n in range(2,10):
        seq=('HHPHHHPPH'*n)[:n]
        paths=sorted({w[:k] for w in oracle_walks(n) for k in range(2,n+1)})
        for path in paths[:25]:
            assert recent_memory_cross_spectrum(seq,path,2)==relaxed_cross_spectrum(seq,path)

def test_square_loops_excluded_by_four_memory():
    path=ROOT
    before=nonbacktracking_mass_upper('P'*10,path,Fraction(1),1,suffix_memory=2)
    after=nonbacktracking_mass_upper('P'*10,path,Fraction(1),1,suffix_memory=4)
    assert after.relaxed_walk_count<before.relaxed_walk_count
    assert after.partition<before.partition
    assert after.event==before.event==0

@pytest.mark.parametrize('seq',['H'*8,'H'*12,'HHPHPHHH','P'*12])
@pytest.mark.parametrize('q',[Fraction(1),Fraction(3,2),Fraction(2)])
@pytest.mark.parametrize('threshold',[0,1,2])
def test_integrated_bound_tightens_monotonically(seq,q,threshold):
    cert=build_frontier(seq,q,threshold,2)
    previous=None
    for memory in (2,4,6):
        bounds=verify_lookahead(cert,lookahead=0,use_matching=True,
                        use_transfer=True,use_monotone=True,
                        transfer_memory=memory)
        if previous is not None:
            assert previous.partition_lower==bounds.partition_lower
            assert bounds.partition_upper<=previous.partition_upper
            assert bounds.event_upper<=previous.event_upper
            assert bounds.complement_upper<=previous.complement_upper
            assert previous.probability_lower<=bounds.probability_lower
            assert bounds.probability_upper<=previous.probability_upper
        previous=bounds

@pytest.mark.parametrize('bad',[True,0,1,9,2.5,'4',None])
def test_invalid_memory_rejected(bad):
    with pytest.raises(ValueError):
        nonbacktracking_mass_upper('HHHH',ROOT,Fraction(2),1,suffix_memory=bad)
    with pytest.raises(ValueError):
        verify_lookahead(build_frontier('HHHH',Fraction(2),1,2),
                         use_transfer=True,transfer_memory=bad)

def test_adaptive_replay_records_memory():
    cert=build_frontier('H'*10,Fraction(2),1,2)
    run=adaptive_refine(cert,steps=3,use_matching=True,use_transfer=True,
                        transfer_memory=4,use_monotone=True,policy='mass')
    assert run.transfer_memory==4
    assert run.bounds==verify_refinement_run(run)
    with pytest.raises(ValueError):
        verify_refinement_run(replace(run,transfer_memory=2))
