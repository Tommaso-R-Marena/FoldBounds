"""Hermetic differential tests for source-only rational GNM flow certificates."""
from __future__ import annotations
import json,sys
from fractions import Fraction
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).parent))
from gnm_flow_independent_replay import integer_flow_check

EDGES=[(0,1),(0,2),(1,2)]
ROOT=2
FLOW=[[1,1,1],[-1,1,1]] # Q=2, exact two-path conservative flows

def energy(u):
    x=[Fraction(i) for i in u]+[Fraction(0)]
    return sum((x[a]-x[b])**2 for a,b in EDGES)

def test_exact_two_path_flow_cert():
    assert integer_flow_check(EDGES,3,ROOT,2,FLOW)==6
    B=Fraction(6,4)
    for a in range(-7,8):
        for b in range(-7,8):
            u=[a,b]
            assert sum(Fraction(x)**2 for x in u) <= B*energy(u)

def test_single_path_conservation_and_weak_bound():
    W=[[0,1,0],[0,0,1]]
    assert integer_flow_check(EDGES,3,ROOT,1,W)==2
    assert Fraction(6,4)<Fraction(2)

@pytest.mark.parametrize('mutator',[
    lambda W: W.__setitem__(0,[1,1,0]),  # lost flow to ground
    lambda W: W.__setitem__(1,[0,0,0]),  # zero missing vertex
    lambda W: W.__setitem__(1,[-1,1,0]), # conservation broken
    lambda W: W.__setitem__(0,[2,0,1]),  # fabricated divergence
])
def test_reject_corrupted_flows(mutator):
    bad=[row[:] for row in FLOW];mutator(bad)
    with pytest.raises(ValueError,match='conservation'):
        integer_flow_check(EDGES,3,ROOT,2,bad)

def test_reject_missing_noninteger_and_huge():
    for changed in ([[1,1],[1,1,1]],[[1,1,1],[1,1,'0']],[[1,1,1],[1,1,1<<70]]):
        with pytest.raises(ValueError):integer_flow_check(EDGES,3,ROOT,2,changed)

def test_reject_bad_denominators():
    for q in (0,-1,1<<30,True):
        with pytest.raises(ValueError):integer_flow_check(EDGES,3,ROOT,q,FLOW)

def test_reject_swapped_edges():
    with pytest.raises(ValueError,match='conservation'):
        integer_flow_check([(0,2),(0,1),(1,2)],3,ROOT,2,FLOW)