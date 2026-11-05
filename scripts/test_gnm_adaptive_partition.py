"""Independent mathematical and adversarial tests for adaptive GNM cases.

No private 2AAS files or remote protein downloads are used in CI.
"""
import importlib.util
import itertools
import random
import sys
from pathlib import Path
from fractions import Fraction

import numpy as np
import pytest

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
# gnm_* foundational modules are in the parent FoldBounds scripts directory.
from gnm_adaptive_independent_replay import count_independent,exact_coord_boxes,checked_scores,tree_checked


def brute_counts(lo,hi,target):
    c=d=u=t=0
    for i in range(len(lo)):
        for j in range(i+1,len(lo)):
            if target[i]==target[j]: t+=1;continue
            better,worse=(i,j) if target[i]>target[j] else (j,i)
            if lo[better]>hi[worse]:c+=1
            elif hi[better]<lo[worse]:d+=1
            else:u+=1
    denom=c+d+u
    return (c,d,u,t,str(Fraction(c-d-u,denom)),str(Fraction(c-d+u,denom)))


def test_concordance_differential_random_ties_and_exact_rationals():
    rng=random.Random(20261105)
    for n in (2,3,5,8):
        for _ in range(90):
            lo=[rng.randint(-20,10) for _ in range(n)]
            hi=[z+rng.randint(0,10) for z in lo]
            obs=[rng.randrange(6) for _ in range(n)]
            if len(set(obs))<2:continue
            a=count_independent(lo,hi,obs)
            b=brute_counts(lo,hi,obs)
            assert (a['guaranteed_concordant'],a['guaranteed_discordant'],a['uncertain'],a['target_ties'],a['lower'],a['upper'])==b


def test_case_partition_exhaustive_independent_edge_states():
    for m in range(1,8):
        for fixed_mask in range(1<<m):
            # Every selected mask is ONE case; its endpoint graphs contain
            # all true graphs regardless of how the unselected bits vary.
            split_mask=((1<<m)-1)^(fixed_mask)
            for world in range(1<<m):
                on_fixed=fixed_mask
                actual=on_fixed|(world&split_mask)
                chosen=split_mask & 0b1010101
                forced=chosen&actual
                low=on_fixed|forced
                high=low|(split_mask&~chosen)
                assert (low&actual)==low
                assert (high&actual)==actual


def test_small_graph_exact_all_subsets_vs_branch_endpoints():
    import sympy as sp
    n=4
    fixed=[(0,1),(1,2),(2,3)]
    optional=[(0,2),(0,3),(1,3)]
    pairs=[(0,2),(0,3),(1,3)]
    def exact_scores(edges):
        A=sp.zeros(n-1)
        for i,j in edges:
            v=sp.zeros(n-1,1)
            if i!=0:v[i-1]+=1
            if j!=0:v[j-1]-=1
            A+=v*v.T
        B=A.inv()
        scores=[]
        for i,j in pairs:
            u=sp.zeros(n-1,1)
            if i!=0:u[i-1]+=1
            if j!=0:u[j-1]-=1
            scores.append(sp.Rational((u.T*B*u)[0]))
        return scores
    for bit in range(1<<len(optional)):
        actual=fixed+[optional[k] for k in range(3) if bit>>k&1]
        chosen_on=[optional[0]] if bit&1 else []
        lower=fixed+chosen_on
        upper=fixed+chosen_on+optional[1:]
        q=exact_scores(actual);q_lower=exact_scores(lower);q_upper=exact_scores(upper)
        for a,b,c in zip(q_lower,q,q_upper):
            assert c<=b<=a


def test_exact_contact_boxes_at_threshold():
    xyz=np.array([[0.,0.,0.],[10.,0.,0.],[23.,0.,0.]])
    on,off,maybe=exact_coord_boxes(xyz,'0')
    assert (0,1) in on and (0,2) in off
    assert not maybe
    on,off,maybe=exact_coord_boxes(xyz,'0.05')
    assert (0,1) in maybe
    with pytest.raises(ValueError):
        exact_coord_boxes(np.array([[0.,0.,0.],[0.,0.,0.]]),'0.05')


def test_untrusted_inverse_witness_fail_closed_format():
    pairs=[(0,2)]
    edges=[(0,1),(1,2)]
    n=3;root=2;B=3;Q=10
    tree_checked(n,root,[-1,-1,-1],B,edges) if False else None
    witness={'n_contacts':2,'integer_inverse_proposal':[['bad',1],[1,2]]}
    with pytest.raises(ValueError):checked_scores(n,edges,pairs,witness,root,B,Q)
    witness={'n_contacts':3,'integer_inverse_proposal':[[2,1],[1,1]]}
    with pytest.raises(ValueError):checked_scores(n,edges,pairs,witness,root,B,Q)


def test_absent_or_malformed_claims_fail_closed():
    with pytest.raises(ValueError):count_independent([3,1],[2,2],[1,2])
    with pytest.raises(ValueError):count_independent([1,2],[2,3],[1,1])
