"""Exhaustive optional-edge state and adversarial coordinate-box tests."""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import random
import numpy as np
import pytest
import sympy as sp
import gnm_topology_boxes as g
import gnm_exact_certificates as ec


def actual_resistance(edges,n,pair):
    neighbors=[set() for _ in range(n)]
    for i,j in edges:neighbors[i].add(j);neighbors[j].add(i)
    root,B,_,_=ec.choose_tree(neighbors)
    L,idx=ec.grounded_laplacian(edges,n,root)
    z=sp.zeros(n-1,1)
    if pair[0]!=root:z[idx[pair[0]]]+=1
    if pair[1]!=root:z[idx[pair[1]]]-=1
    return (z.T*sp.Matrix(L).inv()*z)[0]


def test_extreme_graph_encloses_all_optional_subsets():
    # A connected 5-residue spatial chain with optional contacts near 10A.
    xyz=np.asarray([[3.*i,0,0] for i in range(5)],dtype=float)
    i=np.array([a for a in range(5) for b in range(a+1,5)])
    j=np.array([b for a in range(5) for b in range(a+1,5)])
    epsilon=Fraction(1)
    fixed,off,optional=g.classify_coord_box(xyz,epsilon)
    assert optional and len(optional)<=10
    result=g.robust_scores(xyz,i,j,epsilon,Q=1<<36)
    D=result['denominator']
    for state in range(1<<len(optional)):
        edges=fixed+[e for bit,e in enumerate(optional) if (state>>bit)&1]
        for idx,pair in enumerate(zip(i,j)):
            real=actual_resistance(edges,5,tuple(map(int,pair)))
            assert Fraction(result['lower_score_numerators'][idx],D)<=real
            assert real<=Fraction(result['upper_score_numerators'][idx],D)


def test_uncertainty_contact_classification_exact_boundary():
    x=np.asarray([[0.,0.,0.],[10.,0,0.],[30.,0,0.]])
    fixed,off,uncertain=g.classify_coord_box(x,Fraction(0))
    assert (0,1) in fixed
    assert (0,2) in off
    assert not uncertain
    fixed,off,uncertain=g.classify_coord_box(x,Fraction('0.001'))
    assert (0,1) in uncertain
    with pytest.raises(ValueError,match='coincident'):
        g.classify_coord_box(np.asarray([[0.,0,0],[.2,0,0],[6.,0,0]]),Fraction('0.5'))


def test_persisted_results_when_stage1_available():
    path=Path('/mnt/data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip')
    if not path.exists():pytest.skip('Frozen input not shipped in GitHub CI')
    a,_=ec.open_sealed_stage1(path)
    for eps,possible,all_rank,top10 in [
        ('0.0005',0,3415191,262),('0.005',5,3361308,258),('0.05',108,2270083,98)]:
        result=g.robust_scores(a['reference_coords'],a['pair_i'],a['pair_j'],Fraction(eps))
        assert result['optional_edges']==possible
        assert result['certified_pair_order']==all_rank
        assert result['guaranteed_top_k']==top10


if __name__=='__main__':
    raise SystemExit(pytest.main([__file__,'-q']))