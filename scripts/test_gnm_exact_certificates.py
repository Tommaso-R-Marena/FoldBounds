"""Synthetic exact and adverse corpus for GNM; frozen RNase ZIP optional."""
import io
import json
import random
from fractions import Fraction
from pathlib import Path
import zipfile
import tempfile
import numpy as np
import pytest
import sympy as sp
import gnm_exact_certificates as g


def make_case(n, extra=0):
    # Unit-spaced 1-D line, 3/2-A cutoff makes a connected chain.
    coords=np.column_stack([np.arange(n,dtype=float),np.zeros(n),np.zeros(n)])
    edges,neigh=g.exact_contact_graph(coords,Fraction(9,4))
    assert len(edges)==n-1
    anchor,B,_,_=g.choose_tree(neigh)
    A,idx=g.grounded_laplacian(edges,n,anchor)
    Q=1<<24
    witness=g.untrusted_inverse_proposal(A,Q)
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    result=g.verify_and_score(coords,np.array([i for i,j in pairs]),
           np.array([j for i,j in pairs]),Q,witness,Fraction(9,4))
    true=sp.Matrix(A).inv()
    for pair,(lo,hi) in zip(pairs,result['score_numerators']):
        z=sp.zeros(n-1,1)
        for i,sgn in zip(pair,(1,-1)):
            if i!=anchor:z[idx[i]]+=sgn
        exact=(z.T*true*z)[0]
        assert Fraction(lo,Q*Q)<=Fraction(int(sp.numer(exact)),int(sp.denom(exact)))<=Fraction(hi,Q*Q)
    assert result['strict_adjacent_count'] <= len(pairs)-1


@pytest.mark.parametrize('n',[3,4,5,6,7,8])
def test_exact_gnm_matches_symbolic_solution(n):
    make_case(n)


def test_small_q_abstains_when_residual_large():
    coords=np.array([[0.,0,0],[1,0,0],[2,0,0],[3,0,0]],float)
    edges,nb=g.exact_contact_graph(coords,Fraction(9,4))
    root,B,_,_=g.choose_tree(nb)
    A,idx=g.grounded_laplacian(edges,4,root)
    N=[[0]*3 for _ in range(3)]
    res=g.verify_and_score(coords,np.array([0,0,0]),np.array([1,2,3]),1,N,Fraction(9,4))
    assert res['strict_adjacent_count']==0
    assert all(lo<=hi for lo,hi in res['score_numerators'])


def test_tamper_and_invalid_graph():
    x=np.array([[0.,0,0],[1,0,0],[2,0,0],[3,0,0]],float)
    with pytest.raises(ValueError,match='coincident'):
        g.exact_contact_graph(np.array([[0.,0,0],[0,0,0],[1,0,0]]))
    with pytest.raises(ValueError,match='disconnected'):
        g.choose_tree([set([1]),set([0]),set()])
    with pytest.raises(ValueError):
        g._validate_inputs(x,np.array([1,1]),np.array([2,2]))
    with pytest.raises(ValueError):
        g._validate_inputs(x,np.array([2]),np.array([1]))
    with pytest.raises(ValueError):
        g.verify_integer_witness([[1,0],[0,1]],[[1,2],[3,True]],1)
    with pytest.raises(ValueError):
        g.verify_integer_witness([[1,0],[0,1]],[[1,2],[3,4]],0)


def test_cutoff_exact_boundary():
    x=np.array([[0.,0,0],[3,4,0],[12,0,0]],float)
    edges,_=g.exact_contact_graph(x,Fraction(25))
    assert (0,1) in edges
    assert (0,2) not in edges


def test_real_zip_integrity_if_available():
    p=Path('/mnt/data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip')
    if not p.exists():pytest.skip('Stage1 dataset not part of GitHub CI')
    data,_=g.open_sealed_stage1(p)
    assert len(data['reference_coords'])==124
    with tempfile.TemporaryDirectory() as t:
        dest=Path(t)/'certificate.json'
        c=g.make_certificate(p,dest)
        assert c['strict_adjacent_count']==2613
        g.verify_certificate(p,dest)
        forged=json.loads(dest.read_text())
        forged['score_numerators'][0][0] += 1
        dest.write_text(json.dumps(forged))
        with pytest.raises(ValueError,match='mismatch'):
            g.verify_certificate(p,dest)
        dest.write_text(json.dumps(c))
        forged=json.loads(dest.read_text())
        forged['integer_inverse_proposal'][0][0] += 1
        dest.write_text(json.dumps(forged))
        with pytest.raises(ValueError,match='mismatch'):
            g.verify_certificate(p,dest)


if __name__=='__main__':
    import sys
    raise SystemExit(pytest.main([__file__,'-q']))