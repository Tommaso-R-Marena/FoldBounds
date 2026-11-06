"""Exact grounded-GNM certificates and no-target graph update regression tests."""
from __future__ import annotations
import os,json,copy
from pathlib import Path
import pytest
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix
from gnm_exact_inverse_producer import check_scaled_inverse,grounded,pair_num,create
from gnm_exact_inverse_independent_replay import validate_matrix,strict_bigint,check
from gnm_exact_edge_updates import add_edge_exact,update_certain_to_nominal


def test_exact_integer_inverse_simple_chain():
    A,idx=grounded([(0,1),(1,2)],3,2)
    N,D=[[2,1],[1,1]],1
    assert A==[[1,-1],[-1,2]]
    assert check_scaled_inverse(A,N,D)
    assert pair_num(0,2,N,idx)==2
    assert pair_num(0,1,N,idx)==1
    assert pair_num(1,2,N,idx)==1


def test_exact_edge_update_eliminates_reinversion():
    N,D=([[2,1],[1,1]],1)
    N2,D2,info=add_edge_exact(N,D,0,2,2,n=3)
    assert (N2,D2)==([[2,1],[1,2]],3)
    edges=[(0,1),(1,2),(0,2)]
    A,_=grounded(edges,3,2)
    assert check_scaled_inverse(A,N2,D2)
    assert pair_num(0,2,N2,{0:0,1:1})==2
    assert info['next_D_digits']==1


def test_scaled_inverse_malicious_changes_rejected():
    A=[[2,-1],[-1,2]]
    N=[[2,1],[1,2]]
    assert check_scaled_inverse(A,N,3)
    for Nbad,Dbad in [([[3,1],[1,2]],3),(N,2),(N,-3),([[2,1],[0,2]],3)]:
        with pytest.raises(ValueError,match='ABSTAIN'):
            check_scaled_inverse(A,Nbad,Dbad)


def test_independent_bigint_validation_and_bad_graph():
    for s in ['1.2','01','','-','3 4']:
        with pytest.raises(ValueError,match='ABSTAIN'):
            strict_bigint(s)
    assert strict_bigint('-123')==-123
    assert strict_bigint('0')==0
    with pytest.raises(ValueError,match='ABSTAIN'):
        validate_matrix([(0,1)],[[1]],1,1)


def test_proposal_sympy_is_not_in_checker():
    # The verification code deliberately never imports SymPy.
    src=Path(__file__).with_name('gnm_exact_inverse_independent_replay.py').read_text()
    assert 'import sympy' not in src and 'from sympy' not in src
    assert 'gnm_exact_inverse_producer' not in src
    assert 'gnm_exact_edge_updates' not in src


@pytest.fixture
def real_data():
    stage=os.getenv('FOLDBOUNDS_RNASE_STAGE1')
    nominal=os.getenv('FOLDBOUNDS_EXACT_NOMINAL')
    certain=os.getenv('FOLDBOUNDS_EXACT_CERTAIN')
    if not (stage and nominal and certain):
        pytest.skip('Original pinned Stage1 and exact source-only certificates required')
    return Path(stage),Path(nominal),Path(certain)


def test_real_replay_and_update_no_target(real_data):
    src,nominal,certain=real_data
    n=check(src,nominal);c=check(src,certain)
    assert (n['n_contacts'],n['n_pairs'],n['exact_ties'])==(1020,2614,0)
    assert (c['n_contacts'],c['n_pairs'],c['exact_ties'])==(972,2614,23)
    x=update_certain_to_nominal(src,certain,nominal)
    assert x['n_added_contacts']==48 and x['all_integer_identities']
    assert x['independently_recomputed_exact_inverse_equal']


def test_real_tamper_rejected(real_data,tmp_path):
    src,nominal,_=real_data
    c=json.loads(nominal.read_text())
    variants=[
        ('D',str(int(c['D'])+1)),
        ('N','entry'),
        ('edges',c['edges'][:-1]),
        ('pair_numerators',c['pair_numerators'][:-1]),
        ('canonical_rank_indices',c['canonical_rank_indices'][::-1]),
        ('n_exact_ties',1),
    ]
    for key,value in variants:
        t=copy.deepcopy(c)
        if key=='N':t['N'][0][0]=str(int(t['N'][0][0])+1)
        else:t[key]=value
        dest=tmp_path/(key+'.json');dest.write_text(json.dumps(t))
        with pytest.raises(ValueError,match='ABSTAIN'):
            check(src,dest)


def test_independent_exact_recompute_nominal(real_data,tmp_path):
    src,nominal,_=real_data
    out=tmp_path/'new.json'
    report=create(src,out,'nominal')
    assert out.read_bytes()==nominal.read_bytes(), 'Exact producer not byte reproducible'
    assert report['all_integer_identities_checked']