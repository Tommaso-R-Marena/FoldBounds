"""Tests for the shared-contact PSD cone: exact algebra and real archive tampering."""
from __future__ import annotations
from fractions import Fraction
from pathlib import Path
import json, os, shutil
import pytest
from gnm_coupled_producer import exact_coupled_sign, force_entries
from gnm_coupled_independent_replay import validate


def inverse_for_graph(n,edges,root):
    # SymPy is a *toy test oracle*, never a trusted production verifier.
    from sympy import Matrix
    vertices=[v for v in range(n) if v!=root]
    matrix=[[0 for _ in vertices] for _ in vertices]
    lookup={v:i for i,v in enumerate(vertices)}
    for a,b in edges:
        for u,v,s in ((a,a,1),(b,b,1),(a,b,-1),(b,a,-1)):
            if u!=root and v!=root:matrix[lookup[u]][lookup[v]]+=s
    inverse=Matrix(matrix).inv()
    return [[Fraction(int(inverse[i,j].p),int(inverse[i,j].q)) for j in range(n-1)]
            for i in range(n-1)]


def quadratic(B,force):
    return sum(si*sj*B[i][j] for i,si in force for j,sj in force)


def scalar_transfer(B,f,g):
    return sum(si*sj*B[i][j] for i,si in f for j,sj in g)


def toy_exact_certifies_all_worlds():
    n=6;root=5
    fixed=[(0,1),(1,2),(2,3),(3,4),(4,5)]
    uncertain=[(0,2),(1,4),(0,4)]
    lo=inverse_for_graph(n,fixed,root)
    hi=inverse_for_graph(n,fixed+uncertain,root)
    all_worlds=[inverse_for_graph(n,fixed+[edge for k,edge in enumerate(uncertain) if mask&(1<<k)],root)
                for mask in range(8)]
    M=[[lo[i][j]-hi[i][j] for j in range(n-1)] for i in range(n-1)]
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    def f(pair):
        return [(x,sign) for x,sign in ((pair[0],1),(pair[1],-1)) if x!=root]
    normal=[quadratic(lo,f(p)) for p in pairs]
    # Iterate only comparisons whose exact source nominal ranking has no tie.
    tested=0
    for pi in range(len(pairs)):
        for pj in range(len(pairs)):
            if normal[pi]>=normal[pj]:continue
            f_i,f_j=f(pairs[pi]),f(pairs[pj])
            l_i,l_j=quadratic(hi,f_i),quadratic(hi,f_j)
            u_i,u_j=quadratic(lo,f_i),quadratic(lo,f_j)
            d_i=u_i-l_i;d_j=u_j-l_j
            cross=scalar_transfer(M,f_i,f_j)
            midpoint=(u_j+l_j)-(u_i+l_i)
            coupled=(midpoint>0 and (midpoint**2+4*cross**2>(d_i+d_j)**2))
            if coupled:
                assert all(quadratic(W,f_j)>quadratic(W,f_i) for W in all_worlds)
                tested+=1
    return tested


def test_coupled_cone_exact_on_all_eight_small_networks():
    assert toy_exact_certifies_all_worlds()>0


def test_cone_strict_condition_rejects_zero_margin():
    # Scalar witness with no strictly separated midpoint must abstain.
    lo=[5,5];hi=[8,8]
    M=[[1,0],[0,1]]
    assert not exact_coupled_sign(lo,hi,M,((0,1),),((1,1),),0,1)


def test_cone_discriminant_random_psd_quadratics():
    # A=R^T R, B=S^T S ensure positive semidefinite,
    # and the desired contraction follows for every integer sample.
    import random
    rng=random.Random(88)
    for _ in range(10000):
        R=[[rng.randrange(-5,6) for _ in range(2)] for _ in range(2)]
        S=[[rng.randrange(-5,6) for _ in range(2)] for _ in range(2)]
        puu=sum(row[0]**2 for row in R);pvv=sum(row[1]**2 for row in R)
        puv=sum(row[0]*row[1] for row in R)
        quu=sum(row[0]**2 for row in S);qvv=sum(row[1]**2 for row in S)
        quv=sum(row[0]*row[1] for row in S)
        assert (puv-quv)**2<=(puu+quu)*(pvv+qvv)


def get_real():
    prefix=os.environ.get('FOLDBOUNDS_COUPLED_REAL_DIR')
    if not prefix:pytest.skip('Requires explicit original frozen reference certificates')
    p=Path(prefix)
    return tuple(p/k for k in (
        'FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip',
        'gnm_certain.json','gnm_possible.json','gnm_nominal.json',
        'gnm_source_exact_64cases.json','gnm_coupled_witnesses.json'))


def test_real_archive_full_exact_integer_replay():
    paths=get_real()
    assert validate(*paths)['witnessed_additional_pairs']==849112


def test_real_archive_rejects_forged_witness(tmp_path):
    paths=get_real()
    doc=json.loads(paths[5].read_text())
    doc['cases'][0]['extra_rank_pairs'][0][1]+=1
    forged=tmp_path/'forged.json';forged.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match='ABSTAIN'):
        validate(*paths[:5],forged)


def test_real_archive_rejects_forged_claim(tmp_path):
    paths=get_real()
    doc=json.loads(paths[5].read_text())
    doc['universal_lower_numerator']+=2
    forged=tmp_path/'forged_claim.json';forged.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match='ABSTAIN'):
        validate(*paths[:5],forged)