"""Hermetic and optional source-only exhaustive 64-case effective-resistance tests."""
import copy,json,os,random
from fractions import Fraction
from pathlib import Path
import pytest
from gnm_exact_case_source import count_nominal_pairwise_agreement
from gnm_exact_case_independent_replay import certified_count_offline,replay
from gnm_exact_edge_updates import update_contact_signed
from gnm_exact_inverse_producer import grounded,check_scaled_inverse


def test_exact_low_rank_downdate_reverts_three_vertex_triangle():
    triangleN=[[2,1],[1,2]];D=3
    N2,D2=update_contact_signed(triangleN,D,0,2,2,-1,n=3)
    assert (N2,D2)==([[2,1],[1,1]],1)
    A,_=grounded([(0,1),(1,2)],3,2)
    assert check_scaled_inverse(A,N2,D2)
    with pytest.raises(ValueError):
        update_contact_signed([[1]],1,0,1,1,-1,n=2)


def test_distinct_rank_counting_implementations_agree_on_random_intervals():
    rng=random.Random(5736)
    for n in range(2,120,9):
        for _ in range(15):
            z=list(range(n));rng.shuffle(z)
            lo=[Fraction(rng.randint(1,1000),997) for k in range(n)]
            hi=[v+Fraction(rng.randrange(0,30),101) for v in lo]
            a=count_nominal_pairwise_agreement(z,lo,hi)
            b=certified_count_offline(z,lo,hi)
            assert a==b
            if n<=25:
                correct=reverse=0
                for i in range(n):
                    for j in range(i+1,n):
                        x,y=z[i],z[j]
                        correct+=hi[x]<lo[y]
                        reverse+=lo[x]>hi[y]
                assert a==(correct,reverse)


@pytest.fixture
def real():
    names=('FOLDBOUNDS_RNASE_STAGE1','FOLDBOUNDS_EXACT_CERTAIN','FOLDBOUNDS_EXACT_POSSIBLE','FOLDBOUNDS_EXACT_NOMINAL','FOLDBOUNDS_EXACT_64CASES')
    values=[os.environ.get(s) for s in names]
    if not all(values):pytest.skip('Original frozen reference-only exact graph witnesses not available')
    return values


def test_full_real_64case_independent_replay(real):
    result=replay(*real)
    assert result['n_endpoint_scaled_inverses_checked']==128
    assert result['universal_lower']=='1252145/3415191'
    assert result['target_labels_accessed'] is False


def test_real_64case_tamper_reject(real,tmp_path):
    c=json.loads(Path(real[-1]).read_text())
    variants=[
        ('wrong_root',lambda x:x.__setitem__('selected_contacts',x['selected_contacts'][::-1])),
        ('lost_case',lambda x:x.__setitem__('cases',x['cases'][:-1])),
        ('changed_count',lambda x:x['cases'][0].__setitem__('definitely_correct',x['cases'][0]['definitely_correct']+1)),
    ]
    for name,mutate in variants:
        altered=copy.deepcopy(c);mutate(altered)
        path=tmp_path/(name+'.json');path.write_text(json.dumps(altered))
        with pytest.raises(ValueError,match='ABSTAIN'):
            replay(*real[:-1],path)