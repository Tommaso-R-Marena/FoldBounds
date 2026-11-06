"""Exhaustive boolean-world validation of source-only split certificates."""
from fractions import Fraction
from itertools import product
from random import Random

import pytest
pytest.importorskip('sympy')

from foldbounds.gnm_joint_rank import (
    joint_rank_certificate, untrusted_propose_inverse,
)
from foldbounds.gnm_joint_cases import (
    certify_split_cases,rank_source_only_optional_edges,
)
from foldbounds.gnm_integer_updates import ExactEndpointCache
from experiments.real_1d3z_joint_rank import load_ca, DEFAULT_FIXTURE
from experiments.real_1d3z_joint_cases import evaluate
from foldbounds.gnm_joint_rank import classify_contact_boxes


def exact_oracle_diagonal(n, edges):
    from foldbounds.gnm_joint_rank import scaled_grounded_laplacian
    a=scaled_grounded_laplacian(n,edges)
    m=[[Fraction(a[i][j]) for j in range(n)]+
       [Fraction(int(i==j)) for j in range(n)] for i in range(n)]
    for k in range(n):
        i0=next(i for i in range(k,n) if m[i][k])
        m[k],m[i0]=m[i0],m[k]
        factor=m[k][k];m[k]=[x/factor for x in m[k]]
        for j in range(n):
            if j==k:continue
            c=m[j][k]
            m[j]=[x-c*y for x,y in zip(m[j],m[k])]
    return [n*m[i][i+n]-Fraction(1,n) for i in range(n)]


def test_all_small_boolean_worlds_and_untrusted_updates():
    randomizer=Random(775)
    for n in (5,6,7):
        backbone={(i,i+1) for i in range(n-1)}
        rest=[(i,j) for i in range(n) for j in range(i+2,n)]
        for repeat in range(8):
            randomizer.shuffle(rest)
            certain=tuple(sorted(backbone|set(rest[:2])))
            possible=tuple(sorted(backbone|set(rest[:5])))
            cache=ExactEndpointCache(n,certain)
            baseline=certify_split_cases(n,certain,possible,selected_count=0,update_cache=cache)
            improved=certify_split_cases(n,certain,possible,selected_count=2,update_cache=cache)
            assert improved.case_certified_pairs>=baseline.case_certified_pairs
            assert len(improved.conclusions)==4
            for bits in product((False,True),repeat=3):
                state=tuple(sorted(set(certain)|{e for e,s in zip(rest[2:5],bits) if s}))
                diag=exact_oracle_diagonal(n,state)
                for i,j,sign in improved.signed_orders:
                    assert (diag[i]-diag[j])*sign>0


def test_real_1d3z_four_exhaustive_cases():
    model=evaluate(Fraction(1,20),2,fixed_two=True)
    assert model['all_residue_pairs']==2346
    assert model['optional_nonlocal_edges']==15
    assert model['selected_source_only_optional_edges']==[[43,51],[39,42]]
    assert model['exhaustive_cases']==4
    assert model['verified_integer_endpoint_matrices']==8
    assert model['single_case_joint_certified']==1919
    assert model['all_case_joint_certified']==1941
    assert model['newly_certified']==22


def test_real_trace_sensitivity_selection_no_label_access():
    coords=load_ca(DEFAULT_FIXTURE)
    certain,possible=classify_contact_boxes(coords,Fraction(8),Fraction(1,20))
    r=rank_source_only_optional_edges(len(coords),certain,possible,
                                      untrusted_propose_inverse(len(coords),certain))
    assert len(r)==15
    assert [edge for edge,_ in r[:3]]==[(37,40),(44,65),(48,57)]
    assert all(gain>=0 for _,gain in r)
    assert all(r[i][1]>=r[i+1][1] for i in range(len(r)-1))
    with pytest.raises(ValueError):
        certify_split_cases(len(coords),certain,possible,selected_count=3,
                            selected_override=((0,1),(0,3),(0,3)))


def test_real_1d3z_command_line_import_modes_agree():
    """Prevent pytest masking a failing direct-script experiment import."""
    import json
    from pathlib import Path
    import subprocess
    import sys

    repo = Path(__file__).resolve().parents[1]
    script = repo / 'experiments' / 'real_1d3z_joint_cases.py'
    arguments = ('--epsilon', '0.05', '--selected-count', '2', '--fixed-two')
    direct = subprocess.run(
        [sys.executable, str(script), *arguments], cwd=repo,
        capture_output=True, text=True, check=True, timeout=120,
    )
    module = subprocess.run(
        [sys.executable, '-m', 'experiments.real_1d3z_joint_cases', *arguments],
        cwd=repo, capture_output=True, text=True, check=True, timeout=120,
    )
    d = json.loads(direct.stdout)
    assert d == json.loads(module.stdout)
    assert (d['single_case_joint_certified'], d['all_case_joint_certified']) == (1919, 1941)
    for prefix in (
        [sys.executable, str(repo / 'experiments' / 'greedy_1d3z_cases.py')],
        [sys.executable, '-m', 'experiments.greedy_1d3z_cases'],
    ):
        help_result = subprocess.run(
            [*prefix, '--help'], cwd=repo, capture_output=True,
            text=True, check=True, timeout=30,
        )
        assert '--depth' in help_result.stdout


def test_real_1d3z_completed_greedy_three_exact():
    """Complete the third source-only greedy step and verify all eight cases."""
    r = evaluate(Fraction(1, 20), 3, fixed_three=True)
    assert r['all_residue_pairs'] == 2346
    assert r['optional_nonlocal_edges'] == 15
    assert r['selected_source_only_optional_edges'] == [[43, 51], [39, 42], [7, 13]]
    assert r['exhaustive_cases'] == 8
    assert r['verified_integer_endpoint_matrices'] == 16
    assert r['single_case_joint_certified'] == 1919
    assert r['all_case_joint_certified'] == 1947
    assert r['newly_certified'] == 28
    with pytest.raises(ValueError):
        evaluate(Fraction(1, 20), 2, fixed_three=True)
    with pytest.raises(ValueError):
        evaluate(Fraction(1, 20), 3, fixed_two=True, fixed_three=True)
