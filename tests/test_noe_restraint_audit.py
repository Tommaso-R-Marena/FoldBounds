"""Independent exact certificates and optional full pinned NMR dataset checks."""
from dataclasses import replace
from pathlib import Path
import math
import random
import pytest

from foldbounds.noe_restraint_audit import (
    AtomRef, NOERestraint, record_from_row, classify, audit,
    check_matching_witness, load_pinned_reference,
)


DEF = Path(__file__).resolve().parents[1] / 'data' / 'pinned_unambiguous_noe_2png_6nbn.json'


def rec(identity=1, length=4000, lo=2000, hi=5000):
    return NOERestraint(identity, AtomRef(1,'H'), (0,0,0),
                        AtomRef(2,'HA'), (length,0,0),lo,hi)


def test_exact_integer_boundary():
    assert classify(rec(length=5000)) == 'satisfied'
    assert classify(rec(length=5001)) == 'violated'
    assert classify(rec(length=1999)) == 'violated'
    assert classify(rec(length=2000)) == 'satisfied'
    assert classify(rec(length=4980), 10) == 'satisfied'
    assert classify(rec(length=4981), 10) == 'uncertain'
    assert classify(rec(length=5021), 10) == 'violated'
    assert classify(rec(length=1979), 10) == 'violated'
    assert classify(rec(length=1981), 10) == 'uncertain'


def test_strict_type_validation():
    with pytest.raises(ValueError): classify(rec(), True)
    with pytest.raises(ValueError): classify(rec(), -1)
    with pytest.raises(ValueError): classify(replace(rec(), upper_mA=True))
    with pytest.raises(ValueError): record_from_row([1])
    with pytest.raises(ValueError): audit((rec(),rec()))


def test_independent_perturbation_oracle():
    rng=random.Random(49376)
    for _ in range(240):
        p=tuple(rng.randint(-5000,5000) for _ in range(3))
        q=tuple(rng.randint(-5000,5000) for _ in range(3))
        if p==q: continue
        lo=rng.randint(0,5000)
        hi=lo+rng.randint(0,6000)
        rad=rng.randint(0,200)
        r=NOERestraint(1,AtomRef(1,'H'),p,AtomRef(2,'HA'),q,lo,hi)
        outcome=classify(r,rad)
        for _ in range(12):
            da=[rng.uniform(-1,1) for _ in range(3)]
            db=[rng.uniform(-1,1) for _ in range(3)]
            na=math.sqrt(sum(x*x for x in da)) or 1
            nb=math.sqrt(sum(x*x for x in db)) or 1
            aa=[p[k]+rad*da[k]/max(na,1) for k in range(3)]
            bb=[q[k]+rad*db[k]/max(nb,1) for k in range(3)]
            actual=math.dist(aa,bb)
            if outcome=='satisfied': assert lo-1e-7<=actual<=hi+1e-7
            if outcome=='violated': assert actual<lo-1e-7 or actual>hi+1e-7


def test_matching_proof_rejects_forgery():
    r1=rec(1,length=6000)
    r2=replace(rec(2,length=6200),left=AtomRef(3,'H'),right=AtomRef(4,'H'))
    assert audit([r1,r2],10).minimum_atoms_requiring_more_than_epsilon_movement==2
    assert check_matching_witness([r1,r2],10,[1,2])
    assert not check_matching_witness([r1,r2],10,[1,1])
    assert not check_matching_witness([r1,r2],10,[3])
    assert not check_matching_witness([r1,r2],999,[1])
    r3=replace(rec(3,length=6500),left=r1.left)
    assert not check_matching_witness([r1,r2,r3],10,[1,3])


def test_pinned_model_vs_experiment_provenance():
    if not DEF.exists():
        pytest.skip('source-pinned PDBStat data fixture is available on PR branch')
    j=__import__('json').loads(DEF.read_text())
    known={
        '2png': ('b2b8a878bd61b278bfe76013f63dc510fd7834d9',
                 'fb7e0c40aefc358cf3fb6eeb90e83526eaa5ae86',1902,1732,304),
        '6nbn': ('da9769dde30d07a48ad973f1351be4e85e3931cd',
                 '5b13f49fda78a0d98e018902fab353ae73ce9bf4',1957,1694,553),
    }
    parsed=load_pinned_reference(DEF)
    for name,(nef_sha,pdb_sha,rows,groups,n) in known.items():
        dat=j['samples'][name]
        assert (dat['source_nef_sha'],dat['source_pdb_sha'],
                dat['total_noe_rows'],dat['total_noe_groups'],len(parsed[name]))==(
                nef_sha,pdb_sha,rows,groups,n)


@pytest.mark.parametrize('epsilon,expected',[
    (0, {'2png': (299,5,0,5), '6nbn': (548,5,0,4)}),
    (10,{'2png': (299,4,1,4), '6nbn': (548,4,1,4)}),
    (50,{'2png': (294,0,10,0), '6nbn': (548,1,4,1)}),
    (100,{'2png': (289,0,15,0), '6nbn': (539,1,13,1)}),
    (500,{'2png': (117,0,187,0), '6nbn': (219,0,334,0)}),
])
def test_real_data_fixed_error_audit(epsilon,expected):
    if not DEF.exists():
        pytest.skip('source-pinned NMR restraint dataset is available on PR branch')
    ref=load_pinned_reference(DEF)
    for name, target in expected.items():
        result=audit(ref[name],epsilon)
        assert (result.certainly_satisfied,result.certainly_violated,
                result.undecided,result.minimum_atoms_requiring_more_than_epsilon_movement)==target
        assert check_matching_witness(ref[name],epsilon,result.disjoint_witness_ids)


def test_model_violations_remain_disjoint_at_10_milliangstrom():
    if not DEF.exists():
        pytest.skip('source-pinned dataset is available on PR branch')
    ref=load_pinned_reference(DEF)
    assert audit(ref['2png'],10).disjoint_witness_ids==(662,1407,1701,1710)
    assert audit(ref['6nbn'],10).disjoint_witness_ids==(145,470,507,643)
