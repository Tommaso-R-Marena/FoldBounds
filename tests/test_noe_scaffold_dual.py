"""Exact dual and primal tests for processed 6NBN Val69 NOE geometry."""
from pathlib import Path
from itertools import product
import json,hashlib
import pytest
from foldbounds.noe_scaffold_dual import (parse,verify_no_go,verify_local_geometry_primal,with_delta,check_inherited_restraint)

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/val69_local_noe_dual.json'


def load():
    return parse(json.loads(DATA.read_text()))


def test_real_6nbn_rational_primal_dual_and_delta_threshold():
    d=load()
    assert d['source_pdb_git_blob']=='5b13f49fda78a0d98e018902fab353ae73ce9bf4'
    assert d['source_nef_git_blob']=='da9769dde30d07a48ad973f1351be4e85e3931cd'
    c=verify_no_go(d)
    assert c is not None
    assert c.required_mA2==2637307
    assert c.support_mA2<c.required_mA2
    assert c.delta_mA==103
    assert str(c.support_mA2)=='459873142279194887257392/174477508968339989'
    p=verify_local_geometry_primal(d)
    assert p is not None
    assert p.checked_local_pairs==14
    assert p.worst_distance_change_ticks==211263
    assert p.cap_ticks==211263


@pytest.mark.parametrize('delta,proven',[(0,True),(10,True),(50,True),(79,True),(100,True),(103,True),(104,False),(200,False)])
def test_monotone_exact_infeasibility_gate(delta,proven):
    assert (verify_no_go(with_delta(load(),delta)) is not None)==proven


def test_bad_or_forged_duals_rejected():
    d=load()
    q=json.loads(json.dumps(d))
    q['proof']['H'][0][0]=1
    assert verify_no_go(q) is None
    q=json.loads(json.dumps(d))
    q['proof']['H'][0][3]+=1
    assert verify_no_go(q) is None
    q=json.loads(json.dumps(d))
    q['proof']['HB'][0][4]=0
    assert verify_no_go(q) is None
    q=json.loads(json.dumps(d))
    q['proof']['HB']=q['proof']['HB'][:-1]
    assert verify_no_go(q) is None
    q=json.loads(json.dumps(d))
    q['fixed_heavy_atoms_val69']['CA'][0]+=1
    assert verify_no_go(q) is None
    q=json.loads(json.dumps(d))
    q['processed_restraint_upper_mA']=3600
    assert verify_no_go(q) is None
    q=json.loads(json.dumps(d))
    q['allowed_change_in_each_heavy_to_hydrogen_distance_mA']=103.0
    assert verify_no_go(q) is None
    with pytest.raises(ValueError):with_delta(d,1.5)


def test_bad_or_forged_primal_rejected():
    d=load()
    q=json.loads(json.dumps(d))
    q['constructive_feasible_at_larger_tolerance']['guaranteed_distance_deviation_cap_ticks']=211262
    assert verify_local_geometry_primal(q) is None
    q=json.loads(json.dumps(d))
    q['constructive_feasible_at_larger_tolerance']['moves']['H']=[0,0,0]
    assert verify_local_geometry_primal(q) is None
    q=json.loads(json.dumps(d))
    q['constructive_feasible_at_larger_tolerance']['moves']['HB']=[0,0,0]
    assert verify_local_geometry_primal(q) is None


def test_real_source_restraint_bridge_and_sha():
    d=load()
    fixture=ROOT/'data/pinned_unambiguous_noe_2png_6nbn.json'
    if not fixture.exists():pytest.skip('requires source data from stacked PR #80')
    from foldbounds.noe_restraint_audit import load_pinned_reference
    b=fixture.read_bytes()
    sha=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
    assert sha==d['source_pinned_aggregate_git_blob']
    records=load_pinned_reference(fixture)['6nbn']
    assert check_inherited_restraint(d,records)
    ids=[r.identity for r in records]
    assert ids.count(470)==1
    assert not check_inherited_restraint(d,[r for r in records if r.identity!=470])


def test_direct_perturbation_samples_agree_with_universal_proof():
    d=load()
    h=d['target_hydrogens_val69'];anchors=d['fixed_heavy_atoms_val69']
    def ds(a,b):return sum((x-y)**2 for x,y in zip(a,b))
    n=0
    for dx in product((-40,0,40),repeat=3):
        for dy in product((-40,0,40),repeat=3):
            a=[h['H'][j]+dx[j] for j in range(3)]
            b=[h['HB'][j]+dy[j] for j in range(3)]
            if all(abs(ds(a,v)**.5-ds(h['H'],v)**.5)<=103 for v in anchors.values()) and all(abs(ds(b,v)**.5-ds(h['HB'],v)**.5)<=103 for v in anchors.values()):
                assert ds(a,b)>d['processed_restraint_upper_mA']**2
                n+=1
    assert n>0
