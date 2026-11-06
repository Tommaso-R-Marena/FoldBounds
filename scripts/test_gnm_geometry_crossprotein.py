"""Hermetic and optional actual-source geometric no-good regression tests."""
import json,os
from pathlib import Path
from fractions import Fraction as F
import pytest
from gnm_geometry_independent_replay import certificate_contradiction,independent_count
from gnm_geometry_crossprotein_independent_replay import replay


def test_generic_sat_count_factoring():
    clauses=[([0,1],[0,0]),([1,2],[0,1]),([6,7,8],[1,1,1])]
    count,components=independent_count(12,clauses)
    assert len(components)==8
    # independent brute global 12-bit count
    expected=0
    for z in range(1<<12):
        if all(any(((z>>i)&1)!=s for i,s in zip(ids,states)) for ids,states in clauses):
            expected+=1
    assert count==expected


def test_strict_farkas_rejects_non_separation():
    # One reference edge whose on/off condition is realizable must not
    # become an accepted contradiction through invented zero weights.
    xyz=[[0,0,0],[9,0,0]]
    with pytest.raises(ValueError,match='multipliers'):
        certificate_contradiction(xyz,1,[(0,1)],[1],[F(0)])


def paths():
    folder=os.environ.get('FOLDBOUNDS_GEOM_RELEASE_DIR')
    if not folder:pytest.skip('Run against original pinned source-only release by setting FOLDBOUNDS_GEOM_RELEASE_DIR')
    return Path(folder)


@pytest.mark.parametrize('protein,expected', [('CspB',87.38288879394531),('Crambin',95.07331848144531)])
def test_real_independent_crossprotein(protein,expected):
    folder=paths()
    r=replay(folder/'FoldBounds_Prospective_SEND_BACK.zip',folder/f'{protein.lower()}_geometric_nogoods.json')
    assert r['status']=='PASS' and abs(r['excluded_percent']-expected)<1e-11


def test_tampered_real_crossprotein_rejected(tmp_path):
    folder=paths()
    p=folder/'cspb_geometric_nogoods.json'
    data=json.loads(p.read_text());data['forbidden_patterns'][0]['lambda']=['0']*2
    out=tmp_path/'bad.json';out.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='Farkas|multipliers'):
        replay(folder/'FoldBounds_Prospective_SEND_BACK.zip',out)