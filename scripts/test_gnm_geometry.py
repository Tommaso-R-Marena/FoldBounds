from __future__ import annotations
import hashlib,json,random,sys
from fractions import Fraction as F
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from gnm_geometry_independent_replay import (
    replay,certificate_contradiction,independent_count,checked_source,
    optional_edge_indices)

BASE=Path(__file__).resolve().parent
SOURCE=BASE/'work/data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip'
CERT=BASE/'geometry_nogoods.json'

def test_exact_minimax_farkas_toy():
    # Impossible constraints x <=-1/2 and -x<=-1/2 within [-.05,.05]
    assert F(1)>F(1,10)
    # Component enumeration excludes 1 of 4 independent assignments.
    count,components=independent_count(2,[([0,1],[0,0])])
    assert count==3 and len(components)==1
    assert components[0]['valid_assignments']==3

def test_component_fraction_recomputed():
    patterns=[([0,1],[0,0]),([1,2],[0,0]),([2,3,4],[1,0,1])]
    count,parts=independent_count(108,patterns)
    # Independently check via 5-bit brute force (103 other unconstrained bits)
    good=0
    for assignment in range(1<<5):
        bits=[(assignment>>i)&1 for i in range(5)]
        good+=all(any(bits[i]!=v for i,v in zip(ids,states)) for ids,states in patterns)
    assert count==good*2**103

def test_bad_metadata_rejected(tmp_path):
    c=json.loads(CERT.read_text());c['box_epsilon_A']='1/19'
    p=tmp_path/'bad.json';p.write_text(json.dumps(c))
    with pytest.raises(ValueError,match='metadata'):
        replay(SOURCE,p)

def test_dropped_nogood_rejected(tmp_path):
    c=json.loads(CERT.read_text());c['forbidden_patterns'].pop()
    p=tmp_path/'bad.json';p.write_text(json.dumps(c))
    with pytest.raises(ValueError,match='incomplete'):
        replay(SOURCE,p)

def test_forged_weight_rejected(tmp_path):
    c=json.loads(CERT.read_text());c['forbidden_patterns'][0]['lambda']=['0','0']
    p=tmp_path/'bad.json';p.write_text(json.dumps(c))
    with pytest.raises(ValueError,match='Farkas'):
        replay(SOURCE,p)

def test_forged_excluded_count_rejected(tmp_path):
    c=json.loads(CERT.read_text());c['excluded_fraction']='1'
    p=tmp_path/'bad.json';p.write_text(json.dumps(c))
    with pytest.raises(ValueError,match='forged'):
        replay(SOURCE,p)

def test_real_geometry_certificate():
    result=replay(SOURCE,CERT)
    assert result['status']=='PASS'
    assert result['at_least_excluded_fraction']=='32010617/33554432'
    assert result['n_nogood_certificates']==20

def test_exact_random_grid_structural_sanity():
    # A non-proof sanity check: physically generated coordinate perturbations
    # cannot realize ANY of the 20 exact no-good binary contact patterns.
    coords,den,_=checked_source(SOURCE)
    optional=optional_edge_indices(coords,den)
    c=json.loads(CERT.read_text());rnd=random.Random(20261105)
    delta=den//20
    for _ in range(150):
        perturb=[[rnd.choice((-delta,0,delta)) for j in range(3)] for i in range(124)]
        states=[]
        for i,j in optional:
            d2=sum((coords[i][k]+perturb[i][k]-coords[j][k]-perturb[j][k])**2 for k in range(3))
            states.append(int(d2<=100*den*den))
        assert not any(all(states[i]==value for i,value in zip(pattern['edge_ids'],pattern['forbidden_on_bits']))
                       for pattern in c['forbidden_patterns'])