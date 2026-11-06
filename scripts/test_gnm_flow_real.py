"""Optional real-source tests; no NMR target and no external network.

CI runs hermetic tests, while these tests use the precise previous user archive
only when supplied explicitly. Bad mutations must fail closed.
"""
from __future__ import annotations
import os,sys,json,copy
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).parent))
from gnm_flow_independent_replay import verify, integer_flow_check

S1=os.environ.get('FOLDBOUNDS_RNASE_STAGE1')
CERT=os.environ.get('FOLDBOUNDS_FLOW_CERT')
pytestmark=pytest.mark.skipif(not(S1 and CERT),reason='Requires user pinned Stage1 + flow certificate')

def test_full_original_dyadic_witness():
    r=verify(S1,CERT)
    assert r['B_exact']=='72886878439/4294967296'
    assert r['coercivity_improvement']>15.43
    assert r['flow_entries']==119556

def test_wrong_source_rejected(tmp_path):
    path=tmp_path/'bad.zip';path.write_bytes(Path(S1).read_bytes()+b'x')
    with pytest.raises(ValueError,match='Stage1 SHA'):
        verify(path,CERT)

def test_modified_flow_entry_rejected(tmp_path):
    c=json.loads(Path(CERT).read_text());c['integer_flow_numerators'][0][0]+=1
    out=tmp_path/'tamper.json';out.write_text(json.dumps(c))
    with pytest.raises(ValueError,match='conservation'):
        verify(S1,out)

def test_deleted_contact_rejected(tmp_path):
    c=json.loads(Path(CERT).read_text());c['certain_contact_edges'].pop()
    out=tmp_path/'tamper.json';out.write_text(json.dumps(c))
    with pytest.raises(ValueError,match='graph not exact'):
        verify(S1,out)

def test_false_B_rejected(tmp_path):
    c=json.loads(Path(CERT).read_text());c['coercivity_numerator']//=2
    out=tmp_path/'tamper.json';out.write_text(json.dumps(c))
    with pytest.raises(ValueError,match='fabricated B'):
        verify(S1,out)

def test_non_conserved_root_or_parent_rejected(tmp_path):
    c=json.loads(Path(CERT).read_text());c['tree_parent'][0]=0
    out=tmp_path/'tamper.json';out.write_text(json.dumps(c))
    with pytest.raises(ValueError):verify(S1,out)