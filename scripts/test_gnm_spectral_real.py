"""Optional original-7RSA replay; does not open NMR or network resources."""
from pathlib import Path
import os,sys,json
import pytest
sys.path.insert(0,str(Path(__file__).parent))
from gnm_spectral_sos_independent_replay import verify
S1=os.environ.get('FOLDBOUNDS_RNASE_STAGE1');CERT=os.environ.get('FOLDBOUNDS_SPECTRAL_CERT')
pytestmark=pytest.mark.skipif(not(S1 and CERT),reason='Original exact source and SOS JSON required')

def test_all_source_laplacian_entries():
    r=verify(S1,CERT)
    assert r['B_upper']=='1000000/169821'
    assert r['B_improvement_vs_262']>44.49
    assert r['integer_matrix_entries']==15129

def test_change_source_fails(tmp_path):
    p=tmp_path/'bad.zip';p.write_bytes(Path(S1).read_bytes()+b'tamper')
    with pytest.raises(ValueError,match='source SHA'):
        verify(p,CERT)

def mutate(tmp,change):
    cert=json.loads(Path(CERT).read_bytes());change(cert)
    p=tmp/'mutated.json';p.write_text(json.dumps(cert));return p

def test_wrong_contact_graph_fails(tmp_path):
    p=mutate(tmp_path,lambda c:c['certain_contact_edges'].pop())
    with pytest.raises(ValueError,match='not complete'):verify(S1,p)

def test_false_spectral_scalar_fails(tmp_path):
    p=mutate(tmp_path,lambda c:c.__setitem__('lambda_num',200000))
    with pytest.raises(ValueError,match='diagonally dominant'):verify(S1,p)

def test_wrong_matrix_entry_fails(tmp_path):
    p=mutate(tmp_path,lambda c:c['integer_upper_triangular_W'][0].__setitem__(0,9999999999))
    with pytest.raises(ValueError,match='diagonally dominant'):verify(S1,p)

def test_deleted_witness_entry_fails(tmp_path):
    p=mutate(tmp_path,lambda c:c['integer_upper_triangular_W'][0].pop())
    with pytest.raises(ValueError,match='malformed'):verify(S1,p)

def test_fake_diagonal_margin_fails(tmp_path):
    p=mutate(tmp_path,lambda c:c.__setitem__('min_dd_margin_numerator',-1))
    with pytest.raises(ValueError,match='fake'):verify(S1,p)