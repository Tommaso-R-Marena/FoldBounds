"""Optional full-archive attack tests. Not run in CI (requires private frozen ZIPs).

Set FB_FULL_STAGE1, FB_FULL_STAGE2, FB_FULL_CERT to run, otherwise skipped.
"""
import json,os
from pathlib import Path
import pytest
from gnm_adaptive_independent_replay import replay

@pytest.fixture(scope='module')
def sealed_paths():
    vals=[os.environ.get(k) for k in ('FB_FULL_STAGE1','FB_FULL_STAGE2','FB_FULL_CERT')]
    if any(not v or not Path(v).is_file() for v in vals):
        pytest.skip('full stored original archives required')
    return vals

@pytest.mark.parametrize('attack',[
    'delete_case','fake_concordance','broken_inverse','fake_contact','fake_tree_B',
])
def test_adaptive_certificate_tamper_fails(sealed_paths,attack,tmp_path):
    stage1,stage2,cert_path=sealed_paths
    obj=json.loads(Path(cert_path).read_text())
    if attack=='delete_case':obj['cases'].pop(5)
    elif attack=='fake_concordance':obj['cases'][0]['target_concordance_counts']['lower']='1'
    elif attack=='broken_inverse':obj['models'][0]['integer_inverse_proposal'][5][5]+=10**22
    elif attack=='fake_contact':obj['selected_edges'][0][0]+=1
    elif attack=='fake_tree_B':obj['tree_B']+=100
    test=tmp_path/(attack+'.json');test.write_text(json.dumps(obj))
    with pytest.raises((ValueError,AssertionError)):
        replay(stage1,stage2,test)