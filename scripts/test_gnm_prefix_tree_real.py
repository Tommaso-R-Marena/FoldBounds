"""Optional genuinely independent attacks against the real sealed RNase files.

Only run with FOLDBOUNDS_RNASE_STAGE1, FOLDBOUNDS_RNASE_STAGE2, and
FOLDBOUNDS_PREFIX_TREE_CERT paths in the environment. Hermetic CI skips these
so it never downloads or substitutes a protected/unknown dataset.
"""
import os,json
from pathlib import Path
import pytest
from gnm_prefix_tree_independent_replay import verify

STAGE1=os.environ.get('FOLDBOUNDS_RNASE_STAGE1')
STAGE2=os.environ.get('FOLDBOUNDS_RNASE_STAGE2')
CERT=os.environ.get('FOLDBOUNDS_PREFIX_TREE_CERT')
pytestmark=pytest.mark.skipif(not all((STAGE1,STAGE2,CERT)),
                 reason='Real sealed experiment files not present in this CI')


def attack(tmp_path,mutator):
    data=json.loads(Path(CERT).read_text())
    mutator(data)
    file=tmp_path/'forged.json'
    file.write_text(json.dumps(data))
    with pytest.raises((ValueError,AssertionError,KeyError)):
        verify(STAGE1,STAGE2,file)


def test_unmodified_source_replays():
    report=verify(STAGE1,STAGE2,CERT)
    assert report['lower']=='497/3415191'
    assert report['endpoint_graphs']==26


def test_missing_leaf_rejected(tmp_path):
    attack(tmp_path,lambda x:x['leaves'].pop())


def test_wrong_split_policy_rejected(tmp_path):
    attack(tmp_path,lambda x:x['source_only_selected_indices'].__setitem__(0,0))


def test_false_target_count_rejected(tmp_path):
    attack(tmp_path,lambda x:x['leaves'][0]['target_concordance_counts'].__setitem__('guaranteed_concordant',0))


def test_wrong_tree_parent_rejected(tmp_path):
    attack(tmp_path,lambda x:x.__setitem__('tree_B',x['tree_B']-1))


def test_wrong_model_graph_reference_rejected(tmp_path):
    attack(tmp_path,lambda x:x['leaves'][0].__setitem__('all_unfixed_off_model',x['leaves'][0]['all_unfixed_on_model']))


def test_invalid_inverse_witness_rejected(tmp_path):
    attack(tmp_path,lambda x:x['models'][0]['integer_inverse_proposal'][0].__setitem__(0,1<<257))


def test_forged_global_bound_rejected(tmp_path):
    attack(tmp_path,lambda x:x.__setitem__('lower_concordance','1/2'))