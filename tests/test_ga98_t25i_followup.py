"""Fail-closed genotype identity tests for additional published GA/GB fold switches."""
import pytest
from experiments.ga98_gb98_transfer import GA98,GB98
from experiments.ga98_t25i_followup import (
    GB98_T25I,GB98_T25I_L20A,analyze_followup
)


def test_reciprocal_mutant_sequences_match_published_cascade():
    assert len(GA98)==len(GB98)==len(GB98_T25I)==len(GB98_T25I_L20A)==56
    assert [(i+1,x,y) for i,(x,y) in enumerate(zip(GA98,GB98)) if x!=y]==[(45,"L","Y")]
    assert [(i+1,x,y) for i,(x,y) in enumerate(zip(GB98,GB98_T25I)) if x!=y]==[(25,"T","I")]
    assert [(i+1,x,y) for i,(x,y) in enumerate(zip(GB98_T25I,GB98_T25I_L20A)) if x!=y]==[(20,"L","A")]
    assert [(i+1,x,y) for i,(x,y) in enumerate(zip(GA98,GB98_T25I_L20A)) if x!=y]==[(20,"L","A"),(25,"T","I"),(45,"L","Y")]


def test_reject_untrusted_external_inputs_before_scoring(tmp_path):
    file=tmp_path/"not_reference.json"
    file.write_text('{"codes":{}}')
    with pytest.raises((ValueError,KeyError)):
        analyze_followup(b"wrong",b"wrong",b"wrong",b"wrong",file)
