"""Small independent selective-order coverage tests; no biological labels."""
import pytest
import numpy as np
from model_consensus_audit import concordance


def test_full_consensus_correct():
    r=concordance(np.array([1.,2.,3.]),np.array([1.,2.,3.]),np.array([1.,2.,3.]))
    assert r['n_target_distinct_pairs']==3 and r['n_consensus_pairs']==3
    assert r['consensus_coverage']==1 and r['consensus_accuracy']==1
    assert r['consensus_signed_concordance']==1


def test_complete_disagreement_abstains():
    r=concordance(np.array([3.,2.,1.]),np.array([1.,2.,3.]),np.array([1.,2.,3.]))
    assert r['consensus_coverage']==0 and r['consensus_accuracy'] is None
    assert r['model_disagreement']==1


def test_ties_abstain_and_do_not_boost_accuracy():
    r=concordance(np.array([1.,2.,2.]),np.array([1.,2.,3.]),np.array([1.,2.,3.]))
    assert r['n_consensus_pairs']==2 and r['consensus_accuracy']==1
    assert r['gnm_signed_concordance']==2/3
    assert r['consensus_coverage']==2/3