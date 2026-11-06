"""Actual deposited protein 1D3Z: independent exact integers, never floats."""
import pytest
from fractions import Fraction

pytest.importorskip('sympy')
from foldbounds.gnm_joint_rank import (
    classify_contact_boxes, joint_rank_certificate, untrusted_propose_inverse,
)
from experiments.real_1d3z_joint_rank import load_ca, DEFAULT_FIXTURE, evaluate


@pytest.mark.parametrize('epsilon,lo_edges,hi_edges,baseline,improved', [
    ('0.025', 154, 161, 2189, 2207),
    ('0.05', 152, 167, 1884, 1919),
    ('0.1', 145, 188, 1065, 1131),
    ('0.2', 137, 207, 629, 699),
    ('0.3', 124, 245, 198, 246),
])
def test_real_protein_joint_certificate(epsilon, lo_edges, hi_edges, baseline, improved):
    result=evaluate(Fraction(epsilon))
    assert result['residues']==69
    assert result['total_pairs']==2346
    assert result['certain_nonlocal_contacts']==lo_edges
    assert result['possible_nonlocal_contacts']==hi_edges
    assert result['separate_diagonal_certified']==baseline
    assert result['joint_covariance_certified']==improved
    assert result['strictly_new_pairs']==improved-baseline
    assert result['exact_witnesses_verified']==2


def test_real_structure_corrupt_witness_rejected():
    xyz=load_ca(DEFAULT_FIXTURE)
    certain,possible=classify_contact_boxes(xyz,Fraction(8),Fraction('0.05'))
    a=untrusted_propose_inverse(len(xyz),certain)
    b=untrusted_propose_inverse(len(xyz),possible)
    from foldbounds.gnm_joint_rank import IntegerInverse
    forged=IntegerInverse(b.numerator, b.denominator+1)
    with pytest.raises(ValueError,match='possible-graph'):
        joint_rank_certificate(len(xyz),certain,possible,a,forged)
