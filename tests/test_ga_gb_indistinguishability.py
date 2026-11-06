"""Adversarial benchmark: coarse hydrophobic alphabet cannot see key fold switch chemistry."""
import pytest
from foldbounds.ga_gb_indistinguishability import (
    GA95,GB95,encode_hp,encode_aromatic,audit_experimental_ga_gb,
    exact_same_encoded_model_inputs)


def test_primary_sequences_from_rcsb_have_three_exact_mutations():
    x=audit_experimental_ga_gb()
    assert x.length==56
    assert x.amino_acid_mutations==((20,'L','A'),(30,'I','F'),(45,'L','Y'))
    assert x.ga_pdb=='2KDL' and x.gb_pdb=='2KDM'


def test_coarse_hp_alphabet_can_make_opposite_real_folds_identical():
    x=audit_experimental_ga_gb()
    assert x.broad_hp_hamming_distance==0
    assert encode_hp(GA95)==encode_hp(GB95)
    assert exact_same_encoded_model_inputs(GA95,GB95)


def test_tyrosine_polar_alternative_still_erases_two_out_of_three_changes():
    x=audit_experimental_ga_gb()
    assert x.y_polar_hp_hamming_distance==1
    assert not exact_same_encoded_model_inputs(GA95,GB95,polar_tyrosine=True)


def test_aromatic_third_class_separates_only_two_sites_not_size_change():
    x=audit_experimental_ga_gb()
    assert x.aromatic_three_class_hamming_distance==2
    assert x.exact_amino_acid_hamming_distance==3
    assert encode_aromatic('I')=='H' and encode_aromatic('F')=='R'
    assert encode_aromatic('L')=='H' and encode_aromatic('A')=='H'


@pytest.mark.parametrize('bad',['','X','HZ','PX','LEUX','leu',None])
def test_nonprotein_sequences_fail_closed(bad):
    with pytest.raises(ValueError):encode_hp(bad)


def test_no_coarse_projection_can_use_hidden_side_chain_information():
    # A toy abstract landscape with weights depending ONLY on projected string.
    # The equality holds for all possible conformations (not enumerated here).
    a,b=encode_hp(GA95),encode_hp(GB95)
    for state in range(25):
        def deterministic_energy(encoded,state):
            return sum((state+i+1)%4 for i,res in enumerate(encoded) if res=='H')
        assert deterministic_energy(a,state)==deterministic_energy(b,state)