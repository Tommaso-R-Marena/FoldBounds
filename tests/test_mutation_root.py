from dataclasses import replace
from fractions import Fraction as F
import pytest
from foldbounds.mutation_root import certify_contact_root,verify_contact_root

@pytest.mark.parametrize('seq,target',[
    ('HHPHPHPH',(F(6),F(7))),
    ('HHPHPHPHPH',(F(7),F(8))),
])
def test_single_positive_robust_root_exact_brackets(seq,target):
    c=certify_contact_root(seq,F(3,4),F(5,6),F(6,5),*target,bits=48)
    assert verify_contact_root(c)
    assert c.positive_root_sign_changes==1
    assert c.bracket_high-c.bracket_low==(target[1]-target[0])/2**48
    if len(seq)==8:
        assert F(6081,1000)<c.bracket_low<c.bracket_high<F(6082,1000)
    else:
        assert F(7906,1000)<c.bracket_low<c.bracket_high<F(7907,1000)

@pytest.mark.parametrize('field,value',[
    ('positive_root_sign_changes',0),
    ('polynomial_coefficients',(F(1),)),
    ('bracket_low',F(10)),
    ('bracket_high',F(1)),
])
def test_forged_critical_roots_rejected(field,value):
    c=certify_contact_root('HHPHPHPH',F(3,4),F(5,6),F(6,5),F(6),F(7))
    assert not verify_contact_root(replace(c,**{field:value}))


def test_no_confident_root_when_wrong_sign_bracket():
    with pytest.raises(ValueError):
        certify_contact_root('HHPHPHPH',F(3,4),F(5,6),F(6,5),F(8),F(10))


def test_parent_with_no_energy_response_has_no_root():
    with pytest.raises(ValueError):
        certify_contact_root('PHPHPHPH',F(3,4),F(5,6),F(6,5),F(6),F(10))