from collections import Counter
from dataclasses import replace
from fractions import Fraction as F
from foldbounds.mutation_switch_10 import (
    ORIGINAL,MUTANT,geometry_fingerprints10,exact_spectrum10,
    population10,certify_ten_mutation,verify_ten_mutation)
from foldbounds.ensemble_uncertainty import lower_probability,upper_probability


def test_larger_ten_residue_cubic_ensemble_exact_compression():
    x=geometry_fingerprints10()
    assert len(x)==1852
    assert sum(count for k,count in x)==308981
    assert sum(count for (edges,b,c),count in x if c)==40932
    for seq in (ORIGINAL,MUTANT):
        assert sum(v for key,v in exact_spectrum10(seq))==308981


def test_ten_residue_exact_uniform_robust_switch():
    c=certify_ten_mutation()
    assert verify_ten_mutation(c)
    assert c.original_robust_upper<F(1,2)<c.mutant_robust_lower
    assert c.original_robust_upper<F(21,100)
    assert F(1,2)<c.mutant_robust_lower<F(51,100)
    assert c.original_nominal_lower==population10(ORIGINAL,F(8),F(3,4))
    assert c.original_nominal_upper==population10(ORIGINAL,F(12),F(5,4))
    assert c.mutant_nominal_lower==population10(MUTANT,F(8),F(3,4))
    assert c.mutant_nominal_upper==population10(MUTANT,F(12),F(5,4))


def test_ten_residue_forgery_and_too_large_errors():
    c=certify_ten_mutation()
    assert not verify_ten_mutation(replace(c,mutant_robust_lower=F(1)))
    assert not verify_ten_mutation(replace(c,fingerprints=1851))
    import pytest
    with pytest.raises(ValueError):
        certify_ten_mutation(factor_lower=F(1,2),factor_upper=F(2))


def test_ten_residue_model_is_mutation_invariant_event():
    assert ORIGINAL=='PHPHPHPHPH' and MUTANT=='HHPHPHPHPH'
    assert all(k==0 for (k,b,c),v in exact_spectrum10(ORIGINAL))
    for seq in (ORIGINAL,MUTANT):
        assert sum(v for (k,b,c),v in exact_spectrum10(seq) if c)==40932