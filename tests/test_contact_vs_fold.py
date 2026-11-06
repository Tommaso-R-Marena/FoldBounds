"""Independent same-model dual audit: compactness is not fold selectivity."""
from dataclasses import replace
from fractions import Fraction

import pytest

from foldbounds.contact_vs_fold import certify_contact_vs_fold
from foldbounds.cubic_groundstate_no_go import certify_cubic_groundstate_no_go
from foldbounds.local_rewire_confidence import (
    StructuralModel, hp_pair_weights, make_macro_target
)


def model_from_shape(dims, seq, q, lower=Fraction(1, 2), upper=Fraction(2)):
    path = make_macro_target(*dims)
    assert len(seq) == len(path)
    return StructuralModel(seq, path, hp_pair_weights(Fraction(q)), lower, upper)


def independently_count_hh(seq, path):
    """O(n^2) oracle with no checker scoring helpers."""
    return sum(
        seq[i] == seq[j] == "H"
        and sum(abs(a - b) for a, b in zip(path[i], path[j])) == 1
        for i in range(len(path)) for j in range(i + 2, len(path))
    )


def test_all_h_216_coarse_contact_event_near_one_but_exact_fold_tiny():
    model = model_from_shape((3, 3, 3), "H"*216, 4)
    cert = certify_contact_vs_fold(model, (3, 3, 3), 108)
    assert cert.residues == cert.hydrophobic_residues == 216
    assert cert.polar_residues == 0
    assert cert.target_hh_contacts == independently_count_hh(
        model.sequence, model.target
    )
    assert cert.target_hh_contacts == 325
    assert cert.contact_event_probability_lower > Fraction(99, 100)
    assert cert.target_microstate_probability_upper < Fraction(1, 10**8)
    assert cert.target_rigid_fold_probability_upper < Fraction(1, 10**6)
    assert cert.proves_event_majority_without_unique_fold_majority
    ground = certify_cubic_groundstate_no_go(
        6, q=Fraction(4), correction_lower=Fraction(1,2),
        correction_upper=Fraction(2))
    assert ground.globally_optimal_for_attractive_hh
    assert cert.target_hh_contacts == ground.target_contacts
    assert cert.target_microstate_probability_upper == (
        ground.confidence.microstate_probability_upper
    )


def test_216_sequence_specific_hp_model_consistent_across_certificates():
    seq = "".join("P" if i % 11 == 0 else "H" for i in range(216))
    model = model_from_shape((3, 3, 3), seq, 5)
    cert = certify_contact_vs_fold(model, (3, 3, 3), 54)
    assert cert.hydrophobic_residues == 196
    assert cert.polar_residues == 20
    assert cert.target_hh_contacts == independently_count_hh(
        seq, model.target
    )
    assert cert.contact_event_probability_lower > Fraction(99, 100)
    assert cert.target_microstate_probability_upper < Fraction(1, 10**8)
    assert cert.target_rigid_fold_probability_upper < Fraction(1, 10**6)
    assert cert.structure_certificate.exactly_isoenergetic_conformations > 1000
    assert cert.proves_event_majority_without_unique_fold_majority


def test_512_all_h_certifies_both_same_model_claims():
    model = model_from_shape((4, 4, 4), "H"*512, 4)
    cert = certify_contact_vs_fold(model, (4, 4, 4), 256)
    assert cert.target_hh_contacts == 2 * 512 - 3 * 8**2 + 1
    assert cert.contact_event_probability_lower > Fraction(99, 100)
    assert cert.target_microstate_probability_upper < Fraction(1, 10**18)


@pytest.mark.parametrize("override", [
    {"macro_dimensions": (3, 3, 2), "event": 54},
    {"macro_dimensions": (3, 3, 3.0), "event": 54},
    {"macro_dimensions": (3, 3, 3), "event": 0},
    {"macro_dimensions": (3, 3, 3), "event": True},
    {"macro_dimensions": (3, 3, 3), "event": 100000},
])
def test_model_mismatch_abstains(override):
    model = model_from_shape((3, 3, 3), "H" * 216, 4)
    with pytest.raises(ValueError):
        certify_contact_vs_fold(
            model, override["macro_dimensions"], override["event"]
        )


def test_forged_pair_model_rejected_not_merged_with_hp_certificate():
    model = model_from_shape((3, 3, 3), "H"*216, 4)
    for weights in (
        {("H", "H"): Fraction(4)},
        {("H", "H"): Fraction(4), ("H", "P"): Fraction(3),
         ("P", "P"): Fraction(1)},
        {("H", "H"): Fraction(4), ("H", "P"): Fraction(1),
         ("P", "P"): Fraction(2)},
        {("H", "H"): Fraction(1,2), ("H", "P"): Fraction(1),
         ("P", "P"): Fraction(1)},
        {("H", "H"): 4, ("H", "P"): Fraction(1),
         ("P", "P"): Fraction(1)},
    ):
        with pytest.raises(ValueError):
            certify_contact_vs_fold(
                replace(model, pair_weights=weights), (3, 3, 3), 108
            )
    with pytest.raises(ValueError):
        certify_contact_vs_fold(
            replace(model, sequence="A"*216), (3, 3, 3), 108
        )


def test_target_outside_event_abstains():
    # Non-HH contact event cannot be assigned to a target lacking HH,
    # even if it has a high-scoring local motif family elsewhere.
    model = model_from_shape((3, 3, 3), "P"*216, 4)
    with pytest.raises(ValueError, match="outside"):
        certify_contact_vs_fold(model, (3, 3, 3), 1)
