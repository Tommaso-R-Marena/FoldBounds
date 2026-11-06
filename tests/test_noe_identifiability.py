"""Rigidity rank deficits and exact component translation witnesses."""
from pathlib import Path
from dataclasses import replace

import pytest

from foldbounds.noe_restraint_audit import (
    AtomRef, NOERestraint, load_pinned_reference)
from foldbounds.noe_identifiability import (
    evaluate_rigidity_gap, verify_component_translation_witness)

FIXTURE=Path(__file__).resolve().parents[1]/'data/pinned_unambiguous_noe_2png_6nbn.json'


def make(i,a,b):
    return NOERestraint(i,AtomRef(a,'H'),(a,2*a,3*a),
                        AtomRef(b,'H'),(b,2*b,3*b),1000,8000)


def test_exact_components_and_tamper_rejection():
    rows=[make(1,1,2),make(2,2,3),make(3,5,6)]
    cert=evaluate_rigidity_gap(rows)
    assert (cert.atoms,cert.unique_pair_edges,cert.components,
            cert.component_sizes)==(5,3,2,(3,2))
    assert cert.exact_component_translation_parameters==3
    assert cert.minimum_linearized_internal_zero_modes==6
    assert verify_component_translation_witness(rows,cert)
    assert not verify_component_translation_witness(
        rows,replace(cert,components=1))


def test_duplicate_and_inconsistent_inputs_fail_closed():
    with pytest.raises(ValueError):
        evaluate_rigidity_gap([make(1,1,2),make(1,3,4)])
    with pytest.raises(ValueError):
        evaluate_rigidity_gap([make(1,1,2),make(2,1,2)])
    bad=replace(make(2,1,3),left_position_mA=(999,0,0))
    with pytest.raises(ValueError):
        evaluate_rigidity_gap([make(1,1,2),bad])


def test_real_pinned_restraint_graph_no_go():
    if not FIXTURE.exists():
        pytest.skip("pinned fixture supplied by FoldBounds research branch")
    samples=load_pinned_reference(FIXTURE)
    expected={
      '2png':(170,304,9,(148,6,4,2,2,2,2,2,2),143,200,24,7),
      '6nbn':(381,553,11,(160,150,39,12,8,2,2,2,2,2,2),183,584,30,33),
    }
    for name,e in expected.items():
        c=evaluate_rigidity_gap(samples[name])
        assert (c.atoms,c.unique_pair_edges,c.components,c.component_sizes,
                c.independent_cycle_count,
                c.minimum_linearized_internal_zero_modes,
                c.exact_component_translation_parameters,
                c.long_range_restraints)==e
        assert verify_component_translation_witness(samples[name],c)
