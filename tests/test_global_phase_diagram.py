"""Certify the complete length-14, six-H partition phase diagram for all q>1."""
from fractions import Fraction
from itertools import combinations

import pytest

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.global_phase_diagram import certify_global_phase_diagram
from foldbounds.fixed_composition_dominance import sequence_from_mask


@pytest.fixture(scope="module")
def model14():
    atlas = build_atlas(14)
    # Full independent walk replay, not just an unaudited histogram.
    assert audit_atlas(atlas, max_length=14)
    assert atlas.conformations == 220375
    assert len(atlas.multiplicities) == 5953
    return build_landscape(atlas, require_audit=False)


def test_exact_two_crossing_global_3003_sequence_envelope(model14):
    diagram = certify_global_phase_diagram(
        model14, 6, width=Fraction(1, 100_000)
    )
    assert diagram.population == 3003
    a, b = diagram.candidates
    assert (a, b) == ("HPPHPHPPHPHPPH", "HPPHPPHHPPHPPH")
    assert diagram.cell_winners == ((b,), (a,), (b,))
    assert model14.polynomial(a) == (
        220375, 140328, 41454, 6812, 620, 36, 2
    )
    assert model14.polynomial(b) == (
        220375, 142046, 35966, 5524, 640, 56, 6
    )
    phase = diagram.comparison
    assert phase is not None
    assert phase.difference_t == (0, -1718, 5488, 1288, -20, -20, -4)
    assert phase.positive_root_count == 2
    assert phase.cell_signs == (-1, 1, -1)
    assert Fraction(1) < phase.roots[0].lo + 1 < Fraction(2)
    assert Fraction(7) < phase.roots[1].lo + 1 < Fraction(8)
    for r in phase.roots:
        assert r.hi - r.lo <= Fraction(1, 100_000)

    # Exhaustively validate the winner among all 3003 sequences at three
    # independent rational witnesses. This does NOT replace the Sturm proof.
    ids = tuple(sum(1 << i for i in c) for c in combinations(range(14), 6))
    for q, winner in [
        (Fraction(11, 10), b),
        (Fraction(2), a),
        (Fraction(10), b),
    ]:
        values = model14.all_partitions(q)
        winner_mask = sum(1 << i for i, char in enumerate(winner) if char == "H")
        optimum = values[winner_mask]
        assert all(values[i] < optimum for i in ids if i != winner_mask)
        assert optimum == model14.partition(winner, q)


def test_zero_and_single_global_phase_cases():
    model = build_landscape(build_atlas(12))
    unique = certify_global_phase_diagram(model, 6)
    assert unique.population == 924
    assert unique.candidates == ("HPPHPHHPHPPH",)
    assert unique.cell_winners == (("HPPHPHHPHPPH",),)
    assert unique.comparison is None

    tied = certify_global_phase_diagram(model, 5)
    assert tied.population == 792
    assert len(tied.candidates) == 2
    assert tied.cell_winners == (tied.candidates,)
    assert tied.comparison is not None
    assert tied.comparison.difference_t == (0,)
    assert tied.comparison.positive_root_count == 0
    with pytest.raises(ValueError, match="more than two"):
        certify_global_phase_diagram(model, 4)
    with pytest.raises(ValueError):
        certify_global_phase_diagram(model, 6, width=0.1)
