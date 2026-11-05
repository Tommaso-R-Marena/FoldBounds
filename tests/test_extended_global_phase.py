"""Exhaustive n=16 all-parameter global maximum for 12,870 compositions."""
from fractions import Fraction
from itertools import combinations
from math import comb

import pytest

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.extended_global_phase import certify_extended_global_diagram


@pytest.fixture(scope="module")
def model16():
    atlas = build_atlas(16)
    assert atlas.conformations == 1604149
    assert len(atlas.multiplicities) == 31940
    # Independently enumerate ALL complete anchored walks (not merely a hash).
    assert audit_atlas(atlas, max_length=16)
    landscape = build_landscape(
        atlas, require_audit=False, max_length=16,
        max_edge_subsets=3_000_000
    )
    return atlas, landscape


def test_complete_fixed_eight_h_global_phase_over_all_q(model16):
    atlas, model = model16
    diagram = certify_extended_global_diagram(
        model, 8, width=Fraction(1, 100_000)
    )
    assert diagram.length == 16
    assert diagram.h_count == 8
    assert diagram.population == 12870
    assert len(diagram.coefficientwise_survivors) == 9
    assert diagram.coefficientwise_classes == 5
    assert len(diagram.globally_relevant_classes) == 2

    high_weight = ("HPHHPHPPHHPPHPPH", "HPPHPPHHPPHPHHPH")
    low_weight = ("HPPHPHHPHHPPHPPH", "HPPHPPHHPHHPHPPH")
    assert diagram.globally_relevant_classes == (high_weight, low_weight)
    assert diagram.cell_winners == (low_weight, high_weight)

    crossing = diagram.crossover
    assert crossing is not None
    assert crossing.difference_t == (
        0, -158788, -95134, -23630, -3224, -160, 18, 6
    )
    assert crossing.positive_root_count == 1
    assert crossing.cell_signs == (-1, 1)
    assert Fraction(11) < crossing.roots[0].lo + 1 < Fraction(12)
    assert crossing.roots[0].hi - crossing.roots[0].lo <= Fraction(1, 100_000)

    # Separate exact spectra-to-polynomial reconstruction of the finalists.
    for seq in high_weight + low_weight:
        spectrum = atlas.spectrum(seq)
        poly = tuple(
            sum(count * comb(k, j) for k, count in spectrum.items() if k >= j)
            for j in range(max(spectrum) + 1)
        )
        assert model.polynomial(seq) == poly

    # Exhaustively verify all 12,870 candidates at witness q values. These
    # witness checks complement but do NOT prove continuum root completeness.
    masks = tuple(sum(1 << i for i in s)
                  for s in combinations(range(16), 8))
    for q, winners in [
        (Fraction(2), low_weight),
        (Fraction(10), low_weight),
        (Fraction(20), high_weight),
    ]:
        values = model.all_partitions(q)
        optimum = max(values[mask] for mask in masks)
        maximizing = tuple(
            "".join("H" if mask & (1 << i) else "P" for i in range(16))
            for mask in masks if values[mask] == optimum
        )
        assert maximizing == winners


def test_fail_closed_size_limits():
    model = build_landscape(build_atlas(7))
    with pytest.raises(ValueError):
        certify_extended_global_diagram(model, 3, max_sequences=1)
    with pytest.raises(ValueError):
        certify_extended_global_diagram(model, True)
    with pytest.raises(ValueError):
        certify_extended_global_diagram(model, 3, width=1.0)
    with pytest.raises(ValueError):
        certify_extended_global_diagram(model, 3, max_classes=0)
