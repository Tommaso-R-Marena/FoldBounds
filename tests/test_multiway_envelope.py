"""Adversarial and exhaustive tests for full all-q multiway HP design.

For n=16 every composition is certified: 65,536 sequences, 25 genuine
global changes, and 42 distinct class winners counted by composition.
No floating-point root solving is used by the certificate.
"""
from fractions import Fraction
from itertools import combinations

import pytest

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.multiway_envelope import certify_multiway_envelope


@pytest.fixture(scope="module")
def atlas_landscape_16():
    atlas = build_atlas(16)
    assert atlas.conformations == 1_604_149
    assert len(atlas.multiplicities) == 31_940
    assert audit_atlas(atlas, max_length=16)
    land = build_landscape(
        atlas, require_audit=False, max_length=16,
        max_edge_subsets=3_000_000,
    )
    return atlas, land


# k, number of actual global switches, number of distinct winning
# polynomial classes as q varies through all real q>1.
EXPECTED = (
    (0, 0, 1), (1, 0, 1), (2, 0, 1), (3, 1, 2),
    (4, 0, 1), (5, 3, 4), (6, 1, 2), (7, 2, 3),
    (8, 1, 2), (9, 4, 5), (10, 3, 4), (11, 3, 4),
    (12, 1, 2), (13, 2, 3), (14, 2, 3), (15, 2, 3),
    (16, 0, 1),
)


def _cell_for_q(cert, q):
    t = q - 1
    if any(event.t_lower <= t <= event.t_upper for event in cert.events):
        raise ValueError("chosen witness q lies in an isolating interval")
    return sum(event.t_upper < t for event in cert.events)


def test_complete_all_65536_sequences_at_every_composition(atlas_landscape_16):
    atlas, land = atlas_landscape_16
    certificates = []
    for k, switches, winner_classes in EXPECTED:
        cert = certify_multiway_envelope(land, k)
        assert cert.population == len(list(combinations(range(16), k)))
        assert len(cert.switches) == switches
        assert len(set(cert.root_free_winners)) == winner_classes
        assert len(cert.root_free_winners) == len(cert.events) + 1
        assert all(event.t_upper - event.t_lower <= Fraction(1, 1_000_000)
                   for event in cert.events)
        assert all(event.t_lower >= 0 for event in cert.events)
        certificates.append(cert)

    assert sum(c.population for c in certificates) == 65_536
    assert sum(len(c.switches) for c in certificates) == 25
    assert sum(len(set(c.root_free_winners)) for c in certificates) == 42

    # Exact rational, whole-landscape witness checks independently confirm
    # maxima in all 17 compositions at three q, but do NOT serve as the
    # proof of the continuum claim; Sturm root completeness does that.
    for q in (Fraction(2), Fraction(10), Fraction(30)):
        values = land.all_partitions(q)
        assert len(values) == 65_536
        for cert in certificates:
            masks = tuple(sum(1 << i for i in indices)
                          for indices in combinations(range(16), cert.h_count))
            optimum = max(values[m] for m in masks)
            expected_names = cert.root_free_winners[_cell_for_q(cert, q)]
            actual_names = tuple(
                "".join("H" if m & (1 << i) else "P" for i in range(16))
                for m in masks if values[m] == optimum
            )
            assert actual_names == expected_names


def test_nine_h_five_winner_phases_and_independent_spectra(atlas_landscape_16):
    atlas, land = atlas_landscape_16
    cert = certify_multiway_envelope(land, 9)
    assert cert.population == 11_440
    assert cert.skyline_sequences == 20
    assert cert.skyline_classes == 10
    assert len(cert.final_classes) == 7
    assert len(cert.events) == 21
    assert len(cert.switches) == 4

    expected_representatives = (
        "HPPHHPHHPHHPHPPH",
        "HPPHPHHHPHHPHPPH",
        "HPPHPHHHHPHPHPPH",
        "HPHPHHHHPHPPHPPH",
        "HHHHPHPHPPHPHPPH",
    )
    winning_reps = tuple(dict.fromkeys(c[0] for c in cert.root_free_winners))
    assert winning_reps == expected_representatives
    switch_ranges = (
        (Fraction(3), Fraction(4)),
        (Fraction(11), Fraction(12)),
        (Fraction(12), Fraction(13)),
        (Fraction(18), Fraction(19)),
    )
    for event, (low, high) in zip(cert.switches, switch_ranges):
        assert low < event.t_lower + 1 < event.t_upper + 1 < high
        assert event.before != event.after

    # Reconstruct finalist Z from the full geometry mask scores, not the
    # polynomial transform or the Sturm checker.
    for seq in expected_representatives:
        for q in (Fraction(2), Fraction(5, 2), Fraction(11), Fraction(20)):
            assert land.partition(seq, q) == atlas.partition(seq, q)
        assert land.contact_spectrum(seq) == atlas.spectrum(seq)


class _SyntheticLandscape:
    """Minimal exact-polynomial fixture, NOT a physical HP atlas."""

    length = 4

    def __init__(self, entries):
        self.entries = entries

    def polynomial(self, sequence):
        return self.entries[sequence]


_BASE = {s: (100,) for s in ("HPPP", "PHPP", "PPHP", "PPPH")}


def test_three_way_global_switches_and_irrelevant_pair_root():
    landscape = _SyntheticLandscape(
        _BASE | {
            "HPPP": (100, 10),
            "PHPP": (100, 7, 2),
            "PPHP": (100, 1, 4),
        }
    )
    cert = certify_multiway_envelope(landscape, 1)
    assert cert.population == 4
    assert len(cert.events) == 3
    assert tuple(e.kind for e in cert.events) == (
        "global_switch", "irrelevant", "global_switch"
    )
    assert len(cert.switches) == 2
    assert all(
        e.t_lower <= root <= e.t_upper
        for e, root in zip(cert.switches, (Fraction(3, 2), Fraction(3)))
    )
    assert cert.root_free_winners == (
        ("HPPP",), ("PHPP",), ("PHPP",), ("PPHP",)
    )


def test_exact_tangency_is_not_mislabeled_as_switch():
    landscape = _SyntheticLandscape(
        _BASE | {
            "HPPP": (100, 2, 0, 1),
            "PHPP": (100, 1, 2),
        }
    )
    cert = certify_multiway_envelope(landscape, 1)
    assert not cert.switches
    assert len(cert.events) == 1
    assert cert.events[0].kind == "global_tangent"
    assert cert.events[0].t_lower == cert.events[0].t_upper == 1
    assert cert.root_free_winners == (("HPPP",), ("HPPP",))


def test_coincident_root_and_invalid_budgets_fail_closed():
    # All three nontrivial class differences vanish at t=1. This checker
    # cannot order shared algebraic roots, and must NOT claim a diagram.
    collide = _SyntheticLandscape(
        _BASE | {
            "HPPP": (100, 3),
            "PHPP": (100, 2, 1),
            "PPHP": (100, 1, 2),
        }
    )
    with pytest.raises(ValueError, match="overlapping or coincident"):
        certify_multiway_envelope(collide, 1, max_refinements=2)

    ordinary = _SyntheticLandscape(
        _BASE | {
            "HPPP": (100, 10), "PHPP": (100, 7, 2),
            "PPHP": (100, 1, 4),
        }
    )
    for opts in (
        {"width": 0.1},
        {"max_sequences": 2},
        {"max_classes": 1},
        {"max_events": 1},
        {"max_refinements": 0},
    ):
        with pytest.raises(ValueError):
            certify_multiway_envelope(ordinary, 1, **opts)
    for k in (-1, 5, True, 1.0):
        with pytest.raises(ValueError):
            certify_multiway_envelope(ordinary, k)
