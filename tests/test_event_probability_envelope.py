"""Exact all-q global optimization for 65 model event-design problems.

All 4,096 length-12 H/P sequences, thirteen H-count compositions and
threshold events K>=1,...,5. Proof obligations are checked with exact
integer/rational arithmetic and independent complete-walk auditing.
"""
from fractions import Fraction
from itertools import combinations

import pytest

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.universal_landscape import audit_landscape, build_landscape
from foldbounds.event_probability_envelope import certify_event_probability_envelope


@pytest.fixture(scope="module")
def oracle12():
    atlas = build_atlas(12)
    assert atlas.conformations == 30_073
    assert audit_atlas(atlas)
    model = build_landscape(atlas)
    # Independent full-sequence polynomial oracle, not only a selected
    # mutation case: checks every one of the 4096 HP sequences.
    assert audit_landscape(model, atlas)
    direct = {}
    for mask in range(1 << 12):
        sequence = "".join("H" if mask & (1 << i) else "P" for i in range(12))
        spectrum = atlas.spectrum(sequence)
        assert model.contact_spectrum(sequence) == spectrum
        direct[sequence] = spectrum
    return model, direct


def _probability(spectrum, threshold, q):
    z = sum((n * q ** k for k, n in spectrum.items()), Fraction(0))
    e = sum((n * q ** k for k, n in spectrum.items()
             if k >= threshold), Fraction(0))
    return e / z


def _cell_at(cert, q):
    t = q - 1
    if any(e.t_lower <= t <= e.t_upper for e in cert.events):
        raise AssertionError("q witness lies in a critical isolating bracket")
    return sum(e.t_upper < t for e in cert.events)


def test_all_65_exact_probability_design_problems(oracle12):
    model, direct = oracle12
    envelopes = {}
    for h_count in range(13):
        for threshold in range(1, 6):
            cert = certify_event_probability_envelope(model, h_count, threshold)
            assert cert.population == len(list(combinations(range(12), h_count)))
            assert len(cert.root_free_winners) == len(cert.events) + 1
            assert all(e.t_lower >= 0 and
                       e.t_upper - e.t_lower <= Fraction(1, 1_000_000)
                       for e in cert.events)
            envelopes[h_count, threshold] = cert

    assert len(envelopes) == 65
    assert sum(len(c.switches) for c in envelopes.values()) == 74
    assert sum(len(set(c.root_free_winners)) for c in envelopes.values()) == 136

    # Independent exact all-sequence, all-threshold q-witness comparisons.
    # These are regressions; whole-continuum completeness follows from the
    # exact root-count coverage in the checker, not from finite q sampling.
    for q in (Fraction(2), Fraction(10), Fraction(50)):
        for h_count in range(13):
            seqs = tuple(
                "".join("H" if i in indices else "P" for i in range(12))
                for indices in combinations(range(12), h_count)
            )
            for threshold in range(1, 6):
                cert = envelopes[h_count, threshold]
                values = {
                    s: _probability(direct[s], threshold, q)
                    for s in seqs
                }
                best = max(values.values())
                actual = tuple(s for s in seqs if values[s] == best)
                expected = cert.root_free_winners[_cell_at(cert, q)]
                assert actual == expected


def test_five_successive_globally_optimal_event_families(oracle12):
    model, direct = oracle12
    cert = certify_event_probability_envelope(model, 7, 4)
    assert cert.population == 792
    assert cert.coefficient_survivors == 14
    assert cert.skyline_classes == 7
    assert len(cert.final_classes) == 7
    assert len(cert.events) == 25
    assert len(cert.switches) == 4
    assert len(set(cert.root_free_winners)) == 5

    expected = (
        "HPHHPHHHPHPP",
        "HPHHPPHHPHPH",
        "HPHPHHHPPHPH",
        "HPHHHPHPPHPH",
        "HPHHPPHHHPPH",
    )
    assert tuple(dict.fromkeys(c[0] for c in cert.root_free_winners)) == expected
    for e, (lo, hi) in zip(cert.switches, (
        (Fraction(1), Fraction(2)), (Fraction(4), Fraction(5)),
        (Fraction(9), Fraction(10)), (Fraction(27), Fraction(28))
    )):
        assert lo < e.t_lower + 1 < e.t_upper + 1 < hi
        assert e.before != e.after

    # These are independently scored contact spectra from complete walks.
    assert direct[expected[0]] == {
        0: 15985, 1: 9936, 2: 3218, 3: 770, 4: 164
    }
    assert direct[expected[-1]] == {
        0: 16011, 1: 10200, 2: 3090, 3: 646, 4: 106, 5: 20
    }
    for seq in expected:
        for q in (Fraction(3, 2), Fraction(2), Fraction(5), Fraction(20)):
            direct_probability = _probability(direct[seq], 4, q)
            spectrum_probability = _probability(model.contact_spectrum(seq), 4, q)
            assert direct_probability == spectrum_probability


def test_probability_is_not_equivalent_to_partition_optimality(oracle12):
    model, direct = oracle12
    # At q=1 both objectives are comparable directly without any
    # temperature interpolation. Event probability optimum here is not
    # the partition-optimal all-q six-H sequence.
    cert = certify_event_probability_envelope(model, 6, 4)
    assert len(cert.switches) == 2
    assert len(set(cert.root_free_winners)) == 3
    assert "HPPHPHHPHPPH" not in cert.root_free_winners[0]
    assert "HPHHPPHHPHPP" in cert.root_free_winners[0]


class _Forged:
    length = 4

    def polynomial(self, sequence):
        return (9,)

    def contact_spectrum(self, sequence):
        return {1: 10} if sequence == "HPPP" else {0: 9}


def test_bad_spectra_and_work_budgets_fail_closed(oracle12):
    with pytest.raises(ValueError, match="inconsistent exact event"):
        certify_event_probability_envelope(_Forged(), 1, 1)
    model, _ = oracle12
    for bad in (-1, 13, 1.0, True):
        with pytest.raises(ValueError):
            certify_event_probability_envelope(model, bad, 2)
    for bad in (-1, 1.0, True):
        with pytest.raises(ValueError):
            certify_event_probability_envelope(model, 4, bad)
    for opts in (
        {"width": 0.1}, {"max_sequences": 10}, {"max_classes": 1},
        {"max_events": 1}, {"max_refinements": 0}
    ):
        with pytest.raises(ValueError):
            certify_event_probability_envelope(model, 7, 4, **opts)
    all_zero = certify_event_probability_envelope(model, 3, 5)
    assert all_zero.population == 220
    assert all_zero.skyline_classes == 1
    assert len(all_zero.root_free_winners[0]) == 220
    assert not all_zero.events
