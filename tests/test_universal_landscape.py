"""Independent exact-oracle tests of universal partition-polynomial transforms."""
from dataclasses import replace
from fractions import Fraction
from itertools import combinations
from math import comb

import pytest

from foldbounds.contact_map_atlas import ContactMapAtlas, build_atlas
from foldbounds.universal_landscape import (
    UniversalLandscape, audit_landscape, build_landscape,
)


def as_sequence(mask, n):
    return "".join("H" if mask & (1 << i) else "P" for i in range(n))


def independent_partitions(atlas, seq, q):
    """Score each geometry directly rather than using the polynomial transform."""
    score_counts = {}
    for mask, count in atlas.multiplicities:
        k = sum(1 for bit, (i, j) in enumerate(atlas.pairs)
                if mask & (1 << bit) and seq[i] == seq[j] == 'H')
        score_counts[k] = score_counts.get(k, 0) + count
    z = sum((count * q**k for k, count in score_counts.items()), Fraction(0))
    return z, score_counts


@pytest.mark.parametrize("n", range(2, 9))
def test_all_sequences_all_q_against_independent_geometry(n):
    atlas = build_atlas(n)
    model = build_landscape(atlas)
    for mask in range(1 << n):
        s = as_sequence(mask, n)
        for q in (Fraction(1), Fraction(3, 2), Fraction(3)):
            independent_z, counts = independent_partitions(atlas, s, q)
            assert model.partition(s, q) == independent_z
            assert model.contact_spectrum(s) == {k: counts[k] for k in sorted(counts) if counts[k]}
            z = independent_z
            mean = sum((k * count * q**k for k, count in counts.items()), Fraction(0)) / z
            variance = sum((k*k * count * q**k for k, count in counts.items()), Fraction(0)) / z - mean**2
            assert model.contact_moments(s, q) == (mean, variance)
        assert model.polynomial(s)[0] == atlas.conformations
        expected = tuple(sum((count * comb(k, j) for k, count in counts.items() if k >= j), 0)
                         for j in range(max(counts) + 1))
        assert model.polynomial(s) == expected


@pytest.mark.parametrize("n", [9, 10, 11, 12])
def test_larger_atlases_deterministic_selection(n):
    atlas = build_atlas(n)
    model = build_landscape(atlas)
    for index in range(0, 1 << n, max(1, (1 << n) // 31)):
        s = as_sequence(index, n)
        assert model.contact_spectrum(s) == atlas.spectrum(s)
        assert model.partition(s, Fraction(5, 2)) == atlas.partition(s, Fraction(5, 2))
        assert model.contact_moments(s, Fraction(2)) == atlas.contact_moments(s, Fraction(2))
    assert model.contact_spectrum('H' * n) == atlas.spectrum('H' * n)
    assert model.contact_spectrum('P' * n) == {0: atlas.conformations}
    assert len(model.all_partitions(Fraction(1))) == (1 << n)
    assert set(model.all_partitions(Fraction(1))) == {Fraction(atlas.conformations)}


def test_mixed_derivative_positivity_and_inclusion_exclusion():
    n = 7
    model = build_landscape(build_atlas(n))
    q = Fraction(7, 3)
    for base in range(1 << n):
        available = [i for i in range(n) if not base & (1 << i)]
        s = as_sequence(base, n)
        for k in range(min(3, len(available)) + 1):
            for addition in combinations(available, k):
                diff = model.mixed_difference(s, addition, q)
                assert diff >= 0
                independent = Fraction(0)
                for count in range(k + 1):
                    for chosen in combinations(addition, count):
                        mutant = base | sum(1 << i for i in chosen)
                        independent += (-1)**(k - count) * model.partition(as_sequence(mutant, n), q)
                assert diff == independent


def test_first_order_contact_interactions_recover_edge_frequency():
    atlas = build_atlas(10)
    landscape = build_landscape(atlas)
    for bit, (i, j) in enumerate(atlas.pairs):
        pair = (1 << i) | (1 << j)
        expected = sum(m for contact_mask, m in atlas.multiplicities if contact_mask & (1 << bit))
        assert landscape.interactions[pair][1] == expected
    for subset, poly in enumerate(landscape.interactions):
        if subset != 0:
            assert poly[0] == 0
        if subset.bit_count() <= 1:
            assert not any(poly[1:])
        assert all(c >= 0 for c in poly)
    assert landscape.interaction('P'*10, Fraction(2)) == atlas.conformations


def test_known_hhhh_coefficients_and_epistasis():
    model = build_landscape(build_atlas(4))
    assert model.polynomial('HHHH') == (9, 2)
    assert model.polynomial('PPPP') == (9,)
    assert model.interaction('HPPH', Fraction(5, 2)) == 3
    assert model.mixed_difference('PPPP', (0, 3), Fraction(5, 2)) == 3
    assert model.mixed_difference('PPPH', (0,), Fraction(5, 2)) == 3
    assert model.mixed_difference('PPPP', (1, 2), Fraction(5, 2)) == 0
    assert model.contact_moments('HHHH', Fraction(2)) == (Fraction(4, 11), Fraction(28, 121))
    assert model.partition('HHHH', Fraction(2)) == 11


def test_fail_closed_invalid_inputs_budgets_and_forged_atlas():
    atlas = build_atlas(6)
    model = build_landscape(atlas)
    for q in (Fraction(0), Fraction(1, 2), 2, 2., True):
        with pytest.raises(ValueError):
            model.partition('HHHHHH', q)
    for s in ('', 'HHH', 'XXXXXX', 'hhhhhh', 'HPHHPX'):
        with pytest.raises(ValueError):
            model.partition(s, Fraction(2))
    for additions in ((1, 1), (-1,), (0,), (99,), [1], (1.0,), (True,)):
        with pytest.raises(ValueError):
            model.mixed_difference('HPPPPP', additions, Fraction(2))
    with pytest.raises(ValueError):
        build_landscape(atlas, max_edge_subsets=1)
    with pytest.raises(ValueError):
        build_landscape(atlas, max_length=5)
    with pytest.raises(ValueError):
        build_landscape(atlas, require_audit=1)
    forged = ContactMapAtlas(4, ((0, 1),), 1)
    with pytest.raises(ValueError, match='independent'):
        build_landscape(forged)
    with pytest.raises(ValueError):
        UniversalLandscape(4, 1, ((1,),), ((1,),))


def test_full_independent_polynomial_and_mobius_audit():
    atlas = build_atlas(10)
    landscape = build_landscape(atlas)
    assert audit_landscape(landscape, atlas)
    forged = [list(row) for row in landscape.partitions]
    forged[3][1] += 1
    contaminated = replace(landscape, partitions=tuple(tuple(row) for row in forged))
    assert not audit_landscape(contaminated, atlas)
    with pytest.raises(ValueError):
        audit_landscape(landscape, atlas, max_mask_checks=10)
    with pytest.raises(ValueError):
        audit_landscape(landscape, atlas, max_length=9)
