"""Independent exact-oracle and adversarial tests of reusable geometry atlases."""
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.contact_map_atlas import (
    ContactMapAtlas, audit_atlas, build_atlas, contact_pairs,
)


def brute_spectrum(sequence):
    # Independent walk generator and contact score (no production atlas helpers).
    paths = [((0, 0), (1, 0))]
    for _ in range(2, len(sequence)):
        next_paths = []
        for path in paths:
            x, y = path[-1]
            for dx, dy in ((0, -1), (-1, 0), (0, 1), (1, 0)):
                q = (x + dx, y + dy)
                if q not in path:
                    next_paths.append(path + (q,))
        paths = next_paths
    result = {}
    for path in paths:
        score = 0
        for j in range(len(sequence)):
            for i in range(j - 1):
                if sequence[i] == sequence[j] == "H" and j > i + 1:
                    score += (abs(path[i][0] - path[j][0])
                              + abs(path[i][1] - path[j][1]) == 1)
        result[score] = result.get(score, 0) + 1
    return dict(sorted(result.items()))


@pytest.mark.parametrize("n", range(2, 9))
def test_all_hp_sequences_against_independent_enumeration(n):
    atlas = build_atlas(n)
    assert audit_atlas(atlas)
    for letters in product("HP", repeat=n):
        s = "".join(letters)
        expected = brute_spectrum(s)
        assert atlas.spectrum(s) == expected
        for q in (Fraction(1), Fraction(3, 2), Fraction(3)):
            z = sum((count * q ** k for k, count in expected.items()), Fraction(0))
            assert atlas.partition(s, q) == z
            for threshold in (0, 1, 2, 3):
                mass = sum((count * q ** k for k, count in expected.items()
                            if k >= threshold), Fraction(0))
                assert atlas.event_probability(s, q, threshold) == mass / z


@pytest.mark.parametrize("n,walks", [(2, 1), (3, 3), (4, 9), (5, 25),
                                      (6, 71), (7, 195), (8, 543), (9, 1479)])
def test_known_anchored_saw_counts(n, walks):
    atlas = build_atlas(n)
    assert atlas.conformations == walks
    assert atlas.spectrum("P" * n) == {0: walks}
    assert atlas.multiplicities[0][1] % 2 == 1
    assert all(weight % 2 == 0 for mask, weight in atlas.multiplicities if mask)


def test_hhhh_and_exact_moments():
    atlas = build_atlas(4)
    assert atlas.spectrum("HHHH") == {0: 7, 1: 2}
    assert atlas.partition("HHHH", Fraction(2)) == 11
    assert atlas.event_probability("HHHH", Fraction(2), 1) == Fraction(4, 11)
    assert atlas.contact_moments("HHHH", Fraction(2)) == (Fraction(4, 11), Fraction(28, 121))
    assert atlas.contact_moments("PPPP", Fraction(3)) == (Fraction(0), Fraction(0))
    assert atlas.score_covariance("HHHH", "HHHH", Fraction(2)) == Fraction(28, 121)
    assert atlas.score_covariance("HHHH", "PPPP", Fraction(2)) == 0


def test_derivative_identity_and_mutation_scan():
    atlas = build_atlas(9)
    seq = "HHPPHHPHH"
    spectrum = atlas.spectrum(seq)
    q = Fraction(5, 2)
    z = sum((a * q ** k for k, a in spectrum.items()), Fraction(0))
    zp = sum((k * a * q ** (k - 1) for k, a in spectrum.items() if k), Fraction(0))
    zpp = sum((k * (k - 1) * a * q ** (k - 2) for k, a in spectrum.items() if k >= 2), Fraction(0))
    mean, variance = atlas.contact_moments(seq, q)
    assert mean == q * zp / z
    assert variance == q * zp / z + q * q * zpp / z - mean * mean
    assert variance >= 0
    for i in range(len(seq)):
        changed = seq[:i] + ("P" if seq[i] == "H" else "H") + seq[i + 1:]
        assert atlas.spectrum(changed) == brute_spectrum(changed)


def test_fail_closed_inputs_and_incomplete_atlas():
    for n in (0, 1, 2.0, True):
        with pytest.raises(ValueError):
            contact_pairs(n)
    with pytest.raises(ValueError, match="max_nodes"):
        build_atlas(9, max_nodes=3)
    atlas = build_atlas(5)
    for seq in ("", "HHH", "XXXXX", "hhhhh", "HPHH1"):
        with pytest.raises(ValueError):
            atlas.spectrum(seq)
    for q in (1, 2., Fraction(0), Fraction(1, 2), True):
        with pytest.raises(ValueError):
            atlas.partition("HHHHH", q)
    with pytest.raises(ValueError):
        atlas.event_probability("HHHHH", Fraction(2), 1.0)
    with pytest.raises(ValueError):
        audit_atlas(atlas, max_length=4)
    with pytest.raises(ValueError):
        ContactMapAtlas(5, ((0, 2),), 2)  # violates unique straight fixed point
    with pytest.raises(ValueError):
        ContactMapAtlas(5, ((0, 1), (0, 2)), 3)  # duplicate mask
    with pytest.raises(ValueError):
        ContactMapAtlas(5, ((0, 1), (1 << len(contact_pairs(5)), 2)), 3)
    # Syntactically consistent, but forged: independent audit detects it.
    forged = ContactMapAtlas(4, ((0, 1),), 1)
    assert audit_atlas(forged) is False
