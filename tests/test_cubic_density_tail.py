"""Audit of uniform infinite-size 3D heterogeneous HP contact-density gate."""
from dataclasses import replace
from fractions import Fraction
import pytest

from foldbounds.cubic_density_tail import certify_periodic_cubic_density_tail
from foldbounds.cubic_tiling import (
    CUBE, ROOT, TiledCertificate, check_tiled_family, macro_snake, tile_paths
)


def construct_one_path_through_cube_grid(side):
    """Independently instantiate one complete member of the symbolic family."""
    macro = macro_snake(side, side, side)
    motifs = tile_paths()
    entrance = (0, 0, 0)
    full = []
    for i, cell in enumerate(macro):
        if i == 0:
            # First macro tile must exit through +x for side>=2.
            exit_site = (1, 1, 1)
            candidates = [p for p in motifs if p[:2] == ROOT and p[-1] == exit_site]
        elif i == len(macro) - 1:
            candidates = [p for p in motifs if p[0] == entrance]
        else:
            difference = tuple(b - a for a, b in zip(cell, macro[i + 1]))
            axis = next(j for j in range(3) if difference[j])
            outward = 1 if difference[axis] > 0 else 0
            exit_site = min(p for p in CUBE if sum(p) % 2 == 1 and p[axis] == outward)
            candidates = [
                p for p in motifs if p[0] == entrance and p[-1] == exit_site
            ]
        assert candidates, f"no compatible path in macro tile {i}"
        path = min(candidates)
        full.extend((2 * cell[0] + x, 2 * cell[1] + y, 2 * cell[2] + z)
                    for x, y, z in path)
        if i < len(macro) - 1:
            delta = tuple(b - a for a, b in zip(cell, macro[i + 1]))
            if i == 0:
                exit_site = path[-1]
            entrance = tuple(exit_site[j] - delta[j] for j in range(3))
            assert sum(entrance) % 2 == 0
    return tuple(full)


def test_explicit_global_3d_hamiltonian_skeleton_parity_and_bridges():
    for m in (2, 3, 4):
        p = construct_one_path_through_cube_grid(m)
        assert len(p) == 8 * m**3
        assert p[:2] == ROOT
        assert len(set(p)) == len(p)
        assert all(sum(abs(x - y) for x, y in zip(a, b)) == 1
                   for a, b in zip(p, p[1:]))
        assert set(p) == {
            (x, y, z) for x in range(2 * m)
            for y in range(2 * m) for z in range(2 * m)
        }


def test_infinite_3d_heterogeneous_tail_q5_and_q4():
    first = certify_periodic_cubic_density_tail(11, q=5)
    second = certify_periodic_cubic_density_tail(22, q=4)
    assert first.residues == 10648
    assert first.polar_residues == 968
    assert first.event_threshold == 5324
    assert first.hh_contacts_lower == 15972
    assert second.residues == 85184
    assert second.polar_residues == 7744
    assert second.event_threshold == 42592
    assert second.hh_contacts_lower == 133584
    for result in (first, second):
        assert result.c8 == 308981
        assert result.step == 11
        assert result.certified
        assert result.probability_lower > Fraction(99, 100)
        assert result.family_count_lower < result.all_walk_count_upper
    assert first.failure_odds_upper < Fraction(1, 10**900)
    assert second.failure_odds_upper < Fraction(1, 10**2700)


def test_exact_recurrence_nonincrease_for_further_infinite_tail_members():
    for q, m in ((5, 11), (4, 22)):
        initial = certify_periodic_cubic_density_tail(m, q=q)
        next_case = certify_periodic_cubic_density_tail(m + 11, q=q)
        assert next_case.failure_odds_upper < initial.failure_odds_upper
        assert next_case.probability_lower > initial.probability_lower


def test_large_finite_transfer_rechecks_tail_at_first_side():
    m = 11
    tail = certify_periodic_cubic_density_tail(m, q=5)
    seq = "".join(
        "P" if i % 11 == 0 else "H" for i in range(tail.residues)
    )
    finite = check_tiled_family(TiledCertificate(
        (m, m, m), seq, Fraction(5), tail.event_threshold,
        Fraction(1, 2), Fraction(2)
    ))
    assert finite.family_paths >= tail.family_count_lower
    assert finite.family_min_contacts >= tail.hh_contacts_lower
    assert finite.total_walks_upper <= tail.all_walk_count_upper
    assert finite.probability_lower >= tail.probability_lower


def test_forged_parameters_and_unfavorable_couplings_rejected():
    for bad_m in (0, 10, 12, 21, True, 11.0):
        with pytest.raises(ValueError):
            certify_periodic_cubic_density_tail(bad_m, q=5)
    for bad_q in (1, Fraction(5), 2.0, True):
        with pytest.raises(ValueError):
            certify_periodic_cubic_density_tail(11, q=bad_q)
    for lo, hi in (
        (Fraction(0), Fraction(2)),
        (Fraction(3), Fraction(2)),
        (0.5, Fraction(2)),
        (Fraction(1), 2.0),
    ):
        with pytest.raises(ValueError):
            certify_periodic_cubic_density_tail(
                11, q=5, correction_lower=lo, correction_upper=hi
            )
    with pytest.raises(ValueError):
        certify_periodic_cubic_density_tail(
            11, q=5, goal_probability=1.0
        )
    with pytest.raises(ValueError, match="monotonicity"):
        certify_periodic_cubic_density_tail(11, q=4)
    with pytest.raises(ValueError, match="monotonicity"):
        certify_periodic_cubic_density_tail(22, q=3)


def test_larger_adversarial_uncertainty_preserves_exactness():
    initial = certify_periodic_cubic_density_tail(11, q=5)
    changed = certify_periodic_cubic_density_tail(
        11, q=5, correction_lower=Fraction(1, 3),
        correction_upper=Fraction(3)
    )
    assert changed.failure_odds_upper == initial.failure_odds_upper * Fraction(9, 4)
    assert changed.probability_lower < initial.probability_lower
    assert changed.probability_lower > Fraction(99, 100)
