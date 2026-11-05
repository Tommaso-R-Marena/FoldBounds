"""Independent oracle and fail-closed regression for 3D tiled HP certificates."""
from dataclasses import replace
from fractions import Fraction
from itertools import product
import pytest

from foldbounds.cubic_tiling import (
    CUBE, OFFSETS, ROOT, TiledCertificate, all_saw_count_upper_3d,
    certify_majority_halfline, check_tiled_family, macro_snake,
    minimum_certified_integer_q, tile_paths, _free_saw_counts_3d,
)


def independent_saws(n):
    walks = [ROOT]
    for _ in range(n - 2):
        next_walks = []
        for walk in walks:
            x, y, z = walk[-1]
            for dx, dy, dz in OFFSETS:
                new = (x + dx, y + dy, z + dz)
                if new not in walk:
                    next_walks.append(walk + (new,))
        walks = next_walks
    return tuple(walks)


def independent_hh(seq, walk):
    return sum(
        seq[i] == seq[j] == "H"
        and sum(abs(a - b) for a, b in zip(walk[i], walk[j])) == 1
        for i in range(len(walk))
        for j in range(i + 2, len(walk))
    )


def independent_cube_motifs():
    """Second generator using all eight candidate vertices, not production adjacency."""
    motifs = []
    for start in CUBE:
        def grow(walk):
            if len(walk) == 8:
                motifs.append(walk)
                return
            for candidate in CUBE:
                if candidate not in walk and sum(
                    abs(a - b) for a, b in zip(walk[-1], candidate)
                ) == 1:
                    grow(walk + (candidate,))
        grow((start,))
    return tuple(motifs)


def independent_two_cube_family():
    motifs = independent_cube_motifs()
    combined = []
    for first in motifs:
        if first[:2] != ROOT:
            continue
        for second in motifs:
            walk = first + tuple((x + 2, y, z) for x, y, z in second)
            if sum(abs(a - b) for a, b in zip(walk[7], walk[8])) == 1:
                combined.append(walk)
    return tuple(combined)


def test_all_8_site_motifs_independently_enumerated():
    oracle = independent_cube_motifs()
    assert set(tile_paths()) == set(oracle)
    assert len(tile_paths()) == 144
    assert sum(path[:2] == ROOT for path in tile_paths()) == 6


def test_every_macro_route_covers_disjoint_tiles_with_unit_bridges():
    for dimensions in (
        (1, 1, 1), (2, 1, 1), (1, 2, 3),
        (2, 2, 2), (3, 4, 2), (4, 4, 4)
    ):
        route = macro_snake(*dimensions)
        assert len(route) == dimensions[0] * dimensions[1] * dimensions[2]
        assert route[0] == (0, 0, 0)
        assert len(set(route)) == len(route)
        assert all(
            sum(abs(x - y) for x, y in zip(a, b)) == 1
            for a, b in zip(route, route[1:])
        )


def test_independent_all_8_site_3d_conformations_and_adversarial_extremes():
    walks = independent_saws(8)
    assert len(walks) <= all_saw_count_upper_3d(8)
    for seq in ("H" * 8, "P" * 8, "HPHPHPHP", "HHPHHHPH", "PHHHHHHH"):
        scores = [independent_hh(seq, walk) for walk in walks]
        for event_threshold in (1, 2, 3, 5):
            for q in (Fraction(1), Fraction(3), Fraction(7)):
                e = sum((q ** k for k in scores if k >= event_threshold), Fraction(0))
                n = sum((q ** k for k in scores if k < event_threshold), Fraction(0))
                result = check_tiled_family(TiledCertificate(
                    (1, 1, 1), seq, q, event_threshold, Fraction(1, 2), Fraction(2)
                ))
                lower = (e / 2) / (e / 2 + 2 * n) if e else Fraction(0)
                assert result.probability_lower <= lower
                assert result.event_mass_lower <= e / 2
                assert result.complement_mass_upper >= 2 * n


def test_two_cube_transfer_count_and_heterogeneous_contacts():
    family = independent_two_cube_family()
    assert len(family) == 36
    assert len(set(family)) == 36
    for seq in ("H" * 16, "P" * 16, "HP" * 8,
                "HPHHPHHPPHHPHHHP", "H" * 15 + "P"):
        for threshold in (1, 4, 8):
            result = check_tiled_family(
                TiledCertificate((2, 1, 1), seq, Fraction(3), threshold)
            )
            assert result.family_paths == len(family)
            assert result.family_min_contacts <= min(
                independent_hh(seq, walk) for walk in family
            )
            assert result.total_all_h_contacts == 13
            assert result.family_paths <= result.total_walks_upper


def test_known_short_free_3d_saw_counts():
    assert _free_saw_counts_3d(8) == (
        1, 5, 25, 121, 589, 2821, 13565, 64661, 308981
    )
    # Independent count comparison, without importing a table as ground truth.
    for n in range(2, 9):
        assert len(independent_saws(n)) <= all_saw_count_upper_3d(n)


def test_216_residue_3d_sequence_with_20_polar_residues():
    n = 216
    seq = "".join("P" if i % 11 == 0 else "H" for i in range(n))
    assert seq.count("P") == 20
    cert = TiledCertificate(
        (3, 3, 3), seq, Fraction(5), 54, Fraction(1, 2), Fraction(2)
    )
    result = certify_majority_halfline(cert)
    assert result.family_paths == 76973654136380749302988800
    assert result.family_min_contacts == 244
    assert result.total_all_h_contacts == 325
    assert minimum_certified_integer_q(cert, max_q=12) == 5
    assert not check_tiled_family(replace(cert, q=Fraction(4))).guaranteed_event_majority


def test_3d_cube_scaling_and_512_residue_heterogeneous_model():
    for dims in ((3, 3, 3), (4, 4, 4)):
        n = 8 * dims[0] * dims[1] * dims[2]
        cert = TiledCertificate(
            dims, "H" * n, Fraction(4), n // 2, Fraction(1, 2), Fraction(2)
        )
        result = certify_majority_halfline(cert)
        assert result.family_min_contacts == result.total_all_h_contacts
        assert result.family_paths > 1_000_000
        assert result.probability_lower > Fraction(99, 100)
    n = 512
    seq = "".join("P" if i % 11 == 0 else "H" for i in range(n))
    cert = TiledCertificate(
        (4, 4, 4), seq, Fraction(6), n // 2, Fraction(1, 2), Fraction(2)
    )
    result = certify_majority_halfline(cert)
    assert seq.count("P") == 47
    assert result.family_min_contacts == 644
    assert minimum_certified_integer_q(cert, max_q=8) == 6
    assert not check_tiled_family(replace(cert, q=Fraction(5))).guaranteed_event_majority


def test_contact_lower_valid_for_all_small_cube_hp_sequences():
    complete = [path for path in independent_cube_motifs() if path[:2] == ROOT]
    for letters in product("HP", repeat=8):
        seq = "".join(letters)
        result = check_tiled_family(
            TiledCertificate((1, 1, 1), seq, Fraction(2), 2)
        )
        assert all(
            independent_hh(seq, path) >= result.coordination_lower_contacts
            for path in complete
        )


def test_fail_closed_and_invalid_parameters():
    base = TiledCertificate((2, 2, 2), "H" * 64, Fraction(3), 16)
    invalid = (
        replace(base, sequence="H" * 63),
        replace(base, sequence="H" * 63 + "X"),
        replace(base, q=3),
        replace(base, q=Fraction(0)),
        replace(base, event_at_least=0),
        replace(base, event_at_least=True),
        replace(base, correction_lower=Fraction(0)),
        replace(base, correction_lower=Fraction(3), correction_upper=Fraction(2)),
        replace(base, macro_dimensions=(2, 2)),
        replace(base, macro_dimensions=[2, 2, 2]),
    )
    for case in invalid:
        with pytest.raises(ValueError):
            check_tiled_family(case)
    for dims in ((2, 2, 0), (0, 2, 2), (2, True, 2), (2, 3.0, 2)):
        with pytest.raises(ValueError):
            check_tiled_family(replace(base, macro_dimensions=dims))
    with pytest.raises(ValueError):
        minimum_certified_integer_q(base, max_q=0)
    with pytest.raises(ValueError):
        all_saw_count_upper_3d(1)


def test_exponential_family_improves_same_k_single_witness_baseline():
    cert = TiledCertificate(
        (3, 3, 3), "H" * 216, Fraction(4), 108, Fraction(1, 2), Fraction(2)
    )
    result = check_tiled_family(cert)
    single_event = cert.correction_lower * cert.q ** result.family_min_contacts
    baseline = single_event / (single_event + result.complement_mass_upper)
    assert result.family_paths > 10**25
    assert result.probability_lower > baseline
    assert minimum_certified_integer_q(cert, max_q=6) == 4
