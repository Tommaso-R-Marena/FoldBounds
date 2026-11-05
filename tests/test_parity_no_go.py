"""Independent cubic-lattice parity obstruction and exact-ensemble oracle."""
from fractions import Fraction
from itertools import product
import pytest

from foldbounds.cubic_tiling import ROOT, tile_paths
from foldbounds.local_rewire_confidence import make_macro_target
from foldbounds.parity_no_go import certify_hp_parity_no_go


MOVES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
         (0, -1, 0), (0, 0, 1), (0, 0, -1))


def exhaustive_rooted_paths(n):
    paths = [ROOT]
    for _ in range(n - 2):
        out = []
        for path in paths:
            for dx, dy, dz in MOVES:
                x, y, z = path[-1]
                p = (x + dx, y + dy, z + dz)
                if p not in path:
                    out.append(path + (p,))
        paths = out
    return tuple(paths)


def independent_contacts(path):
    return tuple(
        (i, j) for i in range(len(path))
        for j in range(i + 2, len(path))
        if sum(abs(a-b) for a, b in zip(path[i], path[j])) == 1
    )


def test_independent_every_contact_crosses_sequence_parity():
    for n in range(2, 9):
        for path in exhaustive_rooted_paths(n):
            assert all(i % 2 != j % 2 for i, j in independent_contacts(path))


def test_all_sequence_obstructions_against_full_3d_walk_oracle():
    for n in range(2, 8):
        paths = exhaustive_rooted_paths(n)
        contacts = [independent_contacts(path) for path in paths]
        for choices in product("HP", repeat=n):
            seq = "".join(choices)
            h_parities = {i % 2 for i, ch in enumerate(seq) if ch == "H"}
            if len(h_parities) > 1:
                continue
            assert all(
                not any(seq[i] == seq[j] == "H" for i,j in edges)
                for edges in contacts
            )
            target = paths[0]
            certificate = certify_hp_parity_no_go(
                seq, target, q=Fraction(3),
                correction_lower=Fraction(1, 2),
                correction_upper=Fraction(2))
            assert certificate.all_conformations_have_zero_hh_contacts
            assert certificate.monotone_walks_lower == 3**(n-2)
            assert len(paths) >= certificate.monotone_walks_lower


def test_exact_eight_residue_partition_contains_constructive_lower():
    paths = exhaustive_rooted_paths(8)
    assert len(paths) == 13565
    chosen = next(path for path in tile_paths() if path[:2] == ROOT)
    seq = "HP" * 4
    cert = certify_hp_parity_no_go(
        seq, chosen, q=Fraction(100),
        correction_lower=Fraction(1, 2), correction_upper=Fraction(2))
    assert cert.monotone_walks_lower == 729
    exact_worst_corrected = Fraction(2) / (
        Fraction(2) + Fraction(1,2)*(len(paths)-1)
    )
    assert exact_worst_corrected <= cert.probability_microstate_upper
    assert cert.probability_microstate_upper == Fraction(4, 4 + 728)


def test_positive_monotone_walks_are_injective_and_self_avoiding():
    for n in range(2, 8):
        paths = set()
        for continuation in product(((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                                    repeat=n-2):
            path = list(ROOT)
            for delta in continuation:
                path.append(tuple(a+b for a, b in zip(path[-1], delta)))
            path = tuple(path)
            assert len(path) == len(set(path))
            assert path[:2] == ROOT
            paths.add(path)
        assert len(paths) == 3**(n-2)


def test_long_balanced_heterogeneous_parity_flatlandscape():
    path = make_macro_target(3, 3, 3)
    n = len(path)
    assert n == 216
    seq = "HP" * (n // 2)
    for q in (Fraction(1, 5), Fraction(1), Fraction(4),
              Fraction(1000000)):
        cert = certify_hp_parity_no_go(
            seq, path, q=q,
            correction_lower=Fraction(1, 2),
            correction_upper=Fraction(2)
        )
        assert cert.hydrophobic_residues == 108
        assert cert.compatible_h_index_parity == 0
        assert cert.monotone_walks_lower == 3**214
        assert cert.probability_microstate_upper < Fraction(1, 10**100)
        assert cert.probability_rigid_fold_upper < Fraction(1, 10**98)
        assert cert.rules_out_exact_state_majority


def test_all_polar_model_is_equally_flat_and_odd_parity_h_supported():
    path = make_macro_target(2, 2, 2)
    for seq, hp, parity in (
        ("P" * 64, 0, None),
        ("PH" * 32, 32, 1),
    ):
        cert = certify_hp_parity_no_go(seq, path, q=Fraction(2))
        assert cert.hydrophobic_residues == hp
        assert cert.compatible_h_index_parity == parity
        assert cert.probability_microstate_upper == Fraction(1, 3**62)


def test_reject_mixed_parities_and_non_hp_or_invalid_model():
    path = make_macro_target(2, 2, 2)
    for seq in ("HH" + "P"*62, "HPH" + "H" + "P"*60):
        with pytest.raises(ValueError, match="both"):
            certify_hp_parity_no_go(seq, path, q=Fraction(2))
    with pytest.raises(ValueError):
        certify_hp_parity_no_go("A"*64, path, q=Fraction(2))
    with pytest.raises(ValueError):
        certify_hp_parity_no_go("HP"*32, path, q=2)
    with pytest.raises(ValueError):
        certify_hp_parity_no_go("HP"*32, path, q=Fraction(0))
    with pytest.raises(ValueError):
        certify_hp_parity_no_go("HP"*32, path[:-1], q=Fraction(2))
    with pytest.raises(ValueError):
        certify_hp_parity_no_go(
            "HP"*32, path, q=Fraction(2), correction_lower=Fraction(0))
