"""Independent microstate enumeration and negative controls for local rewiring."""
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.cubic_tiling import CUBE, OFFSETS, ROOT, tile_paths
from foldbounds.local_rewire_confidence import (
    StructuralModel, check_structural_confidence, find_rewire_patches,
    hp_pair_weights, make_macro_target, minimum_selectivity_for_confidence,
)


def independent_complete_walks(n):
    """Separate brute-force 3D SAW generator (small sizes only)."""
    walks = [ROOT]
    for _ in range(n - 2):
        out = []
        for walk in walks:
            x, y, z = walk[-1]
            for dx, dy, dz in OFFSETS:
                pt = (x + dx, y + dy, z + dz)
                if pt not in walk:
                    out.append(walk + (pt,))
        walks = out
    return tuple(walks)


def independent_full_weight(seq, path, weights):
    """Score ALL nonbonded contacts without the production local edge code."""
    result = Fraction(1)
    for i, p in enumerate(path):
        for j in range(i + 2, len(path)):
            if sum(abs(a - b) for a, b in zip(p, path[j])) == 1:
                a, b = sorted((seq[i], seq[j]))
                result *= weights[(a, b)]
    return result


def independent_cube_motifs():
    out = []
    for start in CUBE:
        def grow(walk):
            if len(walk) == 8:
                out.append(walk)
                return
            for nxt in CUBE:
                if nxt in walk:
                    continue
                if sum(abs(a - b) for a, b in zip(walk[-1], nxt)) == 1:
                    grow(walk + (nxt,))
        grow((start,))
    return tuple(out)


def test_cube_library_against_independent_enumerator():
    assert set(independent_cube_motifs()) == set(tile_paths())
    assert len(tile_paths()) == 144
    assert len([p for p in tile_paths() if p[:2] == ROOT and
                p[-1] == (1, 1, 1)]) == 2


def test_macro_generator_is_independently_validated():
    for dims in ((2, 1, 1), (3, 1, 1), (2, 2, 2), (3, 3, 3),
                 (4, 4, 4)):
        path = make_macro_target(*dims)
        assert path[:2] == ROOT
        assert len(path) == 8 * dims[0] * dims[1] * dims[2]
        assert len(set(path)) == len(path)
        assert all(sum(abs(a - b) for a, b in zip(p, q)) == 1
                   for p, q in zip(path, path[1:]))
        assert set(path) == {
            (x, y, z)
            for x in range(2 * dims[0])
            for y in range(2 * dims[1])
            for z in range(2 * dims[2])
        }


def test_independent_exact_local_reweighting_oracle():
    motifs = [p for p in independent_cube_motifs()
              if p[:2] == ROOT and p[-1] == (1, 1, 1)]
    for seq in ("HHHHHHHH", "HPHPHPHP", "HHPPHPHP",
                "ABCDEFGH", "PPPPPPPP"):
        alphabet = sorted(set(seq))
        weights = {(a, b): Fraction(1 + (i + 2 * j) % 7, 3)
                   for i, a in enumerate(alphabet)
                   for j, b in enumerate(alphabet) if i <= j}
        for target in motifs:
            model = StructuralModel(seq, target, weights)
            patches = find_rewire_patches(model)
            assert len(patches) == 1
            patch = patches[0]
            assert set(patch.variants) == set(motifs)
            baseline = independent_full_weight(seq, target, weights)
            oracle_ratios = tuple(
                independent_full_weight(seq, v, weights) / baseline
                for v in patch.variants
            )
            assert patch.relative_weights == oracle_ratios
            assert patch.partition_factor == sum(oracle_ratios)
            assert patch.equal_weight_choices == sum(r == 1 for r in oracle_ratios)


def test_two_spatially_separated_patches_factorize_full_partition():
    path = make_macro_target(3, 1, 1)
    seq = "HPPHHHPHHPHPPHPHPHHPHPPH"
    assert len(seq) == len(path) == 24
    weights = {
        ("H", "H"): Fraction(3),
        ("H", "P"): Fraction(2),
        ("P", "P"): Fraction(1, 2),
    }
    model = StructuralModel(seq, path, weights)
    patches = find_rewire_patches(model)
    first = next(p for p in patches if p.start == 0)
    third = next(p for p in patches if p.start == 16)
    assert max(abs(x-y) for x, y in zip(first.cube_origin, third.cube_origin)) == 4
    baseline = independent_full_weight(seq, path, weights)
    total = Fraction(0)
    seen = set()
    for a, b in product(first.variants, third.variants):
        variant = a + path[8:16] + b
        assert variant[:2] == ROOT
        assert len(set(variant)) == len(variant)
        assert all(sum(abs(x-y) for x, y in zip(u, v)) == 1
                   for u, v in zip(variant, variant[1:]))
        assert variant not in seen
        seen.add(variant)
        total += independent_full_weight(seq, variant, weights) / baseline
    assert len(seen) == first.choices * third.choices
    assert total == first.partition_factor * third.partition_factor


def test_8mer_exact_full_ensemble_adversarial_upper_is_sound():
    walks = independent_complete_walks(8)
    assert len(walks) == 13565
    targets = [p for p in independent_cube_motifs()
               if p[:2] == ROOT and p[-1] == (1, 1, 1)]
    for seq in ("HHHHHHHH", "HPHPHPHP", "HPPHHPPH"):
        for q in (Fraction(1, 2), Fraction(1), Fraction(3)):
            weights = hp_pair_weights(q)
            full_partition = sum(
                (independent_full_weight(seq, p, weights) for p in walks),
                Fraction(0)
            )
            for target in targets:
                model = StructuralModel(
                    seq, target, weights, Fraction(1, 2), Fraction(2))
                cert = check_structural_confidence(model)
                nominal = independent_full_weight(seq, target, weights)
                # Exact extremum of independent [1/2,2] corrections:
                exact_max = 2 * nominal / (
                    2 * nominal + (full_partition - nominal) / 2
                )
                assert exact_max <= cert.microstate_probability_upper
                assert cert.rigid_fold_probability_upper >= cert.microstate_probability_upper


def test_all_h_temperature_universal_exact_isoenergy():
    path = make_macro_target(4, 4, 4)
    for q in (Fraction(1, 3), Fraction(1), Fraction(3),
              Fraction(1000)):
        model = StructuralModel("H" * len(path), path, hp_pair_weights(q))
        result = check_structural_confidence(model)
        assert result.constructive_conformations > 10**12
        assert result.exactly_isoenergetic_conformations == result.constructive_conformations
        assert result.relative_partition_lower == result.constructive_conformations
        assert result.microstate_probability_upper == Fraction(
            1, result.constructive_conformations)
        assert result.rigid_fold_probability_upper < Fraction(1, 10**9)
        assert result.rules_out_single_structure_majority
        assert result.rules_out_rigid_fold_majority


def test_periodic_heterogeneous_sequence_has_nontrivial_isoenergy_core():
    path = make_macro_target(4, 4, 4)
    seq = "".join("P" if i % 11 == 0 else "H" for i in range(len(path)))
    for q in (Fraction(1, 2), Fraction(1), Fraction(5),
              Fraction(100)):
        result = check_structural_confidence(StructuralModel(
            seq, path, hp_pair_weights(q), Fraction(1, 2), Fraction(2)))
        assert result.exactly_isoenergetic_conformations > 1000
        assert result.microstate_probability_upper < Fraction(1, 100)
        assert result.rigid_fold_probability_upper < Fraction(1, 2)
        assert result.relative_partition_lower >= result.exactly_isoenergetic_conformations


def test_arbitrary_sequence_and_pair_model_adversarial_bound():
    path = make_macro_target(3, 3, 3)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    seq = "".join(alphabet[i % len(alphabet)] for i in range(len(path)))
    weights = {
        (a, b): Fraction(2 + (i + 3*j) % 5, 4)
        for i, a in enumerate(sorted(alphabet))
        for j, b in enumerate(sorted(alphabet)) if i <= j
    }
    cert = check_structural_confidence(StructuralModel(
        seq, path, weights, Fraction(1, 2), Fraction(2)))
    assert cert.candidate_cubes >= 27
    assert cert.constructive_conformations > 1_000
    assert cert.relative_partition_lower > 1
    assert cert.microstate_probability_upper < Fraction(1, 2)
    assert cert.rigid_fold_probability_upper <= 1


def test_no_candidate_abstains_without_claiming_nonuniqueness():
    path = tuple((i, 0, 0) for i in range(24))
    model = StructuralModel("HP" * 12, path, hp_pair_weights(Fraction(3)))
    result = check_structural_confidence(model)
    assert not result.selected_cubes
    assert result.microstate_probability_upper == 1
    assert result.rigid_fold_probability_upper == 1
    assert not result.rules_out_single_structure_majority


def test_degeneracy_selectivity_threshold_is_exact():
    assert minimum_selectivity_for_confidence(100, Fraction(9, 10)) == 891
    assert minimum_selectivity_for_confidence(1, Fraction(9, 10)) == 0
    for states, p in ((0, Fraction(9, 10)),
                      (2.0, Fraction(9, 10)),
                      (10, 0.9), (10, Fraction(1))):
        with pytest.raises(ValueError):
            minimum_selectivity_for_confidence(states, p)


@pytest.mark.parametrize("modifier", [
    {"target": ((0, 0, 0), (1, 0, 0))},
    {"target": ((0, 0, 0), (1, 0, 0)) + ((1, 0, 0),)*22},
    {"target": ((0, 0, 0), (2, 0, 0)) + tuple((i, 0, 0) for i in range(2, 24))},
    {"target": [(i, 0, 0) for i in range(24)]},
    {"target": tuple((i, 0, float(i)) for i in range(24))},
    {"sequence": "H" * 23},
    {"sequence": "H" * 23 + "X"},
    {"pair_weights": {("H", "H"): Fraction(3)}},
    {"pair_weights": {("H", "H"): Fraction(3), ("H", "P"): 1,
                      ("P", "P"): Fraction(1)}},
    {"pair_weights": {("P", "H"): Fraction(3), ("H", "H"): Fraction(3),
                      ("P", "P"): Fraction(1)}},
    {"correction_lower": Fraction(0)},
    {"correction_lower": Fraction(3), "correction_upper": Fraction(2)},
    {"correction_upper": 2.0},
])
def test_rejects_invalid_untrusted_proposals(modifier):
    path = tuple((i, 0, 0) for i in range(24))
    base = StructuralModel("HP"*12, path, hp_pair_weights(Fraction(3)))
    with pytest.raises(ValueError):
        check_structural_confidence(replace(base, **modifier))


def test_hp_contact_weights_reject_approximate_and_negative_values():
    for invalid in (0, -1, 3, 3.0, True):
        with pytest.raises(ValueError):
            hp_pair_weights(invalid)
