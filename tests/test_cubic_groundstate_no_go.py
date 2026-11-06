"""Independent tests for 3D lattice isoperimetric ground-state ambiguity."""
from fractions import Fraction

import pytest

from foldbounds.cubic_groundstate_no_go import certify_cubic_groundstate_no_go
from foldbounds.local_rewire_confidence import make_macro_target


def brute_force_anchored_walks(n):
    frontier = [((0, 0, 0), (1, 0, 0))]
    for _ in range(n - 2):
        nxt = []
        for path in frontier:
            x, y, z = path[-1]
            for dx, dy, dz in ((1, 0, 0), (-1, 0, 0),
                               (0, 1, 0), (0, -1, 0),
                               (0, 0, 1), (0, 0, -1)):
                p = (x + dx, y + dy, z + dz)
                if p not in path:
                    nxt.append(path + (p,))
        frontier = nxt
    return tuple(frontier)


def independent_contact_count(path):
    return sum(
        sum(abs(x-y) for x, y in zip(path[i], path[j])) == 1
        for i in range(len(path))
        for j in range(i + 2, len(path))
    )


def test_isoperimetric_cap_is_sharp_on_complete_eight_residue_oracle():
    all_walks = brute_force_anchored_walks(8)
    assert len(all_walks) == 13565
    scores = [independent_contact_count(p) for p in all_walks]
    assert max(scores) == 2 * 8 - 3 * 2**2 + 1 == 5
    assert min(scores) == 0


def test_global_contact_optimality_of_full_cubic_hamiltonian_targets():
    for side in (4, 6, 8):
        result = certify_cubic_groundstate_no_go(side)
        target = make_macro_target(side // 2, side // 2, side // 2)
        assert result.residues == side**3 == len(target)
        assert independent_contact_count(target) == result.target_contacts
        assert result.global_contact_upper == result.target_contacts
        assert result.global_contact_upper == 2 * side**3 - 3 * side**2 + 1
        assert result.globally_optimal_for_attractive_hh
        assert result.confidence.exactly_isoenergetic_conformations > 1
        assert result.confidence.microstate_probability_upper < 1


def test_ground_state_cannot_be_a_unique_fold_even_at_strong_coupling():
    for side in (6, 8):
        result = certify_cubic_groundstate_no_go(
            side, q=Fraction(1000),
            correction_lower=Fraction(1, 2),
            correction_upper=Fraction(2))
        assert result.globally_optimal_for_attractive_hh
        assert result.confidence.exactly_isoenergetic_conformations > 10**4
        assert result.confidence.rules_out_single_structure_majority
        assert result.confidence.rules_out_rigid_fold_majority
        assert result.confidence.microstate_probability_upper < Fraction(1, 1000)


def test_energy_degeneracy_independent_of_contact_weight():
    for q in (Fraction(1, 10), Fraction(1), Fraction(2), Fraction(100)):
        result = certify_cubic_groundstate_no_go(6, q=q)
        assert result.globally_optimal_for_attractive_hh == (q > 1)
        assert result.confidence.relative_partition_lower == (
            result.confidence.exactly_isoenergetic_conformations
        )
    assert certify_cubic_groundstate_no_go(
        6, q=Fraction(1)).globally_optimal_for_attractive_hh is False


@pytest.mark.parametrize("bad_side", (0, 1, 2, 3, 5, -4, True, 4.0))
def test_bad_side_rejected(bad_side):
    with pytest.raises(ValueError):
        certify_cubic_groundstate_no_go(bad_side)


@pytest.mark.parametrize("bad_q", (0, -1, 4, 4.0, True))
def test_bad_contact_weight_rejected(bad_q):
    with pytest.raises(ValueError):
        certify_cubic_groundstate_no_go(4, q=bad_q)


def test_adversarial_corrections_can_only_weaken_ceiling():
    basic = certify_cubic_groundstate_no_go(6, q=Fraction(4))
    robust = certify_cubic_groundstate_no_go(
        6, q=Fraction(4),
        correction_lower=Fraction(1, 2),
        correction_upper=Fraction(2))
    assert robust.confidence.microstate_probability_upper > (
        basic.confidence.microstate_probability_upper)
    assert robust.confidence.rigid_fold_probability_upper > (
        basic.confidence.rigid_fold_probability_upper)
