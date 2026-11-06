"""Independent adversarial tests for rigorous infinite-square-size HP certificates."""
from fractions import Fraction

import pytest

from foldbounds.square_density import (
    _exact_anchored_saw_counts16,
    _majorant,
    certify_square_density_tail,
    check_square_witness,
    exact_anchored_saw16,
    find_first_square_tail,
    square_serpentine,
)


def external_saw_paths(n):
    """Independent conformations, not the production depth-first count code."""
    paths = [((0, 0), (1, 0))]
    for _ in range(n - 2):
        next_paths = []
        for path in paths:
            tip_x, tip_y = path[-1]
            for step in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                position = (tip_x + step[0], tip_y + step[1])
                if position not in path:
                    next_paths.append(path + (position,))
        paths = next_paths
    return paths


def external_hh_contacts(path):
    return sum(abs(path[i][0]-path[j][0]) + abs(path[i][1]-path[j][1]) == 1
               for i in range(len(path)) for j in range(i + 2, len(path)))


def test_known_exact_anchor_counts_and_independent_short_walks():
    counts = _exact_anchored_saw_counts16()
    assert counts == (
        1, 3, 9, 25, 71, 195, 543, 1479, 4067,
        11025, 30073, 81233, 220375, 593611, 1604149,
        4311333, 11616669)
    assert exact_anchored_saw16() == 11616669
    for n in range(2, 10):
        assert len(external_saw_paths(n)) == counts[n - 2]
    assert counts[14] <= counts[16]


def test_square_witness_contacts_independent_oracle():
    for side in range(2, 12):
        walk = square_serpentine(side)
        assert len(walk) == side*side
        assert check_square_witness(side) == (side-1)**2
        assert external_hh_contacts(walk) == (side-1)**2


def test_exact_all_h_short_chain_event_oracle():
    # Small independent example: n=4 has 7 walks with K=0 and 2 with K=1.
    scores = [external_hh_contacts(path) for path in external_saw_paths(4)]
    assert sorted(set(scores)) == [0, 1]
    assert scores.count(0) == 7 and scores.count(1) == 2
    q = Fraction(4)
    exact_event = 2*q / (7 + 2*q)
    assert exact_event == Fraction(8, 15)


def test_square_first_population_threshold_and_infinite_progression():
    cert = find_first_square_tail(q=4, density=Fraction(1, 4),
                                  goal_probability=Fraction(99, 100))
    assert cert is not None
    assert cert.side_min == 124
    assert cert.side_step == 4
    assert cert.event_threshold == 3844
    assert cert.witness_contacts == 15129
    assert cert.c16 == 11616669
    assert cert.probability_lower > Fraction(99998, 100000)
    assert cert.failure_odds_upper < Fraction(1, 50_000)
    assert cert.next_side_ratio < 1
    assert cert.next_ratio_multiplier < 1
    # Sample checks corroborate the induction formula, but the whole-tail
    # theorem follows from the exact geometric progression ratio.
    for i in range(4):
        side = cert.side_min + i*cert.side_step
        other = certify_square_density_tail(side)
        assert other.probability_lower >= cert.probability_lower
        _, _, _, direct = _majorant(side, 4, cert.density,
                                   cert.correction_lower,
                                   cert.correction_upper, cert.c16)
        assert direct == cert.failure_odds_upper * cert.next_side_ratio**i * (
            cert.next_ratio_multiplier**(i*(i-1)//2))
    earlier = certify_square_density_tail(120)
    assert earlier.probability_lower < Fraction(99, 100)


def test_tighter_m128_and_adversarial_weight_uncertainty():
    cert = certify_square_density_tail(128)
    assert cert.side_min**2 == 16_384
    assert cert.event_threshold == 4096
    assert cert.witness_contacts == 16129
    assert cert.probability_lower > Fraction(999_999_999, 1_000_000_000)
    assert cert.failure_odds_upper < Fraction(1, 1_000_000_000)
    changed = certify_square_density_tail(
        128, correction_lower=Fraction(1, 3),
        correction_upper=Fraction(3))
    assert changed.probability_lower < cert.probability_lower
    assert changed.probability_lower > Fraction(99, 100)


def test_tail_is_sufficient_only_and_rejects_bad_gate():
    with pytest.raises(ValueError, match="leading entropy"):
        certify_square_density_tail(124, q=3, density=Fraction(1, 4))
    with pytest.raises(ValueError, match="side_min"):
        certify_square_density_tail(122)
    with pytest.raises(ValueError, match="side_min"):
        certify_square_density_tail(True)
    with pytest.raises(ValueError, match="integer >=2"):
        certify_square_density_tail(124, q=4.0)
    with pytest.raises(ValueError, match="integer >=2"):
        certify_square_density_tail(124, q=True)
    with pytest.raises(ValueError, match="density"):
        certify_square_density_tail(124, density=0.25)
    with pytest.raises(ValueError, match="density"):
        certify_square_density_tail(124, density=Fraction(5, 4))
    with pytest.raises(ValueError, match="correction"):
        certify_square_density_tail(124, correction_lower=Fraction(0))
    with pytest.raises(ValueError, match="correction"):
        certify_square_density_tail(124, correction_lower=Fraction(3),
                                    correction_upper=Fraction(2))
    with pytest.raises(ValueError, match="probability"):
        certify_square_density_tail(124, goal_probability=1.0)
    with pytest.raises(ValueError, match="max_side"):
        find_first_square_tail(max_side=3)
    assert find_first_square_tail(max_side=120) is None


def test_leading_gate_exact_large_integer_check():
    c16 = exact_anchored_saw16()
    # rho=1/4, q=4; 4^12 > c16. The strict inequality makes the
    # ratio of consecutive F(m+s)/F(m) shrink geometrically.
    assert 4**12 > c16
    cert = certify_square_density_tail(124)
    assert cert.next_ratio_multiplier == Fraction(c16**2, 4**24)


def test_second_density_regime_and_no_false_favorable_q():
    # rho=1/3, q=5: 5^32 > c16^3; still a sufficient (not optimal) gate.
    assert 5**32 > exact_anchored_saw16()**3
    c = find_first_square_tail(q=5, density=Fraction(1, 3),
                               goal_probability=Fraction(9, 10), max_side=256)
    assert c is not None
    assert c.density == Fraction(1, 3)
    assert c.probability_lower >= Fraction(9, 10)
    assert c.side_min % c.side_step == 0
    # A negative result of this search has no implication of absence
    # of a physical or model-relative phase transition.
    assert find_first_square_tail(q=2, density=Fraction(1, 3),
                                  max_side=48) is None
