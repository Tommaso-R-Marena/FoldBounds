"""Tests of parameterized infinite-size 3D heterogeneous HP tail inequalities."""
from fractions import Fraction
import pytest

from foldbounds.cubic_general_tail import certify_general_cubic_density_tail
from foldbounds.cubic_tiling import TiledCertificate, check_tiled_family


def test_period_five_20_percent_polar_infinite_family():
    result = certify_general_cubic_density_tail(
        20, polar_period=5, density=Fraction(1, 4), q=5
    )
    assert result.side_step == 20
    assert result.residues == 64000
    assert result.polar_residues == 12800
    assert result.hh_contacts_lower == 72000
    assert result.event_threshold == 16000
    assert result.guaranteed_majority
    assert result.failure_odds_upper < Fraction(1, 10**9)
    assert result.probability_lower > Fraction(99, 100)
    assert result.all_walks_upper > result.witness_family_lower


def test_period_four_25_percent_polar_infinite_family():
    result = certify_general_cubic_density_tail(
        24, polar_period=4, density=Fraction(1, 4), q=8
    )
    assert result.side_step == 4
    assert result.residues == 110592
    assert result.polar_residues == 27648
    assert result.hh_contacts_lower == 103680
    assert result.event_threshold == 27648
    assert result.guaranteed_majority
    assert result.failure_odds_upper < Fraction(1, 10**9)


def test_infinite_family_ratio_decreases_between_eligible_sides():
    first = certify_general_cubic_density_tail(
        24, polar_period=4, density=Fraction(1, 4), q=8
    )
    second = certify_general_cubic_density_tail(
        28, polar_period=4, density=Fraction(1, 4), q=8
    )
    assert second.failure_odds_upper < first.failure_odds_upper


def test_independent_finite_transfer_dominates_analytic_lower_family():
    analytic = certify_general_cubic_density_tail(
        20, polar_period=5, density=Fraction(1, 4), q=5
    )
    seq = "".join(
        "P" if i % 5 == 0 else "H" for i in range(analytic.residues)
    )
    finite = check_tiled_family(TiledCertificate(
        (20, 20, 20), seq, Fraction(5), analytic.event_threshold,
        Fraction(1, 2), Fraction(2)
    ))
    assert finite.family_paths >= analytic.witness_family_lower
    assert finite.family_min_contacts >= analytic.hh_contacts_lower
    assert finite.total_walks_upper <= analytic.all_walks_upper
    assert finite.probability_lower >= analytic.probability_lower


def test_unfavorable_parameter_gates_abstain():
    with pytest.raises(ValueError, match="monotonicity"):
        certify_general_cubic_density_tail(
            20, polar_period=5, density=Fraction(1, 4), q=4
        )
    with pytest.raises(ValueError, match="monotonicity"):
        certify_general_cubic_density_tail(
            24, polar_period=4, density=Fraction(1, 4), q=7
        )
    with pytest.raises(ValueError):
        certify_general_cubic_density_tail(
            4, polar_period=2, density=Fraction(1, 4), q=10
        )


def test_adversarial_reweighting_and_fail_closed_input_types():
    default = certify_general_cubic_density_tail(
        20, polar_period=5, density=Fraction(1, 4), q=5
    )
    perturbed = certify_general_cubic_density_tail(
        20, polar_period=5, density=Fraction(1, 4), q=5,
        correction_lower=Fraction(1, 3),
        correction_upper=Fraction(3)
    )
    assert perturbed.failure_odds_upper == default.failure_odds_upper * Fraction(9, 4)
    assert perturbed.probability_lower < default.probability_lower
    assert perturbed.probability_lower > Fraction(99, 100)
    for invalid_period in (0, 1, True, 5.0, Fraction(5)):
        with pytest.raises(ValueError):
            certify_general_cubic_density_tail(
                20, polar_period=invalid_period, density=Fraction(1, 4), q=5
            )
    for invalid_density in (0.25, Fraction(0), Fraction(2)):
        with pytest.raises(ValueError):
            certify_general_cubic_density_tail(
                20, polar_period=5, density=invalid_density, q=5
            )
    for invalid_q in (1, True, Fraction(5), 5.0):
        with pytest.raises(ValueError):
            certify_general_cubic_density_tail(
                20, polar_period=5, density=Fraction(1, 4), q=invalid_q
            )
    for invalid_side in (0, 19, 21, True, 20.0):
        with pytest.raises(ValueError):
            certify_general_cubic_density_tail(
                invalid_side, polar_period=5, density=Fraction(1, 4), q=5
            )
    with pytest.raises(ValueError):
        certify_general_cubic_density_tail(
            20, polar_period=5, density=Fraction(1, 4), q=5,
            correction_lower=Fraction(2), correction_upper=Fraction(1)
        )
