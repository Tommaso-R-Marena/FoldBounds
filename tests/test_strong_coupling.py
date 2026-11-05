"""Adversarial regression for independently checked strong-coupling certificates."""
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.strong_coupling import (
    StrongCouplingCertificate, anchored_walk_count_upper,
    certify_majority_halfline, check_strong_coupling,
    minimum_integer_majority_q, serpentine_witness,
)


def independent_paths(n):
    """Separate exhaustive oracle: no production search or scoring helpers."""
    frontier = [((0, 0), (1, 0))]
    for _ in range(n - 2):
        children = []
        for path in frontier:
            x, y = path[-1]
            for point in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
                if point not in path:
                    children.append(path + (point,))
        frontier = children
    return tuple(frontier)


def independent_score(seq, path):
    result = 0
    for i in range(len(seq)):
        for j in range(i + 2, len(seq)):
            if seq[i] == "H" and seq[j] == "H":
                x, y = path[i]
                xx, yy = path[j]
                if abs(x - xx) + abs(y - yy) == 1:
                    result += 1
    return result


@pytest.mark.parametrize("n", range(2, 9))
def test_block_count_upper_encloses_independent_complete_walk_count(n):
    oracle = independent_paths(n)
    for block in (1, 2, 4, 8):
        cap = anchored_walk_count_upper(n, block_steps=block)
        assert len(oracle) <= cap <= 3 ** (n - 2)
    assert anchored_walk_count_upper(n) == len(oracle)


def test_independent_exhaustive_adversarial_oracle():
    # Iterate all H/P sequences and both sides of an arbitrary per-state
    # correction interval. No independence between state corrections assumed.
    q_values = (Fraction(1), Fraction(3, 2), Fraction(3))
    for n in range(2, 8):
        paths = independent_paths(n)
        sample = (paths[:2] + paths[-2:]) if len(paths) >= 4 else paths[:1]
        sample = tuple(dict.fromkeys(sample))
        for letters in product("HP", repeat=n):
            seq = "".join(letters)
            scores = [independent_score(seq, p) for p in paths]
            for threshold in (0, 1, 2, n):
                for q in q_values:
                    event = sum((q**k for k in scores if k >= threshold), Fraction(0))
                    comp = sum((q**k for k in scores if k < threshold), Fraction(0))
                    cert = StrongCouplingCertificate(
                        seq, q, threshold, sample, Fraction(1, 2), Fraction(2))
                    b = check_strong_coupling(cert)
                    # The exact extremal ratio over all independent corrections
                    # is achieved by choosing lower on one class, upper on other.
                    exact_lo = ((event / 2) / (event / 2 + 2 * comp))
                    exact_hi = ((2 * event) / (2 * event + comp / 2))
                    assert b.probability_lower <= exact_lo <= exact_hi <= b.probability_upper
                    assert b.event_lower <= event / 2 <= 2 * event <= b.event_upper
                    assert b.complement_lower <= comp / 2 <= 2 * comp <= b.complement_upper


def test_complete_oracle_witnesses_give_exact_result_under_no_correction():
    for seq, t, q in (("HHHH", 1, Fraction(2)),
                      ("HPHHPHH", 2, Fraction(5, 2)),
                      ("PPPPPPP", 1, Fraction(1))):
        paths = independent_paths(len(seq))
        spec = [independent_score(seq, p) for p in paths]
        e = sum((q**k for k in spec if k >= t), Fraction(0))
        n = sum((q**k for k in spec if k < t), Fraction(0))
        b = check_strong_coupling(StrongCouplingCertificate(seq, q, t, paths))
        assert b.unseen_walk_count_upper == 0
        assert b.event_lower == b.event_upper == e
        assert b.complement_lower == b.complement_upper == n
        assert b.probability_lower == b.probability_upper == e / (e + n)


def test_large_fifty_residue_fourfold_energy_robustness_and_threshold():
    path = serpentine_witness(10, 5)
    assert len(path) == 50
    assert independent_score("H" * 50, path) == 36
    cert = StrongCouplingCertificate(
        "H" * 50, Fraction(21), 20, (path,), Fraction(1, 2), Fraction(2))
    b = certify_majority_halfline(cert)
    assert b.witness_count == b.event_witness_count == 1
    assert b.all_walk_count_upper == 4067**6
    assert b.guaranteed_event_majority
    assert b.probability_lower > Fraction(3, 5)
    assert minimum_integer_majority_q(cert, max_q=100) == 21
    with pytest.raises(ValueError):
        certify_majority_halfline(replace(cert, q=Fraction(20)))
    for q in (Fraction(21), Fraction(43, 2), Fraction(22),
              Fraction(50), Fraction(100)):
        assert check_strong_coupling(replace(cert, q=q)).guaranteed_event_majority


def test_large_hundred_residue_fourfold_energy_robustness_and_threshold():
    path = serpentine_witness(10, 10)
    assert len(path) == 100
    assert independent_score("H" * 100, path) == 81
    cert = StrongCouplingCertificate(
        "H" * 100, Fraction(12), 40, (path,), Fraction(1, 2), Fraction(2))
    b = certify_majority_halfline(cert)
    assert b.all_walk_count_upper == 4067**12 * 9
    assert b.guaranteed_event_majority
    assert b.probability_lower > Fraction(7, 10)
    assert minimum_integer_majority_q(cert, max_q=100) == 12
    assert not check_strong_coupling(replace(cert, q=Fraction(11))).guaranteed_event_majority


@pytest.mark.parametrize("bad_path", [
    ((0, 0), (1, 0), (0, 0), (0, 1)),
    ((0, 0), (1, 0), (3, 0), (3, 1)),
    ((0, 0), (0, 1), (1, 1), (1, 2)),
    ((0, 0), (1, 0), (2, 0)),
    ((0, 0), (True, 0), (1, 1), (0, 1)),
])
def test_reject_malformed_witness(bad_path):
    cert = StrongCouplingCertificate("HHHH", Fraction(2), 1, (bad_path,))
    with pytest.raises(ValueError):
        check_strong_coupling(cert)


def test_reject_duplicate_and_bad_energy_uncertainty():
    good = serpentine_witness(2, 2)
    base = StrongCouplingCertificate("HHHH", Fraction(2), 1, (good,))
    with pytest.raises(ValueError):
        check_strong_coupling(replace(base, witnesses=(good, good)))
    for invalid_q in (0, Fraction(1, 2), 2.0, True, "2"):
        with pytest.raises(ValueError):
            check_strong_coupling(replace(base, q=invalid_q))
    for l, u in ((Fraction(0), Fraction(1)),
                 (Fraction(2), Fraction(1)),
                 (0.5, Fraction(2)),
                 (Fraction(1), 2.0)):
        with pytest.raises(ValueError):
            check_strong_coupling(replace(base, correction_lower=l, correction_upper=u))
    with pytest.raises(ValueError):
        check_strong_coupling(replace(base, event_at_least=True))
    with pytest.raises(ValueError):
        check_strong_coupling(replace(base, sequence="HPHx"))
    with pytest.raises(ValueError):
        check_strong_coupling(replace(base, witnesses=list(base.witnesses)))
    for invalid_steps in (0, 9, 2.0, True):
        with pytest.raises(ValueError):
            anchored_walk_count_upper(50, block_steps=invalid_steps)


def test_empty_witness_cannot_claim_majority():
    cert = StrongCouplingCertificate("H" * 50, Fraction(1000), 20, ())
    b = check_strong_coupling(cert)
    assert b.probability_lower == 0
    assert not b.guaranteed_event_majority
    assert minimum_integer_majority_q(cert, max_q=10) is None


def test_vacuous_events_are_correctly_classified():
    path = serpentine_witness(2, 2)
    for t, p in ((0, Fraction(1)), (100, Fraction(0))):
        c = StrongCouplingCertificate("HHHH", Fraction(2), t, (path,))
        b = check_strong_coupling(c)
        assert b.probability_lower == b.probability_upper == p
