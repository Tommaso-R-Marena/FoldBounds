"""Regression + independent small-universe oracle for thermal response bounds."""
from collections import Counter
from fractions import Fraction
from itertools import product
from random import Random

import pytest

from foldbounds.lattice import (Certificate, CertificateError, Node,
                                build_frontier, exact_spectrum, refine_frontier)
from foldbounds.thermal_response import (_ratio_extreme,
                                         thermal_response_bounds)


def _independent_counts(sequence):
    """Enumerate anchored SAWs independently of production lattice helpers."""
    paths = [((0, 0), (1, 0))]
    for _ in range(2, len(sequence)):
        next_paths = []
        for walk in paths:
            x, y = walk[-1]
            for point in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
                if point not in walk:
                    next_paths.append(walk + (point,))
        paths = next_paths
    histogram = Counter()
    for walk in paths:
        score = sum(
            sequence[i] == sequence[j] == 'H'
            and abs(walk[i][0] - walk[j][0]) + abs(walk[i][1] - walk[j][1]) == 1
            for i in range(len(sequence))
            for j in range(i + 2, len(sequence))
        )
        histogram[score] += 1
    return dict(histogram)


def _exact_response(histogram, q, threshold):
    z = sum((multiplicity * q ** k for k, multiplicity in histogram.items()), Fraction(0))
    m1 = sum((k * multiplicity * q ** k for k, multiplicity in histogram.items()), Fraction(0))
    m2 = sum((k * k * multiplicity * q ** k for k, multiplicity in histogram.items()), Fraction(0))
    e = sum((multiplicity * q ** k for k, multiplicity in histogram.items() if k >= threshold), Fraction(0))
    em = sum((k * multiplicity * q ** k for k, multiplicity in histogram.items() if k >= threshold), Fraction(0))
    return z, m1 / z, m2 / z, m2 / z - (m1 / z) ** 2, (em * z - e * m1) / z ** 2


@pytest.mark.parametrize('n', range(2, 7))
def test_every_short_hp_sequence_and_extreme_thresholds(n):
    for chars in product('HP', repeat=n):
        sequence = ''.join(chars)
        spectrum = _independent_counts(sequence)
        assert exact_spectrum(sequence) == dict(sorted(spectrum.items()))
        for q in (Fraction(1), Fraction(3, 2), Fraction(3)):
            for threshold in (-1, 1, 3, 99):
                exact = _exact_response(spectrum, q, threshold)
                for depth in sorted({2, n - 1, n}):
                    if depth < 2:
                        continue
                    result = thermal_response_bounds(build_frontier(sequence, q, threshold, depth))
                    z, mean, m2, variance, event_sensitivity = exact
                    assert result.partition_lower <= z <= result.partition_upper
                    assert result.mean_lower <= mean <= result.mean_upper
                    assert result.second_moment_lower <= m2 <= result.second_moment_upper
                    assert result.variance_lower <= variance <= result.variance_upper
                    assert (result.event_susceptibility_lower <= event_sensitivity
                            <= result.event_susceptibility_upper)
                    if depth == n:
                        assert result.exact
                        assert result.partition_lower == result.partition_upper == z
                        assert result.mean_lower == result.mean_upper == mean
                        assert result.second_moment_lower == result.second_moment_upper == m2
                        assert result.variance_lower == result.variance_upper == variance
                        assert (result.event_susceptibility_lower ==
                                result.event_susceptibility_upper == event_sensitivity)


def test_longer_mixed_frontiers_against_independent_oracle():
    rng = Random(20261105)
    for n in (7, 8, 9):
        sequences = {'H' * n, 'P' * n, 'HP' * (n // 2) + ('H' if n % 2 else '')}
        sequences.update(''.join(rng.choice('HP') for _ in range(n)) for _ in range(6))
        for sequence in sorted(sequences):
            spectrum = _independent_counts(sequence)
            for q in (Fraction(1), Fraction(5, 4), Fraction(2)):
                for threshold in (1, 2, 3):
                    exact = _exact_response(spectrum, q, threshold)
                    for depth in (2, 4, n - 1):
                        cert = build_frontier(sequence, q, threshold, depth)
                        result = thermal_response_bounds(cert)
                        _, mean, m2, variance, event_sensitivity = exact
                        assert result.mean_lower <= mean <= result.mean_upper
                        assert result.second_moment_lower <= m2 <= result.second_moment_upper
                        assert result.variance_lower <= variance <= result.variance_upper
                        assert (result.event_susceptibility_lower <= event_sensitivity
                                <= result.event_susceptibility_upper)


def test_refinement_and_mixed_unknown_exact_nodes():
    cert = build_frontier('HHHHHHHH', Fraction(3, 2), 1, 3)
    baseline = thermal_response_bounds(cert)
    assert not baseline.exact
    target = cert.frontier[0].path
    refined = refine_frontier(cert, target)
    after = thermal_response_bounds(refined)
    assert after.partition_upper <= baseline.partition_upper
    assert after.mean_lower >= baseline.mean_lower
    assert after.mean_upper <= baseline.mean_upper
    assert after.variance_lower >= baseline.variance_lower
    assert after.variance_upper <= baseline.variance_upper


def test_terminal_bound_nodes_are_exact_according_to_base_verifier():
    cert = build_frontier('HHHHH', Fraction(2), 1, 5)
    relabeled = Certificate(cert.sequence, cert.q, cert.event_at_least,
                            tuple(Node(node.path, 'bound') for node in cert.frontier))
    expected = thermal_response_bounds(cert)
    assert thermal_response_bounds(relabeled) == expected


def test_rejects_incomplete_duplicate_invalid_and_wrong_q():
    cert = build_frontier('HHHHHH', Fraction(2), 1, 3)
    for invalid in (
        Certificate(cert.sequence, cert.q, cert.event_at_least, cert.frontier[:-1]),
        Certificate(cert.sequence, cert.q, cert.event_at_least, cert.frontier + (cert.frontier[0],)),
        Certificate(cert.sequence, Fraction(1, 2), cert.event_at_least, cert.frontier),
        Certificate(cert.sequence, cert.q, cert.event_at_least,
                    (Node(((0, 0), (1, 0), (0, 0)), 'bound'),) + cert.frontier),
    ):
        with pytest.raises((CertificateError, ValueError)):
            thermal_response_bounds(invalid)


def test_ratio_solver_against_all_subsets_and_fractional_samples():
    rng = Random(1789)
    for _ in range(100):
        n = rng.randrange(0, 6)
        z = Fraction(rng.randrange(1, 8))
        m = z * rng.randrange(0, 9)
        candidates = tuple((rng.randrange(0, 9), Fraction(rng.randrange(1, 9)))
                           for _ in range(n))
        corners = []
        for choice in product((0, 1), repeat=n):
            numerator = m + sum((v * w for (v, w), chosen in zip(candidates, choice) if chosen), Fraction(0))
            denominator = z + sum((w for (_, w), chosen in zip(candidates, choice) if chosen), Fraction(0))
            corners.append(numerator / denominator)
        assert _ratio_extreme(z, m, candidates, minimize=True) == min(corners)
        assert _ratio_extreme(z, m, candidates, minimize=False) == max(corners)
    assert _ratio_extreme(Fraction(0), Fraction(0), ((2, Fraction(9)), (8, Fraction(1))), minimize=True) == 2
    assert _ratio_extreme(Fraction(0), Fraction(0), ((2, Fraction(9)), (8, Fraction(1))), minimize=False) == 8


def test_informative_variance_and_event_response_lower_when_known_classes_exist():
    cert = build_frontier('HHHHHH', Fraction(2), 1, 6)
    result = thermal_response_bounds(cert)
    assert result.variance_lower > 0
    assert result.event_susceptibility_lower > 0
    assert result.exact


def _full_walks(n):
    walks = [((0, 0), (1, 0))]
    for _ in range(2, n):
        following = []
        for path in walks:
            x, y = path[-1]
            for next_site in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
                if next_site not in path:
                    following.append(path + (next_site,))
        walks = following
    return tuple(walks)


def test_validated_witnesses_strengthen_partial_moment_enclosures():
    sequence = 'HHHHHHH'
    q = Fraction(3, 2)
    t = 1
    walks = _full_walks(len(sequence))
    spectrum_scores = []
    for path in walks:
        k = sum(sequence[i] == sequence[j] == 'H'
                and abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) == 1
                for i in range(len(sequence)) for j in range(i + 2, len(sequence)))
        spectrum_scores.append((path, k))
    subset = tuple([path for path, k in spectrum_scores if k == 0][:12]
                   + [path for path, k in spectrum_scores if k > 0][:12])
    cert = build_frontier(sequence, q, t, 2)
    base = thermal_response_bounds(cert)
    improved = thermal_response_bounds(cert, witnesses=subset)
    assert improved.witness_paths == 24
    assert improved.exact_paths == base.exact_paths == 0
    assert improved.partition_lower > base.partition_lower
    assert improved.partition_upper == base.partition_upper
    assert improved.mean_lower >= base.mean_lower
    assert improved.mean_upper <= base.mean_upper
    assert improved.variance_lower > base.variance_lower
    assert improved.variance_upper <= base.variance_upper
    assert improved.event_susceptibility_lower > base.event_susceptibility_lower
    assert improved.event_susceptibility_upper <= base.event_susceptibility_upper
    z, mean, m2, variance, event_response = _exact_response(_independent_counts(sequence), q, t)
    assert improved.mean_lower <= mean <= improved.mean_upper
    assert improved.second_moment_lower <= m2 <= improved.second_moment_upper
    assert improved.variance_lower <= variance <= improved.variance_upper
    assert improved.event_susceptibility_lower <= event_response <= improved.event_susceptibility_upper


def test_witnesses_must_be_distinct_anchored_and_in_unknown_subtree():
    path = _full_walks(6)[0]
    cert = build_frontier('HHHHHH', Fraction(2), 1, 2)
    invalid = (
        (path, path),
        (path[:-1],),
        (path[:-1] + (path[-2],),),
        (path[:-1] + ((200, 200),),),
        (((1, 0), (0, 0)) + path[2:],),
        (list(path),),
    )
    for witnesses in invalid:
        with pytest.raises(ValueError):
            thermal_response_bounds(cert, witnesses=witnesses)
    complete = build_frontier('HHHHHH', Fraction(2), 1, 6)
    with pytest.raises(ValueError):
        thermal_response_bounds(complete, witnesses=(path,))


def test_witnesses_can_exhaust_unexpanded_upper_cap_and_make_result_exact():
    n = 4
    cert = build_frontier('HHHH', Fraction(1), 1, 2)
    paths = _full_walks(n)
    assert len(paths) == 9
    got = thermal_response_bounds(cert, witnesses=paths)
    expected = thermal_response_bounds(build_frontier('HHHH', Fraction(1), 1, n))
    assert got.exact
    assert got.bounded_subtrees == 0
    assert got.witness_paths == 9
    for field in ('mean_lower', 'mean_upper', 'second_moment_lower',
                  'second_moment_upper', 'variance_lower', 'variance_upper',
                  'event_susceptibility_lower', 'event_susceptibility_upper',
                  'partition_lower', 'partition_upper', 'probability_lower',
                  'probability_upper'):
        assert getattr(got, field) == getattr(expected, field)
