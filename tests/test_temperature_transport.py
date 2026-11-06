"""Independent exact enumeration and adversarial temperature-transport tests.

The oracle deliberately avoids FoldBounds' contact/scoring/enumeration helpers.
"""
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.lattice import build_frontier
from foldbounds.temperature_transport import (
    certify_temperature_window, transport_certificate, transport_certificates,
)


def oracle_spectrum(sequence):
    paths = [((0, 0), (1, 0))]
    for _ in range(2, len(sequence)):
        next_paths = []
        for walk in paths:
            x, y = walk[-1]
            for candidate in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
                if candidate not in walk:
                    next_paths.append(walk + (candidate,))
        paths = next_paths
    counts = {}
    for path in paths:
        score = sum(sequence[i] == sequence[j] == "H"
                    and abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) == 1
                    for i in range(len(sequence))
                    for j in range(i + 2, len(sequence)))
        counts[score] = counts.get(score, 0) + 1
    return counts


def exact_masses(spectrum, q, threshold):
    e = sum((count * q**k for k, count in spectrum.items() if k >= threshold), Fraction(0))
    n = sum((count * q**k for k, count in spectrum.items() if k < threshold), Fraction(0))
    return e + n, e, n, e / (e + n)


def assert_encloses(bound, exact):
    z, e, n, p = exact
    assert bound.partition_lower <= z <= bound.partition_upper
    assert bound.event_lower <= e <= bound.event_upper
    assert bound.complement_lower <= n <= bound.complement_upper
    assert bound.probability_lower <= p <= bound.probability_upper


def test_independent_exhaustive_short_sequences_and_reverse_transport():
    source_targets = ((Fraction(1), Fraction(5, 2)),
                      (Fraction(2), Fraction(1)),
                      (Fraction(3, 2), Fraction(3)),
                      (Fraction(3), Fraction(4, 3)))
    for n in range(2, 7):
        for letters in product("HP", repeat=n):
            sequence = "".join(letters)
            spectrum = oracle_spectrum(sequence)
            for threshold in (-1, 0, 1, 2, n):
                for q_from, q_target in source_targets:
                    cert = build_frontier(sequence, q_from, threshold, 2)
                    got = transport_certificate(cert, q_target)
                    assert_encloses(got, exact_masses(spectrum, q_target, threshold))


def test_exact_same_parameter_and_terminal_source_are_preserved():
    for sequence in ("HHHH", "HPHHPH", "HHHHHHH", "PPPPPPP"):
        spectrum = oracle_spectrum(sequence)
        for threshold in (0, 1, 2, 99):
            cert = build_frontier(sequence, Fraction(3, 2), threshold, len(sequence))
            result = transport_certificate(cert, Fraction(3, 2))
            z, e, n, p = exact_masses(spectrum, Fraction(3, 2), threshold)
            assert (result.partition_lower, result.partition_upper) == (z, z)
            assert (result.event_lower, result.event_upper) == (e, e)
            assert (result.complement_lower, result.complement_upper) == (n, n)
            assert (result.probability_lower, result.probability_upper) == (p, p)


def test_intersection_of_independent_temperature_anchors():
    sequence = "HHPHHPHH"
    threshold = 1
    anchors = (build_frontier(sequence, Fraction(3, 2), threshold, 3),
               build_frontier(sequence, Fraction(3), threshold, 4))
    target = Fraction(2)
    a, b = (transport_certificate(c, target, lookahead=1) for c in anchors)
    merged = transport_certificates(anchors, target, lookahead=1)
    assert_encloses(merged, exact_masses(oracle_spectrum(sequence), target, threshold))
    for candidate in (a, b):
        assert merged.partition_lower >= candidate.partition_lower
        assert merged.partition_upper <= candidate.partition_upper
        assert merged.event_lower >= candidate.event_lower
        assert merged.event_upper <= candidate.event_upper
        assert merged.complement_lower >= candidate.complement_lower
        assert merged.complement_upper <= candidate.complement_upper
        assert merged.probability_lower >= candidate.probability_lower
        assert merged.probability_upper <= candidate.probability_upper


def test_continuum_majority_stability_without_sampling():
    # HHHH has 7 noncontact walks and 2 contact walks. At q=2, P=4/11;
    # at q=3, P=6/13. Thus the entire interval is a certified minority.
    initial = build_frontier("HHHH", Fraction(2), 1, 4)
    window = certify_temperature_window((initial,), Fraction(2), Fraction(3))
    assert window.guaranteed_majority == "complement"
    assert window.probability_upper == Fraction(6, 13)
    spectrum = oracle_spectrum("HHHH")
    for q in (Fraction(2), Fraction(9, 4), Fraction(5, 2), Fraction(11, 4), Fraction(3)):
        z, _, _, p = exact_masses(spectrum, q, 1)
        assert window.partition_lower <= z <= window.partition_upper
        assert window.probability_lower <= p <= window.probability_upper
    high = certify_temperature_window(
        (build_frontier("HHHH", Fraction(4), 1, 4),), Fraction(4), Fraction(6))
    assert high.guaranteed_majority == "event"
    assert high.probability_lower == Fraction(8, 15)
    across_crossover = certify_temperature_window((initial,), Fraction(2), Fraction(4))
    assert across_crossover.guaranteed_majority == "undetermined"


def test_fail_closed_on_invalid_target_or_incompatible_anchors():
    cert = build_frontier("HHHHH", Fraction(2), 1, 3)
    for bad in (Fraction(0), Fraction(1, 2), 1.5, 0, True, "2"):
        with pytest.raises(ValueError):
            transport_certificate(cert, bad)
    with pytest.raises(ValueError):
        transport_certificates((), Fraction(2))
    with pytest.raises(ValueError):
        transport_certificates((cert, replace(cert, event_at_least=2)), Fraction(2))
    with pytest.raises(ValueError):
        transport_certificates((cert, replace(cert, sequence="PPPPP")), Fraction(2))
    with pytest.raises(ValueError):
        certify_temperature_window((cert,), Fraction(3), Fraction(2))


def test_malformed_source_certificate_is_never_trusted():
    cert = build_frontier("HHHHH", Fraction(2), 1, 3)
    tampered = replace(cert, frontier=cert.frontier[:-1])
    with pytest.raises(ValueError):
        transport_certificate(tampered, Fraction(5, 2))
