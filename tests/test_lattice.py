from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.lattice import (
    Certificate, CertificateError, Node, ROOT, build_frontier,
    exact_spectrum, refine_frontier, verify_certificate,
)


def exact_values(sequence, q, threshold):
    spectrum = exact_spectrum(sequence)
    z = sum((count * q ** score for score, count in spectrum.items()), Fraction(0))
    event = sum((count * q ** score for score, count in spectrum.items()
                 if score >= threshold), Fraction(0))
    return z, event / z


@pytest.mark.parametrize("sequence", ["HP", "HHHH", "HPHHP", "HHHHHH", "PPHHPPH"])
@pytest.mark.parametrize("q", [Fraction(1), Fraction(3, 2), Fraction(2)])
@pytest.mark.parametrize("threshold", [0, 1, 2])
def test_bounded_probabilities_enclose_exact(sequence, q, threshold):
    z, probability = exact_values(sequence, q, threshold)
    for depth in range(2, len(sequence) + 1):
        result = verify_certificate(build_frontier(sequence, q, threshold, depth))
        assert result.partition_lower <= z <= result.partition_upper
        assert result.probability_lower <= probability <= result.probability_upper
        if depth == len(sequence):
            assert result.exact
            assert result.partition_lower == z
            assert result.probability_lower == probability == result.probability_upper


def test_exhaustive_small_sequences():
    for n in range(2, 8):
        for residues in product("HP", repeat=n):
            sequence = "".join(residues)
            q = Fraction(3, 2)
            z, probability = exact_values(sequence, q, 1)
            for depth in (2, n):
                result = verify_certificate(build_frontier(sequence, q, 1, depth))
                assert result.partition_lower <= z <= result.partition_upper
                assert result.probability_lower <= probability <= result.probability_upper


def test_refinement_preserves_soundness_and_shrinks_upper():
    original = build_frontier("HHHHHHHH", Fraction(2), 1, 2)
    target = original.frontier[0].path
    refined = refine_frontier(original, target)
    before = verify_certificate(original)
    after = verify_certificate(refined)
    exact_z, exact_p = exact_values(original.sequence, original.q, original.event_at_least)
    assert before.partition_lower <= after.partition_lower <= exact_z
    assert exact_z <= after.partition_upper <= before.partition_upper
    assert after.probability_lower <= exact_p <= after.probability_upper


def test_known_contact_probabilities():
    # 4-residue HHHH can make one noncovalent H-H contact around a unit square.
    spectrum = exact_spectrum("HHHH")
    assert spectrum == {0: 7, 1: 2}
    result = verify_certificate(build_frontier("HHHH", Fraction(2), 1, 4))
    assert result.partition_lower == 11
    assert result.probability_lower == Fraction(4, 11)


def test_reject_missing_coverage():
    cert = build_frontier("HHHHH", Fraction(2), 1, 3)
    with pytest.raises(CertificateError):
        verify_certificate(replace(cert, frontier=cert.frontier[:-1]))


def test_reject_duplicate_and_shadowed_nodes():
    cert = build_frontier("HHHHH", Fraction(2), 1, 3)
    with pytest.raises(CertificateError):
        verify_certificate(replace(cert, frontier=cert.frontier + cert.frontier[:1]))
    with pytest.raises(CertificateError):
        verify_certificate(replace(cert, frontier=cert.frontier + (Node(ROOT, "bound"),)))


def test_reject_invalid_paths_and_claims():
    valid = build_frontier("HHHHH", Fraction(2), 1, 3)
    for path in (((0, 0), (1, 0), (0, 0)), ((0, 0), (1, 0), (3, 0)),
                 ((1, 0), (0, 0), (0, 1))):
        with pytest.raises(CertificateError):
            verify_certificate(replace(valid, frontier=(Node(path, "bound"),)))
    with pytest.raises(CertificateError):
        verify_certificate(replace(valid, frontier=(Node(ROOT, "exact"),)))
    with pytest.raises(CertificateError):
        verify_certificate(replace(valid, frontier=(Node(ROOT, "claim"),)))


@pytest.mark.parametrize("bad", [Fraction(0), Fraction(1, 2), 2.0])
def test_reject_invalid_q(bad):
    with pytest.raises(ValueError):
        build_frontier("HHHH", bad, 1, 2)


def test_reject_invalid_sequences_and_depth():
    for sequence in ("", "H", "HXH", "hph"):
        with pytest.raises(ValueError):
            build_frontier(sequence, Fraction(2), 1, 2)
    for depth in (1, 8, 2.0):
        with pytest.raises(ValueError):
            build_frontier("HHHH", Fraction(2), 1, depth)


@pytest.mark.parametrize("n, count", [(2, 1), (3, 3), (4, 9), (5, 25), (6, 71), (7, 195)])
def test_reference_self_avoiding_walk_counts(n, count):
    assert sum(exact_spectrum("P" * n).values()) == count


def test_non_event_event_certificate_vacuity_elimination():
    # All sequences of length 3 have zero noncovalent contacts in this model.
    result = verify_certificate(build_frontier("HHH", Fraction(2), 1, 2))
    assert result.probability_lower == result.probability_upper == 0
