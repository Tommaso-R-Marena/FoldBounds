"""Regression tests for exact multi-crossover and continuum dominance certificates.

The full sequence data come from a separately audited lattice contact atlas.
No floating-point roots or external CAS are required by the repository tests.
"""
from fractions import Fraction
from random import Random

import pytest

from foldbounds.contact_map_atlas import build_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.exact_phase import (
    _eval, _product, certify_polynomial_signs,
    certify_sequence_comparison, isolate_positive_roots,
)
from foldbounds.fixed_composition_dominance import (
    certify_fixed_composition_dominance, dominates, sequence_from_mask,
)


@pytest.fixture(scope="module")
def model12():
    return build_landscape(build_atlas(12))


def test_multiple_genuine_probability_crossovers(model12):
    # Exactly four H residues in both sequences, in the SAME n=12 universe.
    a, b = "PPPPPHHHPPHP", "PHPPHHPPPPHP"
    assert a.count("H") == b.count("H") == 4
    assert model12.contact_spectrum(a) == {0: 26063, 1: 3800, 2: 210}
    assert model12.contact_spectrum(b) == {0: 25639, 1: 4322, 2: 108, 3: 4}
    phase = certify_sequence_comparison(
        model12, a, b, event_at_least=1, width=Fraction(1, 10_000)
    )
    assert phase.objective == "event_probability"
    assert phase.positive_root_count == 2
    assert phase.cell_signs == (-1, 1, -1)
    assert Fraction(9) < phase.roots[0].lo + 1 < Fraction(10)
    assert Fraction(14) < phase.roots[1].lo + 1 < Fraction(15)
    for interval in phase.roots:
        assert interval.hi - interval.lo <= Fraction(1, 10_000)
    for q, sign in [
        (Fraction(1), -1), (Fraction(2), -1),
        (Fraction(10), 1), (Fraction(12), 1),
        (Fraction(16), -1), (Fraction(20), -1),
    ]:
        # Independent ratio-of-contact-sums, not the Sturm polynomial:
        pa = model12.event_probability(a, q, 1) if hasattr(model12, "event_probability") else None
        def probability(sequence):
            spec = model12.contact_spectrum(sequence)
            z = sum((count * q ** k for k, count in spec.items()), Fraction(0))
            e = sum((count * q ** k for k, count in spec.items()
                     if k >= 1), Fraction(0))
            return e / z
        result = probability(a) - probability(b)
        assert (result > 0) - (result < 0) == sign


def test_known_exact_single_partition_crossing():
    model = build_landscape(build_atlas(10))
    phase = certify_sequence_comparison(
        model, "HHPHPPPPPP", "PHPHPPHPPP"
    )
    assert phase.difference_t == (0, 150, -24)
    assert phase.positive_root_count == 1
    assert phase.roots[0].lo == phase.roots[0].hi == Fraction(25, 4)
    assert phase.cell_signs == (1, -1)


def test_exact_tangencies_and_multiple_roots():
    cases = [
        ((1, -2, 1), 1, (1, 1)),  # (t-1)^2, touch without sign reversal
        ((2, -3, 1), 2, (1, -1, 1)),
        ((2, -7, 7, -2), 3, (1, -1, 1, -1)),
        ((1, 0, 1), 0, (1,)),
        ((0, 150, -24), 1, (1, -1)),
    ]
    for p, n, signs in cases:
        c = certify_polynomial_signs(p)
        assert c.positive_root_count == n
        assert c.cell_signs == signs
        assert len(c.roots) == n
        assert all(r.hi - r.lo <= Fraction(1, 1024) for r in c.roots)
    assert certify_polynomial_signs((0,)).cell_signs == (0,)


def test_random_products_with_known_positive_real_roots():
    rng = Random(1407)
    # These have exact known factorizations, including multiple roots, complex
    # conjugate pairs and negative roots. They are independent of the Sturm
    # algorithm; root counts and sign parities follow from the factors.
    for _ in range(200):
        roots = rng.sample(range(1, 9), rng.randrange(0, 5))
        p = (1,)
        for root in roots:
            for _ in range(rng.randrange(1, 4)):
                p = _product(p, (-root, 1))
        for _ in range(rng.randrange(0, 3)):
            p = _product(p, (rng.randrange(1, 6), 1))
        if rng.choice((True, False)):
            p = _product(p, (1, 0, 1))
        certificate = certify_polynomial_signs(p, width=Fraction(1, 4096))
        assert certificate.positive_root_count == len(roots)
        for v in roots:
            assert any(r.lo <= v <= r.hi for r in certificate.roots)
        for r in certificate.roots:
            assert r.hi - r.lo <= Fraction(1, 4096)


@pytest.mark.parametrize("k,total,winners,survivors", [
    (3, 220, (), 4),
    (4, 495, (), 12),
    (5, 792, ("HPPHPHPPHPPH", "HPPHPPHPHPPH"), 2),
    (6, 924, ("HPPHPHHPHPPH",), 1),
])
def test_all_q_fixed_composition_winners(model12, k, total, winners, survivors):
    cert = certify_fixed_composition_dominance(model12, k)
    assert cert.sequences_checked == total
    assert cert.universally_dominant == winners
    assert len(cert.coefficientwise_nondominated) == survivors
    assert cert.dominated_count == total - survivors
    for winner in winners:
        a = model12.polynomial(winner)
        for mask in range(1 << 12):
            if mask.bit_count() == k:
                b = model12.polynomial(sequence_from_mask(mask, 12))
                assert dominates(a, b)


def test_validation_and_failure_to_claim_unsupported_results(model12):
    for bad in (Fraction(0), 0, 0.5, True):
        with pytest.raises(ValueError):
            isolate_positive_roots((1, -1), width=bad)
    with pytest.raises(ValueError):
        isolate_positive_roots((0, 1))
    with pytest.raises(ValueError):
        isolate_positive_roots((1, -3, 2), max_steps=1)
    with pytest.raises(ValueError):
        certify_sequence_comparison(model12, "P" * 12, "H" * 12, event_at_least=1.5)
    for k in (-1, 13, 4.0, True):
        with pytest.raises(ValueError):
            certify_fixed_composition_dominance(model12, k)
    for invalid in [((1, -2), (1, 0)), ((1.0,), (1,)), ((1,), (1, 2.0))]:
        with pytest.raises(ValueError):
            dominates(*invalid)
    assert not dominates((1, 2), (1, 3))  # undecided != reverse crossing
