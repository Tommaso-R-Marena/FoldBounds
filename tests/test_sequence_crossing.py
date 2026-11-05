"""Exact continuum sign certificates (no floating-point roots)."""
from fractions import Fraction

import pytest

from foldbounds.contact_map_atlas import build_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.sequence_crossing import certify_sequence_order


def test_exact_fixed_composition_crossover():
    land = build_landscape(build_atlas(10))
    a, b = 'HHPHPPPPPP', 'PHPHPPHPPP'
    assert a.count('H') == b.count('H') == 3
    assert land.polynomial(a) == (4067, 692)
    assert land.polynomial(b) == (4067, 542, 24)
    cert = certify_sequence_order(land, a, b)
    assert cert.kind == 'unique_crossing'
    assert cert.difference_in_q_minus_one == (0, 150, -24)
    assert (cert.sign_below, cert.sign_above) == (1, -1)
    assert cert.q_lower == cert.q_upper == Fraction(29, 4)
    assert land.partition(a, Fraction(2)) > land.partition(b, Fraction(2))
    assert land.partition(a, Fraction(29, 4)) == land.partition(b, Fraction(29, 4))
    assert land.partition(a, Fraction(10)) < land.partition(b, Fraction(10))


def test_irrational_crossing_isolated_with_exact_rationals():
    land = build_landscape(build_atlas(10))
    def seq(mask):
        return ''.join('H' if mask & (1 << i) else 'P' for i in range(10))
    a, b = seq(980), seq(317)
    cert = certify_sequence_order(land, a, b, width=Fraction(1, 10_000))
    assert cert.kind == 'unique_crossing'
    assert cert.sign_below == -1 and cert.sign_above == 1
    assert 0 < cert.q_upper - cert.q_lower <= Fraction(1, 10_000)
    assert land.partition(a, cert.q_lower) < land.partition(b, cert.q_lower)
    assert land.partition(a, cert.q_upper) > land.partition(b, cert.q_upper)
    assert 5 < cert.q_lower < 6


def test_identical_and_fixed_all_q_order():
    land = build_landscape(build_atlas(6))
    cert = certify_sequence_order(land, 'HHHHHH', 'HHHHHH')
    assert cert.kind == 'identical' and cert.q_lower is None
    cert = certify_sequence_order(land, 'HHHHHH', 'PPPPPP')
    assert cert.kind == 'fixed' and (cert.sign_below, cert.sign_above) == (1, 1)
    assert cert.q_lower is None
    assert land.partition('HHHHHH', Fraction(3)) > land.partition('PPPPPP', Fraction(3))


def test_fail_closed_invalid_width_or_resource_limit():
    land = build_landscape(build_atlas(10))
    a, b = 'HHPHPPPPPP', 'PHPHPPHPPP'
    for bad in (0, Fraction(0), Fraction(-1), 1.0):
        with pytest.raises(ValueError):
            certify_sequence_order(land, a, b, width=bad)
    for bad in (0, 1.0, True, 2049):
        with pytest.raises(ValueError):
            certify_sequence_order(land, a, b, max_bisections=bad)
    with pytest.raises(ValueError):
        certify_sequence_order(land, 'H', b)
    with pytest.raises(ValueError):
        certify_sequence_order(land, a, b, max_bisections=1, width=Fraction(1, 100000))
