"""Exact rational certificates for globally unique sequence-partition crossovers.

Only partition functions from the same fixed-length model are compared.
They all agree at q=1. Remove that known factor and count sign variations
of the remaining integer coefficient list in t=q-1:

* 0 variations: fixed nonzero sign for all q>1.
* 1 variation: exactly one positive real root with sign reversal.
* >1 variations: ABSTAIN; never guess uniqueness.

Descartes' rule of signs and the intermediate value theorem establish
the sufficiency. Integer coefficients and exact rational brackets are
returned, not floating-point roots.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .universal_landscape import UniversalLandscape, _sequence_to_bits


def _sign(x: int | Fraction) -> int:
    return (x > 0) - (x < 0)


def _eval(poly: tuple[int, ...], x: Fraction) -> Fraction:
    ans = Fraction(0)
    for coeff in reversed(poly):
        ans = ans * x + coeff
    return ans


def _trim(poly: tuple[int, ...]) -> tuple[int, ...]:
    if not poly:
        return (0,)
    last = len(poly)
    while last > 1 and poly[last-1] == 0:
        last -= 1
    return poly[:last]


def _sign_variations(poly: tuple[int, ...]) -> int:
    signs = [_sign(c) for c in poly if c]
    return sum(a != b for a, b in zip(signs, signs[1:]))


@dataclass(frozen=True)
class SequenceOrderCertificate:
    """All-q classification of model partition-function ordering.

    kind: identical / fixed / unique_crossing. In the fixed case,
    sign_below and sign_above are signs of Z_a-Z_b for all q>1.
    In the crossing case, they specify the signs on either side
    of the unique root bracket [q_lower,q_upper].

    Exact integer/rational calculations; not Lean-kernel checked.
    """
    first: str
    second: str
    kind: Literal["identical", "fixed", "unique_crossing"]
    difference_in_q_minus_one: tuple[int, ...]
    sign_below: int
    sign_above: int
    q_lower: Fraction | None
    q_upper: Fraction | None


def certify_sequence_order(
    landscape: UniversalLandscape, first: str, second: str, *,
    width: Fraction = Fraction(1, 1000), max_bisections: int = 256,
) -> SequenceOrderCertificate:
    """Classify a sequence pair over all q>1, or explicitly fail closed.

    Bracket endpoints are tested with exact rational arithmetic.
    Sign-variation count 0 implies fixed sign; count 1 implies exactly one
    positive root by Descartes' rule plus the intermediate value theorem.
    Multiple variations are unsupported, not assumed unique.
    """
    if not isinstance(landscape, UniversalLandscape):
        raise TypeError("expected UniversalLandscape")
    _sequence_to_bits(first, landscape.length)
    _sequence_to_bits(second, landscape.length)
    if type(width) is not Fraction or width <= 0:
        raise ValueError("width must be a positive exact Fraction")
    if type(max_bisections) is not int or not 1 <= max_bisections <= 2048:
        raise ValueError("max_bisections must be an integer within [1,2048]")
    a = landscape.polynomial(first)
    b = landscape.polynomial(second)
    difference = _trim(tuple((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
                             for i in range(max(len(a), len(b)))))
    if difference == (0,):
        return SequenceOrderCertificate(first, second, "identical", difference, 0, 0, None, None)
    if difference[0] != 0:
        raise ValueError("incompatible model: unequal q=1 partition totals")
    first_nonzero = next(i for i, c in enumerate(difference) if c)
    reduced = difference[first_nonzero:]
    variations = _sign_variations(reduced)
    sign_below = _sign(reduced[0])
    sign_above = _sign(reduced[-1])
    if variations == 0:
        return SequenceOrderCertificate(first, second, "fixed", difference,
                                        sign_below, sign_above, None, None)
    if variations > 1:
        raise ValueError("cannot certify a unique positive crossing: multiple sign variations")
    assert sign_below == -sign_above
    low = Fraction(0)
    high = Fraction(1)
    for _ in range(max_bisections):
        value = _eval(reduced, high)
        if value == 0:
            q = high + 1
            return SequenceOrderCertificate(first, second, "unique_crossing", difference,
                                            sign_below, sign_above, q, q)
        if _sign(value) == sign_above:
            break
        high *= 2
    else:
        raise ValueError("failed to bracket unique root within search budget")
    for _ in range(max_bisections):
        if high - low <= width:
            return SequenceOrderCertificate(first, second, "unique_crossing", difference,
                                            sign_below, sign_above, low + 1, high + 1)
        mid = (low + high) / 2
        value = _eval(reduced, mid)
        if value == 0:
            q = mid + 1
            return SequenceOrderCertificate(first, second, "unique_crossing", difference,
                                            sign_below, sign_above, q, q)
        if _sign(value) == sign_below:
            low = mid
        else:
            high = mid
    raise ValueError("failed to isolate unique root to requested precision")
