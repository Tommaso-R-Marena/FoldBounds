"""Exact continuum phase atlas for pairs of finite HP-model observables.

Uses only stdlib Fraction arithmetic. Every open interval returned is proven
root-free by Sturm's theorem; all positive real roots are isolated, including
tangencies and cases with many coefficient sign variations.

This Python program is a research checker, not a Lean-refined executable.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import comb, gcd
from typing import Literal

Polynomial = tuple[int, ...]  # lowest exponent first


def _trim(a):
    a = tuple(a)
    i = len(a)
    while i > 1 and a[i - 1] == 0:
        i -= 1
    return a[:i] or (0,)


def _poly_deriv(a):
    return _trim(tuple(i * a[i] for i in range(1, len(a))))


def _eval(p, x):
    result = Fraction(0)
    for a in reversed(p):
        result = result * x + a
    return result


def _divmod(a, b):
    a = list(_trim(Fraction(x) for x in a))
    b = _trim(Fraction(x) for x in b)
    if b == (0,):
        raise ZeroDivisionError("zero polynomial divisor")
    if len(a) < len(b):
        return (Fraction(0),), _trim(a)
    quotient = [Fraction(0)] * (len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        i = len(a) - len(b)
        factor = a[-1] / b[-1]
        quotient[i] = factor
        for j in range(len(b)):
            a[i + j] -= factor * b[j]
        a = list(_trim(a))
    return _trim(quotient), _trim(a)


def _monic_gcd(a, b):
    a = _trim(Fraction(x) for x in a)
    b = _trim(Fraction(x) for x in b)
    while b != (0,):
        _, rem = _divmod(a, b)
        a, b = b, rem
    return _trim(x / a[-1] for x in a) if a != (0,) else (Fraction(0),)


def _squarefree(p):
    p = _trim(p)
    if len(p) < 2:
        return p
    divisor = _monic_gcd(p, _poly_deriv(p))
    quotient, remainder = _divmod(p, divisor)
    if remainder != (0,):
        raise AssertionError("polynomial gcd division failed")
    return quotient


def _primitive(p):
    """Positive scale normalization; sign of p(x) unchanged."""
    p = _trim(p)
    if p == (0,):
        return (0,)
    import math
    l = 1
    for v in p:
        l = math.lcm(l, Fraction(v).denominator)
    ints = [int(v * l) for v in p]
    g = 0
    for i in ints:
        g = gcd(g, abs(i))
    return _trim(x // g for x in ints)


def sturm_sequence(p):
    """Signed Euclidean remainder sequence of squarefree p over rationals."""
    p = _trim(p)
    if p == (0,):
        raise ValueError("zero polynomial has infinitely many roots")
    if len(p) == 1:
        return (p,)
    sq = _primitive(_squarefree(p))
    first = _trim(sq)
    second = _primitive(_poly_deriv(first))
    seq = [first, second]
    while len(seq[-1]) > 1:
        _, rem = _divmod(seq[-2], seq[-1])
        if rem == (0,):
            raise ArithmeticError("squarefree polynomial has nontrivial gcd")
        seq.append(_primitive(tuple(-r for r in rem)))
    return tuple(seq)


def variations(sturm, x):
    signs = []
    for p in sturm:
        v = _eval(p, x)
        if v:
            signs.append(1 if v > 0 else -1)
    return sum(a != b for a, b in zip(signs, signs[1:]))


def root_count_open(sturm, a, b):
    if not (a < b):
        raise ValueError("interval must have positive width")
    if _eval(sturm[0], a) == 0 or _eval(sturm[0], b) == 0:
        raise ValueError("root endpoints require explicit handling")
    return variations(sturm, a) - variations(sturm, b)


def positive_bound(p):
    """Strict rational Cauchy upper bound for absolute root magnitudes."""
    p = _trim(p)
    if len(p) < 2:
        return Fraction(1)
    leading = abs(Fraction(p[-1]))
    bound = max((abs(Fraction(c)) / leading for c in p[:-1]), default=Fraction(0))
    return Fraction(bound.numerator // bound.denominator + 2)


@dataclass(frozen=True)
class RootInterval:
    lo: Fraction  # t=q-1
    hi: Fraction
    sign_left: int
    sign_right: int

    @property
    def tangent(self):
        return self.sign_left == self.sign_right


@dataclass(frozen=True)
class PhaseCertificate:
    difference_t: Polynomial  # polynomial in t=q-1
    roots: tuple[RootInterval, ...]
    cell_signs: tuple[int, ...]
    positive_root_count: int
    width: Fraction
    objective: Literal["partition", "event_probability"]
    threshold: int | None


def isolate_positive_roots(p, *, width=Fraction(1, 1024), max_steps=100_000):
    """Return disjoint rational intervals containing *every* positive real root.

    For an exactly hit rational root, return (t,t). Other intervals have
    non-root rational endpoints and contain precisely one distinct root.
    Raises if computation budget is exceeded, never returns partial work.
    """
    if type(width) is not Fraction or width <= 0:
        raise ValueError("width must be a positive exact Fraction")
    if type(max_steps) is not int or max_steps < 1:
        raise ValueError("max_steps must be positive integer")
    p = _trim(p)
    if any(type(c) is not int for c in p) or p == (0,):
        raise ValueError("nonzero integer polynomial required")
    if len(p) == 1:
        return ()
    if p[0] == 0:
        raise ValueError("factor zero roots before isolating positive roots")
    seq = sturm_sequence(p)
    high = positive_bound(p)
    assert _eval(seq[0], 0) and _eval(seq[0], high)
    expected = root_count_open(seq, Fraction(0), high)
    pending = [(Fraction(0), high, expected)]
    found = []
    steps = 0
    while pending:
        steps += 1
        if steps > max_steps:
            raise ValueError("root-isolation subdivision budget exceeded")
        a, b, count = pending.pop()
        if count == 0:
            continue
        if count == 1:
            middle = (a + b) / 2
            if _eval(p, middle) == 0:
                found.append((middle, middle))
                continue
            if b - a <= width:
                found.append((a, b))
                continue
        cut = None
        for denom in range(2, 100):
            candidate = a + (b - a) / denom
            if _eval(p, candidate) != 0 and _eval(seq[0], candidate) != 0:
                cut = candidate
                break
        if cut is None:
            raise ValueError("no suitable rational split found")
        left = root_count_open(seq, a, cut)
        right = root_count_open(seq, cut, b)
        if left + right != count:
            raise ArithmeticError("Sturm count conservation failure")
        if right:
            pending.append((cut, b, right))
        if left:
            pending.append((a, cut, left))
    found.sort()
    if len(found) != expected or any(
        found[i][1] >= found[i + 1][0] for i in range(len(found) - 1)
    ):
        raise ArithmeticError("root isolation count/order mismatch")
    return tuple(found)


def certify_polynomial_signs(p, *, width=Fraction(1, 1024), max_steps=100_000,
                             objective="partition", threshold=None):
    """Certify full positive-t sign topology, including repeated roots."""
    p = _trim(p)
    if p == (0,):
        return PhaseCertificate(p, (), (0,), 0, width, objective, threshold)
    if p[0] == 0:
        i = next(i for i, v in enumerate(p) if v)
        reduced = p[i:]
    else:
        reduced = p
    intervals = isolate_positive_roots(reduced, width=width, max_steps=max_steps)
    bound = positive_bound(reduced)
    samples = []
    start = Fraction(0)
    for a, b in intervals:
        if not start <= a:
            raise ArithmeticError("root ordering failed")
        samples.append((start + a) / 2 if a > start else Fraction(0))
        start = b
    samples.append(max(bound, start + 1))
    signs = tuple(
        (1 if v > 0 else -1 if v < 0 else 0)
        for v in (_eval(reduced, x) for x in samples)
    )
    if any(s == 0 for s in signs):
        raise ArithmeticError("cell sign sample fell on a polynomial root")
    roots = tuple(
        RootInterval(a, b, signs[i], signs[i + 1])
        for i, (a, b) in enumerate(intervals)
    )
    return PhaseCertificate(p, roots, signs, len(roots), width, objective, threshold)


def _product(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            out[i + j] += u * v
    return _trim(out)


def _diff(a, b):
    return _trim(
        (a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    )


def q_coefficients_to_t(spectrum, threshold=None):
    """Convert sum_{k>=threshold} N_k q^k to t=q-1 basis."""
    if not isinstance(spectrum, dict) or not spectrum:
        raise ValueError("invalid spectrum")
    if any(
        type(k) is not int or type(v) is not int or k < 0 or v < 0
        for k, v in spectrum.items()
    ):
        raise ValueError("spectrum must be nonnegative integer counts")
    max_k = max(spectrum)
    return _trim(tuple(
        sum(n * comb(k, j) for k, n in spectrum.items()
            if k >= j and (threshold is None or k >= threshold))
        for j in range(max_k + 1)
    ))


def certify_sequence_comparison(landscape, first, second, *,
                                event_at_least=None, width=Fraction(1, 1024),
                                max_steps=100_000):
    """Prove every ranking region for Z or threshold-event probability.

    For event probabilities, denominators are strictly positive for q>=1;
    the sign of P_A-P_B is the sign of E_A Z_B-E_B Z_A. All components
    are exact integer polynomials in t=q-1.
    """
    a = landscape.polynomial(first)
    b = landscape.polynomial(second)
    if event_at_least is None:
        p = _diff(a, b)
        objective = "partition"
    else:
        if type(event_at_least) is not int:
            raise ValueError("event_at_least must be integer")
        ea = q_coefficients_to_t(landscape.contact_spectrum(first), event_at_least)
        eb = q_coefficients_to_t(landscape.contact_spectrum(second), event_at_least)
        p = _diff(_product(ea, b), _product(eb, a))
        objective = "event_probability"
    return certify_polynomial_signs(
        p, width=width, max_steps=max_steps,
        objective=objective, threshold=event_at_least
    )
