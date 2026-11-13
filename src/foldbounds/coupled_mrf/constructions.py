"""Soft-separator enclosure (Thm 2.7), the sharpness exterior (Thm 2.2) and a
counterexample showing the single-boundary-variable hypothesis of Thm 2.7 is needed.

Everything is exact (integers / Fractions).
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import prod

from .certificate import LocalData, all_local_ratios, decode_state, local_tables
from .model import CoupledModel, Factor
from .separator import Separator, split_factors


class SoftHypothesisError(ValueError):
    pass


# ------------------------------------------------------------- Theorem 2.7
def soft_sensitivities(model: CoupledModel, sep: Separator) -> list[Fraction]:
    """Rational c_v >= max_{o} max_a f_v(a,o) / min_a f_v(a,o) for each boundary variable.

    Requires every exterior factor to touch at most ONE boundary variable (the E_s
    hypothesis); otherwise raises SoftHypothesisError.  For each exterior factor g
    touching boundary variable v we multiply in
    max over assignments of the other scope variables of (max_a g / min_a g).
    A zero minimum with a positive maximum makes c_v infinite -> error.
    """
    bset = set(sep.boundary)
    _, ext_a = split_factors(model.shared + model.only_a, sep)
    _, ext_b = split_factors(model.shared + model.only_b, sep)
    ext = ext_a  # exterior factors are shared (S2), so ext_a == ext_b as factor sets
    c = {v: Fraction(1) for v in sep.boundary}
    for f in ext:
        touch = [v for v in f.scope if v in bset]
        if len(touch) > 1:
            raise SoftHypothesisError(f"exterior factor {f.name or f.scope} touches boundary variables {touch}")
        if not touch:
            continue
        v = touch[0]
        pos = f.scope.index(v)
        dims = [model.domains_a[w] for w in f.scope]
        strides = [1] * len(dims)
        for j in range(len(dims) - 2, -1, -1):
            strides[j] = strides[j + 1] * dims[j + 1]
        worst = Fraction(1)
        others = [j for j in range(len(dims)) if j != pos]
        for idx in product(*(range(dims[j]) for j in others)):
            base = sum(i * strides[j] for i, j in zip(idx, others))
            vals = [f.table[base + a * strides[pos]] for a in range(dims[pos])]
            lo, hi = min(vals), max(vals)
            if lo == 0:
                if hi != 0:
                    raise SoftHypothesisError("exterior factor with a hard zero: infinite sensitivity")
                continue
            worst = max(worst, hi / lo)
        c[v] *= worst
    return [c[v] for v in sep.boundary]


def soft_enclosure_from_tables(phi_a, phi_b, dims, c):
    """Exact [min R(F), max R(F)] over vertex fields F_v(a) in {1, c_v}  (Thm 2.7).

    ``phi_a``/``phi_b``: positive-support Fraction tables (row-major over dims),
    ``c``: list of rationals >= 1.  Cost O(3^n 2^n)-ish for n binary variables
    via progressive contraction.  Returns (lo, hi, lo_field, hi_field).
    """
    n = len(dims)
    den = 1
    for x in list(phi_a) + list(phi_b):
        den = _lcm(den, Fraction(x).denominator)
    A = [int(Fraction(x) * den) for x in phi_a]
    B = [int(Fraction(x) * den) for x in phi_b]
    options = []
    for v in range(n):
        p, q = c[v].numerator, c[v].denominator   # fields scaled by q: {q, p}
        opts = set()
        for bits in product((0, 1), repeat=dims[v]):
            if all(bits):
                bits = (0,) * dims[v]
            opts.add(bits)
        options.append([(bits, [p if b else q for b in bits]) for bits in sorted(opts)])
    best = {"lo": None, "hi": None}

    def rec(k, a, b, chosen):
        if k == n:
            za, zb = a[0], b[0]
            if za == 0:
                return
            if best["lo"] is None or zb * best["lo"][1] < best["lo"][0] * za:
                best["lo"] = (zb, za, tuple(chosen))
            if best["hi"] is None or zb * best["hi"][1] > best["hi"][0] * za:
                best["hi"] = (zb, za, tuple(chosen))
            return
        rest = len(a) // dims[k]
        for bits, F in options[k]:
            na = [sum(F[x] * a[x * rest + j] for x in range(dims[k])) for j in range(rest)]
            nb = [sum(F[x] * b[x * rest + j] for x in range(dims[k])) for j in range(rest)]
            rec(k + 1, na, nb, chosen + [bits])
    rec(0, A, B, [])
    lo, hi = best["lo"], best["hi"]
    return Fraction(lo[0], lo[1]), Fraction(hi[0], hi[1]), lo[2], hi[2]


def _lcm(a, b):
    from math import lcm
    return lcm(a, b)


def soft_enclosure(model: CoupledModel, sep: Separator, c=None, scale: Fraction | None = None):
    """Soft-separator enclosure of Z_B/Z_A for the exterior class E_s.

    ``c``: per-boundary-variable rational bounds (default: derived from the model's
    exterior factors).  ``scale``: if given, c_v := c_v ** scale is NOT used; instead
    pass explicit c.  Returns (lo, hi).
    """
    if c is None:
        c = soft_sensitivities(model, sep)
    phi_a, phi_b = all_local_ratios(model, sep, max_states=2 ** 14)
    lo, hi, _, _ = soft_enclosure_from_tables(phi_a, phi_b, [model.domains_a[v] for v in sep.boundary], c)
    return lo, hi


def soft_rational_c(s: float, max_den: int = 10 ** 6) -> Fraction:
    """A rational >= e^s (rounding up), so the enclosure stays sound."""
    import math
    x = Fraction(math.exp(s)).limit_denominator(max_den)
    while x < Fraction(math.exp(s)) * (1 + Fraction(1, 10 ** 9)):
        x += Fraction(1, max_den)
    return max(x, Fraction(1))


# ----------------------------------------------------------- Theorem 2.2
def sharp_exterior_model(base: CoupledModel, sep: Separator, theta: Fraction):
    """Append the two-boundary-state exterior of Thm 2.2 to ``base``.

    ``base`` must contain only interior + boundary variables (binary boundary).  A new
    variable ``h`` with domain |T| carries hard couplings [h_j = sigma_j] to every
    boundary variable and the unary Psi(h): mass theta/Phi_A(tau-) on the argmin
    state, (1 - theta)/Phi_A(tau+) on the argmax state.  Then Z_A = 1 up to the common
    scale and Z_B/Z_A = theta*min + (1 - theta)*max exactly.
    Returns (model, separator_with_exterior_h, lam).
    """
    from .certificate import exhaustive_enclosure
    enc = exhaustive_enclosure(base, sep)
    dims = [base.domains_a[v] for v in sep.boundary]
    nT = prod(dims)
    phi_a, phi_b = all_local_ratios(base, sep)
    kminus = _encode(enc.lower_witness, dims)
    kplus = _encode(enc.upper_witness, dims)
    psi = [Fraction(0)] * nT
    psi[kminus] += Fraction(theta) / phi_a[kminus]
    psi[kplus] += (1 - Fraction(theta)) / phi_a[kplus]
    h = "h"
    factors = [Factor((h,), tuple(psi), "psi")]
    for j, v in enumerate(sep.boundary):
        tab = []
        for k in range(nT):
            st = decode_state(k, dims)
            for x in range(dims[j]):
                tab.append(Fraction(1) if st[j] == x else Fraction(0))
        factors.append(Factor((h, v), tuple(tab), f"couple_{v}"))
    da = dict(base.domains_a); da[h] = nT
    db = dict(base.domains_b); db[h] = nT
    model = CoupledModel(base.variables + (h,), da, db, base.shared + tuple(factors),
                         base.only_a, base.only_b, {})
    sep2 = Separator(sep.interior, sep.boundary, frozenset({h}))
    lam = Fraction(theta) * enc.lower + (1 - Fraction(theta)) * enc.upper
    return model, sep2, lam


def _encode(state, dims):
    k = 0
    for x, d in zip(state, dims):
        k = k * d + x
    return k


# ---------------------------------------- hypothesis-needed counterexample
def two_boundary_counterexample(c: Fraction = Fraction(2)):
    """Exterior factor touching TWO boundary variables that escapes the vertex-field enclosure.

    Two binary boundary variables with Phi_A = 1 everywhere and Phi_B = [1, 0, 0, 1]
    (B prefers agreement, vanishing on disagreement).  Any exterior factor g(t1, t2)
    changes by at most a factor c when ONE boundary variable flips, so each variable
    has sensitivity c, and the vertex-field enclosure (Thm 2.7) with c_1 = c_2 = c is
    claimed for it.  Take g = [[c, 1], [1, c]]... (single flips change by c, so valid
    log-Lipschitz).  Then Z_B/Z_A = 2c/(2c + 2) while the soft enclosure for
    product exteriors is strictly narrower than that.  Returns a dict of exact values.
    """
    phi_a = [Fraction(1)] * 4
    phi_b = [Fraction(1), Fraction(0), Fraction(0), Fraction(1)]
    lo, hi, _, _ = soft_enclosure_from_tables(phi_a, phi_b, [2, 2], [c, c])
    g = [c, Fraction(1), Fraction(1), c]
    ratio = sum(b * x for b, x in zip(phi_b, g)) / sum(a * x for a, x in zip(phi_a, g))
    g2 = [Fraction(1), c, c, Fraction(1)]
    ratio2 = sum(b * x for b, x in zip(phi_b, g2)) / sum(a * x for a, x in zip(phi_a, g2))
    return {"c": c, "soft_lo": lo, "soft_hi": hi, "ratio_agree_exterior": ratio,
            "ratio_disagree_exterior": ratio2,
            "violates": ratio > hi or ratio2 < lo or ratio < lo or ratio2 > hi}
