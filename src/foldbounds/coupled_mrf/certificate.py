"""Separator certificates for Z_B / Z_A (CEC_THEORY.md section 2), producer side.

For a separator (I, dB, E) the local weights are

    Phi_X(tau) = sum_{x_I} prod_{local factors f of X} f(x_I, tau)      (X = A, B)

for every boundary state tau.  On T+ = {tau : Phi_A(tau) > 0} the local ratio
is rho(tau) = Phi_B(tau) / Phi_A(tau).  Theorem 2.1: if Z_A > 0 and no tau has
Phi_A(tau) = 0 < Phi_B(tau), then  min rho <= Z_B/Z_A <= max rho  for EVERY
exterior.  If some tau has Phi_A(tau) = 0 < Phi_B(tau) the upper end is +inf
(an exterior concentrated near tau makes the ratio arbitrarily large).

All quantities are exact: Phi_X = scale_X * (integer table), so comparisons of
ratios are integer cross-multiplications.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import product
from math import prod
from typing import Mapping, Sequence

from .io import frac_to_str, model_sha256
from .model import CoupledModel
from .separator import (Separator, graph_distances, radius_separator, require_separator,
                        split_factors)
from .tables import Plan, plan_elimination, scaled_tables, sum_product

CERT_FORMAT = "foldbounds.coupled_mrf.certificate"
CERT_VERSION = 1
DEFAULT_MAX_STATES = 2 ** 20


class TooManyBoundaryStates(ValueError):
    pass


@dataclass(frozen=True)
class Enclosure:
    """Sound enclosure of Z_B/Z_A valid for every exterior with Z_A > 0."""

    lower: Fraction | None          # None only when T+ is empty
    upper: Fraction | None          # None means +infinity (or T+ empty)
    lower_witness: tuple[int, ...] | None
    upper_witness: tuple[int, ...] | None
    unbounded_witness: tuple[int, ...] | None
    n_states: int
    n_support: int
    method: str = "exhaustive"
    tight: bool = True
    seconds: float = 0.0
    extra: dict = field(default_factory=dict)

    @property
    def unbounded(self) -> bool:
        return self.unbounded_witness is not None or (self.upper is None and self.lower is not None)

    @property
    def empty_support(self) -> bool:
        return self.lower is None

    def width(self) -> Fraction | None:
        if self.lower is None or self.upper is None:
            return None
        return self.upper - self.lower

    def contains(self, value: Fraction) -> bool:
        if self.lower is None:
            return False
        if value < self.lower:
            return False
        return self.upper is None or value <= self.upper


@dataclass
class LocalData:
    """Local factor lists of A and B for a verified separator."""

    model: CoupledModel
    sep: Separator
    local_a: list
    local_b: list

    @staticmethod
    def build(model: CoupledModel, sep: Separator) -> "LocalData":
        require_separator(model, sep)
        local_a, _ = split_factors(model.shared + model.only_a, sep)
        local_b, _ = split_factors(model.shared + model.only_b, sep)
        return LocalData(model, sep, local_a, local_b)

    def boundary_dims(self) -> tuple[int, ...]:
        return tuple(self.model.domains_a[v] for v in self.sep.boundary)

    def n_states(self) -> int:
        return prod(self.boundary_dims())


def _priority(model: CoupledModel, sep: Separator) -> dict[str, int]:
    dist = graph_distances(model.adjacency(), model.perturbation_support())
    return {v: dist.get(v, 10 ** 9) for v in sep.interior}


class _RatioScan:
    """Running exact min / max of b/a over integer pairs (a > 0)."""

    def __init__(self):
        self.min_num = self.min_den = None
        self.max_num = self.max_den = None
        self.min_idx = self.max_idx = None
        self.unbounded_idx = None
        self.support = 0

    def scan(self, phi_a: Sequence[int], phi_b: Sequence[int], offset: int):
        mn, md, mx, xd = self.min_num, self.min_den, self.max_num, self.max_den
        mi, xi = self.min_idx, self.max_idx
        support = 0
        for k, (a, b) in enumerate(zip(phi_a, phi_b)):
            if a == 0:
                if b != 0 and self.unbounded_idx is None:
                    self.unbounded_idx = offset + k
                continue
            support += 1
            if mn is None:
                mn, md, mx, xd = b, a, b, a
                mi = xi = offset + k
                continue
            if b * md < mn * a:
                mn, md, mi = b, a, offset + k
            if b * xd > mx * a:
                mx, xd, xi = b, a, offset + k
        self.min_num, self.min_den, self.max_num, self.max_den = mn, md, mx, xd
        self.min_idx, self.max_idx = mi, xi
        self.support += support


def decode_state(index: int, dims: Sequence[int]) -> tuple[int, ...]:
    out = []
    for d in reversed(dims):
        out.append(index % d)
        index //= d
    return tuple(reversed(out))


def encode_state(state: Sequence[int], dims: Sequence[int]) -> int:
    k = 0
    for x, d in zip(state, dims):
        k = k * d + x
    return k


def local_tables(local: LocalData, which: str, clamp: Mapping[str, int] | None = None,
                 plan: Plan | None = None) -> tuple[Fraction, list[int], Plan]:
    """Exact Phi_X over the (unclamped) boundary variables, row-major."""
    model, sep = local.model, local.sep
    factors = local.local_a if which == "A" else local.local_b
    domains = model.domains_a if which == "A" else model.domains_b
    keep = [v for v in sep.boundary if not clamp or v not in clamp]
    return sum_product(factors, domains, keep=keep, clamp=clamp,
                       eliminate_vars=sep.interior, priority=_priority(model, sep), plan=plan)


def exhaustive_enclosure(model: CoupledModel, sep: Separator,
                         max_states: int = DEFAULT_MAX_STATES,
                         chunk_states: int = 2 ** 16) -> Enclosure:
    """Exact [min rho, max rho] by enumerating every boundary state.

    The boundary is processed in chunks: the leading boundary variables are
    clamped, so each chunk's table has at most ``chunk_states`` entries.
    """
    t0 = time.perf_counter()
    local = LocalData.build(model, sep)
    dims = local.boundary_dims()
    n = prod(dims)
    if n > max_states:
        raise TooManyBoundaryStates(f"{n} boundary states exceed the exhaustive limit {max_states}")
    k = 0
    rest = n
    while rest > chunk_states and k < len(dims):
        rest //= dims[k]
        k += 1
    head = sep.boundary[:k]
    scan = _RatioScan()
    plans: dict[str, Plan | None] = {"A": None, "B": None}
    scales = {}
    for c_index, cvals in enumerate(product(*(range(d) for d in dims[:k]))):
        clamp = dict(zip(head, cvals))
        tabs = {}
        for which in "AB":
            scale, vals, plan = local_tables(local, which, clamp or None, plans[which])
            plans[which] = plan
            scales[which] = scale
            tabs[which] = vals
        scan.scan(tabs["A"], tabs["B"], c_index * rest)
    factor = scales["B"] / scales["A"] if scales else Fraction(1)
    seconds = time.perf_counter() - t0
    unb = decode_state(scan.unbounded_idx, dims) if scan.unbounded_idx is not None else None
    extra = {"plan_A": plans["A"], "plan_B": plans["B"], "chunk_vars": k}
    if scan.min_num is None:
        return Enclosure(None, None, None, None, unb, n, 0, "exhaustive", True, seconds, extra)
    lower = factor * Fraction(scan.min_num, scan.min_den)
    upper = factor * Fraction(scan.max_num, scan.max_den)
    return Enclosure(lower, None if unb is not None else upper,
                     decode_state(scan.min_idx, dims), decode_state(scan.max_idx, dims), unb,
                     n, scan.support, "exhaustive", True, seconds, extra)


def local_ratio(model: CoupledModel, sep: Separator, state: Sequence[int]) -> tuple[Fraction, Fraction]:
    """Exact (Phi_A(tau), Phi_B(tau)) for one boundary state."""
    local = LocalData.build(model, sep)
    clamp = dict(zip(sep.boundary, state))
    out = []
    for which in "AB":
        scale, vals, _ = local_tables(local, which, clamp)
        out.append(scale * vals[0])
    return out[0], out[1]


def all_local_ratios(model: CoupledModel, sep: Separator,
                     max_states: int = 2 ** 16) -> tuple[list[Fraction], list[Fraction]]:
    """Full exact Phi_A, Phi_B tables (small boundaries only)."""
    local = LocalData.build(model, sep)
    if local.n_states() > max_states:
        raise TooManyBoundaryStates("boundary too large for full tables")
    res = []
    for which in "AB":
        scale, vals, _ = local_tables(local, which)
        res.append([scale * v for v in vals])
    return res[0], res[1]


# ----------------------------------------------------------- certificates
def _state_json(state):
    return None if state is None else list(state)


def certificate_dict(model: CoupledModel, sep: Separator, enc: Enclosure,
                     method_data: dict | None = None,
                     elimination_order: Sequence[str] | None = None) -> dict:
    """Serialise an enclosure as a checkable certificate (JSON-ready dict)."""
    claim = {
        "lower": None if enc.lower is None else frac_to_str(enc.lower),
        "upper": ("inf" if enc.upper is None and enc.lower is not None
                  else (None if enc.upper is None else frac_to_str(enc.upper))),
        "tight": bool(enc.tight),
    }
    cert = {
        "format": CERT_FORMAT,
        "version": CERT_VERSION,
        "model_sha256": model_sha256(model),
        "separator": sep.to_dict(model),
        "method": enc.method,
        "claim": claim,
        "witness": {"lower": _state_json(enc.lower_witness),
                    "upper": _state_json(enc.upper_witness),
                    "unbounded": _state_json(enc.unbounded_witness)},
    }
    if elimination_order is None and enc.extra.get("plan_A") is not None:
        elimination_order = enc.extra["plan_A"].order
    if elimination_order is not None:
        cert["hints"] = {"elimination_order": list(elimination_order)}
    if method_data:
        cert[enc.method] = method_data
    return cert


# --------------------------------------------------------------- hierarchy
@dataclass
class HierarchyLevel:
    radius: int
    separator: Separator | None
    enclosure: Enclosure | None
    note: str = ""


def shell_hierarchy(model: CoupledModel, radii: Sequence[int],
                    max_states: int = DEFAULT_MAX_STATES,
                    enclosure_fn=None) -> list[HierarchyLevel]:
    """Enclosures for the graph-distance shell separators of the given radii.

    Radii whose separator is invalid (a domain-changed variable on the
    boundary at R = 0) or too large are reported with a note.  Theorem 2.4
    predicts nesting, which ``check_nesting`` verifies exactly.
    """
    enclosure_fn = enclosure_fn or (lambda m, s: exhaustive_enclosure(m, s, max_states))
    adj = model.adjacency()
    out = []
    for R in radii:
        sep = radius_separator(model, R, adjacency=adj)
        try:
            enc = enclosure_fn(model, sep)
            out.append(HierarchyLevel(R, sep, enc))
        except TooManyBoundaryStates as exc:
            out.append(HierarchyLevel(R, sep, None, f"skipped: {exc}"))
        except ValueError as exc:
            out.append(HierarchyLevel(R, sep, None, f"invalid: {exc}"))
    return out


def check_nesting(levels: Sequence[HierarchyLevel]) -> list[str]:
    """Return violations of [lo', hi'] subset [lo, hi] for consecutive valid levels."""
    problems = []
    prev = None
    for lev in levels:
        enc = lev.enclosure
        if enc is None or enc.lower is None:
            continue
        if prev is not None:
            p = prev.enclosure
            if enc.lower < p.lower:
                problems.append(f"R={lev.radius}: lower {enc.lower} < previous {p.lower}")
            if p.upper is not None and (enc.upper is None or enc.upper > p.upper):
                problems.append(f"R={lev.radius}: upper {enc.upper} > previous {p.upper}")
        prev = lev
    return problems


# ---------------------------------------------------- contraction (Thm 2.6)
@dataclass
class ContractionDiagnostic:
    inner_radius: int
    outer_radius: int
    dobrushin: float | Fraction | None
    exact: bool
    osc_inner: Fraction
    osc_outer: Fraction
    observed_ratio: float | None
    identity_ok: bool
    bound_ok: bool | None


def annulus_contraction(model: CoupledModel, inner: Separator, outer: Separator,
                        exact_limit: int = 2_000_000, float_limit: int = 2 * 10 ** 8
                        ) -> ContractionDiagnostic:
    """Annulus kernel kappa, its Dobrushin coefficient and Theorem 2.6 check.

    Requires ``outer.interior == inner.interior | inner.boundary``.  The kernel
    K(sigma, tau') is the product of the shared annulus factors (local to the
    outer separator but not to the inner one).  We verify exactly that
    Phi'_X(tau') = sum_sigma Phi_X(sigma) K(sigma, tau') for X = A, B, then
    compute delta(kappa) = max TV(kappa(tau'_1; .), kappa(tau'_2; .)) exactly
    (small cases) or in floating point (diagnostic only), and test
    osc(rho') <= delta * osc(rho).
    """
    if outer.interior != inner.interior | set(inner.boundary):
        raise ValueError("outer separator must absorb exactly the inner boundary")
    loc_in = LocalData.build(model, inner)
    loc_out = LocalData.build(model, outer)
    inner_keys = {id(f) for f in loc_in.local_a}
    annulus = [f for f in loc_out.local_a if id(f) not in inner_keys]
    shared_ids = {id(f) for f in model.shared}
    if any(id(f) not in shared_ids for f in annulus):
        raise ValueError("annulus contains a system-specific factor")
    dims_in = loc_in.boundary_dims()
    dims_out = loc_out.boundary_dims()
    n_in, n_out = prod(dims_in), prod(dims_out)
    if n_in * n_out > 4 * 10 ** 6:
        raise TooManyBoundaryStates("annulus kernel too large")
    phi = {}
    for which, loc in (("A", loc_in), ("B", loc_in)):
        s, v, _ = local_tables(loc, which)
        phi[("in", which)] = [s * x for x in v]
    for which in "AB":
        s, v, _ = local_tables(loc_out, which)
        phi[("out", which)] = [s * x for x in v]
    # kernel K(sigma, tau') via the annulus factors
    domains = model.domains_a
    ann_tables, ann_scale = scaled_tables(annulus, domains)
    scope = list(inner.boundary) + list(outer.boundary)
    from .tables import joint_table
    kern = joint_table(ann_tables, scope, domains)   # row-major: sigma major, tau' minor
    identity_ok = True
    for which in "AB":
        pin, pout = phi[("in", which)], phi[("out", which)]
        for t in range(n_out):
            acc = Fraction(0)
            for s in range(n_in):
                k = kern[s * n_out + t]
                if k:
                    acc += pin[s] * k
            if acc * ann_scale != pout[t]:
                identity_ok = False
    pa_in, pb_in = phi[("in", "A")], phi[("in", "B")]
    pa_out, pb_out = phi[("out", "A")], phi[("out", "B")]
    rho_in = [pb_in[s] / pa_in[s] for s in range(n_in) if pa_in[s] > 0]
    rho_out = [pb_out[t] / pa_out[t] for t in range(n_out) if pa_out[t] > 0]
    osc_in = max(rho_in) - min(rho_in)
    osc_out = max(rho_out) - min(rho_out)
    rows = [t for t in range(n_out) if pa_out[t] > 0]
    pairs_cost = len(rows) * len(rows) * n_in // 2
    delta = None
    exact = False
    if pairs_cost <= exact_limit:
        kap = []
        for t in rows:
            den = pa_out[t] / ann_scale
            kap.append([pa_in[s] * kern[s * n_out + t] / den for s in range(n_in)])
        best = Fraction(0)
        for i in range(len(kap)):
            ki = kap[i]
            for j in range(i + 1, len(kap)):
                kj = kap[j]
                tv = sum((x - y for x, y in zip(ki, kj) if x > y), Fraction(0))
                if tv > best:
                    best = tv
        delta, exact = best, True
    elif pairs_cost <= float_limit:
        kapf = []
        for t in rows:
            den = pa_out[t] / ann_scale
            kapf.append([float(pa_in[s] * kern[s * n_out + t] / den) for s in range(n_in)])
        best = 0.0
        for i in range(len(kapf)):
            ki = kapf[i]
            for j in range(i + 1, len(kapf)):
                tv = sum(x - y for x, y in zip(ki, kapf[j]) if x > y)
                if tv > best:
                    best = tv
        delta = best
    bound_ok = None
    if delta is not None:
        if exact:
            bound_ok = osc_out <= delta * osc_in
        else:
            bound_ok = float(osc_out) <= float(delta) * float(osc_in) * (1 + 1e-9) + 1e-300
    ratio = float(osc_out / osc_in) if osc_in > 0 else None
    return ContractionDiagnostic(inner_radius=-1, outer_radius=-1, dobrushin=delta, exact=exact,
                                 osc_inner=osc_in, osc_outer=osc_out, observed_ratio=ratio,
                                 identity_ok=identity_ok, bound_ok=bound_ok)
