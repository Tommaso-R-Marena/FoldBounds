"""Exact coupled linear-fractional enclosures (CEC spec §1: T1 and T1').

Every function here uses exact ``fractions.Fraction`` arithmetic; floats and
bools are rejected. Notation follows ``docs/coupled/CEC_THEORY.md`` §1.

A *region* k has A-mass bounds ``0 <= L_k <= U_k`` and ratio bounds
``0 <= r_k <= R_k``. The admissible set is::

    Adm = {(alpha, beta) : L_k <= alpha_k <= U_k,
                           r_k alpha_k <= beta_k <= R_k alpha_k,
                           sum alpha > 0}

and the target is ``Q = sum beta / sum alpha``. ``Adm`` is nonempty iff
some ``U_k > 0``. An empty region list, or one where every ``U_k`` is zero,
raises :class:`EmptyAdmissibleSet`.

Implemented results (candidate theorems in the spec; proofs are sketched
there and tested in ``tests/test_coupled_ratio.py``):

* ``g_minus`` / ``g_plus``: the linear-time exact tests of Thms 1.1/1.2.
* ``lower_enclosure`` / ``upper_enclosure``: the threshold sort of Thm 1.3.
  They return the optimal value and an optimal vertex configuration.
* ``separate_lower`` / ``separate_upper``: the baseline ``S-`` and ``S+``
  of Thm 1.4 (``None`` stands for ``+inf``).
* ``verify_lp_lower_certificate`` / ``verify_lp_upper_certificate``: the
  T1' dual-certificate verifier for ``P = {alpha >= 0 : G alpha <= h}``.
* ``box_polytope`` / ``box_dual_certificate``: express T1 data as a T1'
  polytope and build its explicit dual certificate.

No novelty is claimed. These are Dinkelbach-style linear-fractional and LP
duality facts specialised to box/ratio data.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Sequence

__all__ = [
    "Region", "Enclosure", "EmptyAdmissibleSet", "InvalidRegion",
    "as_fraction", "validate_regions", "admissible_nonempty",
    "g_minus", "g_plus", "lower_enclosure", "upper_enclosure",
    "separate_lower", "separate_upper", "objective",
    "verify_lp_lower_certificate", "verify_lp_upper_certificate",
    "box_polytope", "box_dual_certificate", "reciprocal_regions",
]


class InvalidRegion(ValueError):
    """Region data violate ``0 <= L <= U`` or ``0 <= r <= R``, or are not exact."""


class EmptyAdmissibleSet(ValueError):
    """No admissible configuration exists: every ``U_k`` is zero, or no region is given."""


def as_fraction(value: object, name: str = "value") -> Fraction:
    """Convert an int or Fraction to Fraction. Reject bool, float and everything else."""
    if isinstance(value, bool):
        raise InvalidRegion(f"{name} must be an exact rational, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    raise InvalidRegion(f"{name} must be int or Fraction, got {type(value).__name__}")


@dataclass(frozen=True)
class Region:
    """T1 data for one region: A-mass in [L, U] and ratio beta/alpha in [r, R]."""
    L: Fraction
    U: Fraction
    r: Fraction
    R: Fraction

    def __post_init__(self) -> None:
        L = as_fraction(self.L, "L")
        U = as_fraction(self.U, "U")
        r = as_fraction(self.r, "r")
        R = as_fraction(self.R, "R")
        if not (0 <= L <= U):
            raise InvalidRegion(f"need 0 <= L <= U, got L={L}, U={U}")
        if not (0 <= r <= R):
            raise InvalidRegion(f"need 0 <= r <= R, got r={r}, R={R}")
        object.__setattr__(self, "L", L)
        object.__setattr__(self, "U", U)
        object.__setattr__(self, "r", r)
        object.__setattr__(self, "R", R)


@dataclass(frozen=True)
class Enclosure:
    """An optimal endpoint and a vertex configuration that attains it.

    ``alpha`` and ``beta`` are listed in the caller's region order.
    ``at_upper`` marks the regions that receive ``U`` (the threshold set).
    """
    value: Fraction
    alpha: tuple[Fraction, ...]
    beta: tuple[Fraction, ...]
    at_upper: tuple[bool, ...]

    @property
    def total_alpha(self) -> Fraction:
        return sum(self.alpha, Fraction(0))

    @property
    def total_beta(self) -> Fraction:
        return sum(self.beta, Fraction(0))


def validate_regions(regions: Iterable[Region]) -> tuple[Region, ...]:
    out = tuple(regions)
    for reg in out:
        if not isinstance(reg, Region):
            raise InvalidRegion("every region must be a Region instance")
    return out


def admissible_nonempty(regions: Sequence[Region]) -> bool:
    """``Adm`` is nonempty iff some region can carry positive A-mass."""
    return any(reg.U > 0 for reg in validate_regions(regions))


def _require_nonempty(regions: Sequence[Region]) -> tuple[Region, ...]:
    regs = validate_regions(regions)
    if not any(reg.U > 0 for reg in regs):
        raise EmptyAdmissibleSet("admissible set is empty (no region has U > 0)")
    return regs


def g_minus(regions: Sequence[Region], lam: Fraction) -> Fraction:
    """g-(lam) = sum_k (r_k - lam) * (L_k if r_k >= lam else U_k)."""
    lam = as_fraction(lam, "lambda")
    total = Fraction(0)
    for reg in validate_regions(regions):
        diff = reg.r - lam
        total += diff * (reg.L if diff >= 0 else reg.U)
    return total


def g_plus(regions: Sequence[Region], lam: Fraction) -> Fraction:
    """g+(lam) = sum_k (R_k - lam) * (U_k if R_k >= lam else L_k)."""
    lam = as_fraction(lam, "lambda")
    total = Fraction(0)
    for reg in validate_regions(regions):
        diff = reg.R - lam
        total += diff * (reg.U if diff >= 0 else reg.L)
    return total


def objective(alpha: Sequence[Fraction], beta: Sequence[Fraction]) -> Fraction:
    """Q = sum beta / sum alpha, for a configuration with positive total A-mass."""
    den = sum(alpha, Fraction(0))
    if den <= 0:
        raise ValueError("configuration has zero total A-mass")
    return sum(beta, Fraction(0)) / den


def _threshold(regs: tuple[Region, ...], key: str, descending: bool) -> Enclosure:
    """Threshold sort shared by the lower (key 'r') and upper (key 'R') enclosures.

    Regions are grouped by their ratio key. Only thresholds between distinct
    key values need to be evaluated: at the optimum lambda*, regions with key
    strictly beyond lambda* take U and those strictly before take L. For a tie
    group whose key equals lambda*, either choice gives the same value.
    Taking U for the whole group never empties the configuration.
    """
    groups: dict[Fraction, list[int]] = {}
    for k, reg in enumerate(regs):
        groups.setdefault(getattr(reg, key), []).append(k)
    ordered = sorted(groups, reverse=descending)
    # Start with every region at L, then move whole groups to U in order.
    num = sum((getattr(reg, key) * reg.L for reg in regs), Fraction(0))
    den = sum((reg.L for reg in regs), Fraction(0))
    best_value: Fraction | None = None
    best_t = -1
    if den > 0:
        best_value, best_t = num / den, 0
    for t, value in enumerate(ordered, start=1):
        for k in groups[value]:
            gap = regs[k].U - regs[k].L
            num += value * gap
            den += gap
        if den > 0:
            cand = num / den
            better = (best_value is None or
                      (cand > best_value if descending else cand < best_value))
            if better:
                best_value, best_t = cand, t
    if best_value is None:  # pragma: no cover - excluded by _require_nonempty
        raise EmptyAdmissibleSet("admissible set is empty")
    upper_keys = set(ordered[:best_t])
    at_upper = tuple(getattr(reg, key) in upper_keys for reg in regs)
    alpha = tuple(reg.U if up else reg.L for reg, up in zip(regs, at_upper))
    beta = tuple(getattr(reg, key) * a for reg, a in zip(regs, alpha))
    enc = Enclosure(best_value, alpha, beta, at_upper)
    if objective(alpha, beta) != best_value:  # pragma: no cover - internal invariant
        raise AssertionError("threshold configuration does not attain its value")
    return enc


def lower_enclosure(regions: Sequence[Region]) -> Enclosure:
    """lambda- = min over Adm of Q, by the threshold sort of Theorem 1.3.

    Regions with the smallest r receive U; all others receive L; and
    beta = r * alpha. The returned configuration is admissible and attains
    the minimum. It satisfies g_minus(value) == 0.
    """
    return _threshold(_require_nonempty(regions), "r", descending=False)


def upper_enclosure(regions: Sequence[Region]) -> Enclosure:
    """lambda+ = max over Adm of Q. Regions with the largest R receive U, and beta = R * alpha."""
    return _threshold(_require_nonempty(regions), "R", descending=True)


def separate_lower(regions: Sequence[Region]) -> Fraction:
    """S- = sum r_k L_k / sum U_k: the lower bound from separate Z_B, Z_A bounds."""
    regs = _require_nonempty(regions)
    return (sum((reg.r * reg.L for reg in regs), Fraction(0)) /
            sum((reg.U for reg in regs), Fraction(0)))


def separate_upper(regions: Sequence[Region]) -> Fraction | None:
    """S+ = sum R_k U_k / sum L_k, or ``None`` (+inf) when every L_k = 0."""
    regs = _require_nonempty(regions)
    den = sum((reg.L for reg in regs), Fraction(0))
    if den == 0:
        return None
    return sum((reg.R * reg.U for reg in regs), Fraction(0)) / den


def reciprocal_regions(regions: Sequence[Region]) -> tuple[Region, ...]:
    """Re-express T1 data with B as the reference ensemble, for bounding Z_A/Z_B.

    Requires ``r_k > 0`` for every region with ``U_k > 0``. The B-mass of
    region k then lies in ``[r_k L_k, R_k U_k]``, and alpha/beta lies in
    ``[1/R_k, 1/r_k]``. Regions with ``U_k = 0`` are inert and are dropped.
    This uses no mass bound beyond the original data. Callers that know a
    tighter B cap should build the reciprocal regions directly.
    """
    out = []
    for reg in validate_regions(regions):
        if reg.U == 0:
            continue
        if reg.r <= 0:
            raise InvalidRegion("reciprocal data need r > 0 on every massive region")
        out.append(Region(reg.r * reg.L, reg.R * reg.U, 1 / reg.R, 1 / reg.r))
    return tuple(out)


# ---------------------------------------------------------------- T1' (LP duals)

Matrix = Sequence[Sequence[Fraction]]


def _check_lp_shapes(G: Matrix, h: Sequence[Fraction], coeffs: Sequence[Fraction],
                     y: Sequence[Fraction]) -> tuple[list[list[Fraction]], list[Fraction],
                                                    list[Fraction], list[Fraction]]:
    K = len(coeffs)
    Gm = [[as_fraction(v, "G entry") for v in row] for row in G]
    hv = [as_fraction(v, "h entry") for v in h]
    cv = [as_fraction(v, "ratio") for v in coeffs]
    yv = [as_fraction(v, "y entry") for v in y]
    if len(Gm) != len(hv) or len(yv) != len(hv) or any(len(row) != K for row in Gm):
        raise ValueError("inconsistent T1' certificate dimensions")
    return Gm, hv, cv, yv


def verify_lp_lower_certificate(G: Matrix, h: Sequence[Fraction], r: Sequence[Fraction],
                                lam: Fraction, y: Sequence[Fraction]) -> bool:
    """Check a T1' dual certificate that lam <= Q on {alpha >= 0 : G alpha <= h}.

    Accept iff ``y >= 0``, ``(r - lam*1) + G^T y >= 0`` componentwise, and
    ``h^T y <= 0``. Soundness: for alpha in P,
    ``sum (r_k - lam) alpha_k >= -(G^T y)^T alpha >= -y^T h >= 0``, and
    ``beta_k >= r_k alpha_k``. A True result therefore proves ``Q >= lam``
    at every admissible point with ``sum alpha > 0``. False means only that
    this y does not certify lam.
    """
    lam = as_fraction(lam, "lambda")
    Gm, hv, rv, yv = _check_lp_shapes(G, h, r, y)
    if any(v < 0 for v in yv):
        return False
    for k in range(len(rv)):
        if rv[k] - lam + sum((Gm[i][k] * yv[i] for i in range(len(yv))), Fraction(0)) < 0:
            return False
    return sum((hv[i] * yv[i] for i in range(len(yv))), Fraction(0)) <= 0


def verify_lp_upper_certificate(G: Matrix, h: Sequence[Fraction], R: Sequence[Fraction],
                                lam: Fraction, y: Sequence[Fraction]) -> bool:
    """Check a T1' dual certificate that Q <= lam on {alpha >= 0 : G alpha <= h}.

    Accept iff ``y >= 0``, ``(lam*1 - R) + G^T y >= 0`` and ``h^T y <= 0``.
    Then ``sum (lam - R_k) alpha_k >= 0``, and ``beta_k <= R_k alpha_k``.
    """
    lam = as_fraction(lam, "lambda")
    Gm, hv, Rv, yv = _check_lp_shapes(G, h, R, y)
    if any(v < 0 for v in yv):
        return False
    for k in range(len(Rv)):
        if lam - Rv[k] + sum((Gm[i][k] * yv[i] for i in range(len(yv))), Fraction(0)) < 0:
            return False
    return sum((hv[i] * yv[i] for i in range(len(yv))), Fraction(0)) <= 0


def box_polytope(regions: Sequence[Region]) -> tuple[list[list[Fraction]], list[Fraction]]:
    """Box constraints as ``G alpha <= h``: the rows are alpha_k <= U_k, then -alpha_k <= -L_k."""
    regs = validate_regions(regions)
    K = len(regs)
    G: list[list[Fraction]] = []
    h: list[Fraction] = []
    for k, reg in enumerate(regs):
        row = [Fraction(0)] * K
        row[k] = Fraction(1)
        G.append(row)
        h.append(reg.U)
    for k, reg in enumerate(regs):
        row = [Fraction(0)] * K
        row[k] = Fraction(-1)
        G.append(row)
        h.append(-reg.L)
    return G, h


def box_dual_certificate(regions: Sequence[Region], lam: Fraction,
                         side: str = "lower") -> list[Fraction]:
    """Explicit dual vector for the box polytope of :func:`box_polytope`.

    For the lower side, ``y^U_k = max(0, lam - r_k)`` and
    ``y^L_k = max(0, r_k - lam)``. The reduced costs are then all zero and
    ``h^T y = -g_minus(lam)``. So the certificate verifies iff
    ``g_minus(lam) >= 0``, which is exactly Theorem 1.1.
    The upper side is symmetric, with ``h^T y = g_plus(lam)``.
    """
    lam = as_fraction(lam, "lambda")
    regs = validate_regions(regions)
    if side == "lower":
        yu = [max(Fraction(0), lam - reg.r) for reg in regs]
        yl = [max(Fraction(0), reg.r - lam) for reg in regs]
    elif side == "upper":
        yu = [max(Fraction(0), reg.R - lam) for reg in regs]
        yl = [max(Fraction(0), lam - reg.R) for reg in regs]
    else:
        raise ValueError("side must be 'lower' or 'upper'")
    return yu + yl
