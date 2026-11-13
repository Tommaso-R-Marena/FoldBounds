"""Independent exact checker for separator certificates (standard library, Fractions only).

It re-derives the rational weights from the integer/decimal-string spec (it does not import
model.py or the proposer) and verifies, for a claimed bound lambda on Q_i = Z_B/Z_A:
  * the boundary order is exactly the set of non-interior sites coupled to the interior;
  * the decision tree is a valid partition of {0,1}^boundary (every path fixes distinct
    variables, depth <= len(order));
  * for every leaf cube the T1 test (CEC_THEORY Thm 1.1 / 1.2) holds EXACTLY:
      lower: sum_xi (rmin_xi - lam) * (L_xi if rmin_xi >= lam else U_xi) >= 0
      upper: sum_xi (rmax_xi - lam) * (U_xi if rmax_xi >= lam else L_xi) <= 0
"""
from __future__ import annotations

import hashlib
import itertools
import json
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction


class CertificateError(Exception):
    pass


def _p10(e: int) -> Fraction:
    with localcontext() as c:
        c.prec = 60
        x = Decimal(10) ** (Decimal(e) / Decimal(1000))
        c.prec = 15
        c.rounding = ROUND_HALF_EVEN
        x = +x
    return Fraction(x)


def spec_hash(spec: dict) -> str:
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def check(spec: dict, cert: dict) -> bool:
    """Return True iff the certificate is valid; raise CertificateError otherwise."""
    if cert["spec_hash"] != spec_hash(spec):
        raise CertificateError("spec hash mismatch")
    sites = spec["sites"]
    n = len(sites)
    ph = Fraction(cert["ph"])
    if (ph * 1000).denominator != 1:
        raise CertificateError("pH not a multiple of 1/1000")
    lam = Fraction(cert["lambda"])
    side = cert["side"]
    if side not in ("lower", "upper"):
        raise CertificateError("bad side")
    i0 = cert["site"]
    S = list(cert["interior"])
    if S[0] != i0 or len(set(S)) != len(S) or not all(0 <= s < n for s in S):
        raise CertificateError("bad interior")
    acid = [s["acid"] for s in sites]
    nbr = [dict() for _ in range(n)]
    for key, v in spec["pairs_milli"].items():
        a, b = map(int, key.split(","))
        nbr[a][b] = v
        nbr[b][a] = v
    Sset = set(S)
    bd = sorted({t for s in S for t in nbr[s] if t not in Sset})
    order = list(cert["order"])
    if sorted(order) != bd:
        raise CertificateError("boundary order is not the boundary")

    def u(s, sg):
        if not sg:
            return Fraction(1)
        e = (Fraction(sites[s]["pka"]) - ph) * 1000
        return _p10(int(e))

    def phi(a, b, sa, sb):
        k = nbr[a].get(b)
        if not k:
            return Fraction(1)
        return _p10(-k * (sa - acid[a]) * (sb - acid[b]))

    m = len(S) - 1
    regions = list(itertools.product((0, 1), repeat=m))
    # interior weights Psi(x, xi)
    def psi(x, xi):
        cfg = (x,) + xi
        w = Fraction(1)
        for a in range(len(S)):
            w *= u(S[a], cfg[a])
            for b in range(a + 1, len(S)):
                w *= phi(S[a], S[b], cfg[a], cfg[b])
        return w
    d0 = [psi(0, xi) for xi in regions]
    n0 = [psi(1, xi) for xi in regions]  # numerator mass
    # boundary factor tables: h[x][t][tau] per region
    def hval(t, x, xi, tau):
        cfg = (x,) + xi
        w = Fraction(1)
        for a, s in enumerate(S):
            if t in nbr[s]:
                w *= phi(s, t, cfg[a], tau)
        return w
    H0 = [[[hval(t, 0, xi, tau) for tau in (0, 1)] for t in order] for xi in regions]
    H1 = [[[hval(t, 1, xi, tau) for tau in (0, 1)] for t in order] for xi in regions]
    R, B = len(regions), len(order)
    one = Fraction(1)
    # suffix products of min/max of h0 and of the ratio h1/h0 over tau, per region
    sufLo = [[one] * (B + 1) for _ in range(R)]
    sufHi = [[one] * (B + 1) for _ in range(R)]
    sufRlo = [[one] * (B + 1) for _ in range(R)]
    sufRhi = [[one] * (B + 1) for _ in range(R)]
    for r in range(R):
        for k in range(B - 1, -1, -1):
            a, b = H0[r][k]
            sufLo[r][k] = sufLo[r][k + 1] * min(a, b)
            sufHi[r][k] = sufHi[r][k + 1] * max(a, b)
            q = (H1[r][k][0] / a, H1[r][k][1] / b)
            sufRlo[r][k] = sufRlo[r][k + 1] * min(q)
            sufRhi[r][k] = sufRhi[r][k + 1] * max(q)

    def leaf_ok(k, dfix, nfix):
        g = Fraction(0)
        for r in range(R):
            L = dfix[r] * sufLo[r][k]
            U = dfix[r] * sufHi[r][k]
            rat = nfix[r] / dfix[r]
            if side == "lower":
                rmin = rat * sufRlo[r][k]
                g += (rmin - lam) * (L if rmin >= lam else U)
            else:
                rmax = rat * sufRhi[r][k]
                g += (rmax - lam) * (U if rmax >= lam else L)
        return g >= 0 if side == "lower" else g <= 0

    def walk(node, k, dfix, nfix):
        if node == 0:
            if not leaf_ok(k, dfix, nfix):
                raise CertificateError(f"leaf at depth {k} violates the T1 test")
            return
        if not (isinstance(node, list) and len(node) == 2) or k >= B:
            raise CertificateError("malformed tree")
        for tau in (0, 1):
            d2 = [dfix[r] * H0[r][k][tau] for r in range(R)]
            n2 = [nfix[r] * H1[r][k][tau] for r in range(R)]
            walk(node[tau], k + 1, d2, n2)

    walk(cert["tree"], 0, d0, n0)
    return True
