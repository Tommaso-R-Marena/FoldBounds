"""Exact brute-force marginals (all 2^n states, integer arithmetic) for small systems."""
from __future__ import annotations

import math
from fractions import Fraction

from .model import Model


def exact_marginals(spec: dict, ph) -> list[Fraction]:
    """Return exact Q_i = P(s_i=1)/P(s_i=0) for every site by enumerating all 2^n states.

    Depth-first, one site per level; every state weight is a product of exactly n unary and
    m pair factors, all scaled to integers over one global denominator (cancels in Q).
    """
    mod = Model(spec, Fraction(ph))
    n = mod.n
    nb_prev = [sorted(j for j in mod.nbr[k] if j < k) for k in range(n)]
    vals = []
    for k in range(n):
        vals += [mod.u(k, 0), mod.u(k, 1)]
        for j in nb_prev[k]:
            vals += [mod.phi(k, j, a, b) for a in (0, 1) for b in (0, 1)]
    D = 1
    for v in vals:
        D = D * v.denominator // math.gcd(D, v.denominator)
    def I(v: Fraction) -> int:
        x = v * D
        assert x.denominator == 1
        return int(x)
    U = [[I(mod.u(k, 0)), I(mod.u(k, 1))] for k in range(n)]
    P = []
    for k in range(n):
        P.append([(j, [[I(mod.phi(k, j, a, b)) for b in (0, 1)] for a in (0, 1)]) for j in nb_prev[k]])
    M = [0] * n
    sig = [0] * n

    def rec(k: int, p: int) -> int:
        if k == n:
            return p
        tot = 0
        for s in (0, 1):
            q = p * U[k][s]
            for j, tab in P[k]:
                q *= tab[s][sig[j]]
            sig[k] = s
            sub = rec(k + 1, q)
            if s:
                M[k] += sub
            tot += sub
        return tot

    Z = rec(0, 1)
    return [Fraction(M[i], Z - M[i]) for i in range(n)]


def float_marginals_all_ph(spec: dict, phs):
    """Float64 numpy brute force of theta_i(pH); NOT part of the trusted path."""
    import numpy as np
    n = len(spec["sites"])
    states = ((np.arange(1 << n)[:, None] >> np.arange(n)[None, :]) & 1).astype(np.float64)
    acid = np.array([s["acid"] for s in spec["sites"]], float)
    pka = np.array([float(s["pka"]) for s in spec["sites"]])
    q = states - acid[None, :]
    e = np.zeros(1 << n)  # log10 weight pair part
    for key, k in spec["pairs_milli"].items():
        i, j = map(int, key.split(","))
        e -= (k / 1000.0) * q[:, i] * q[:, j]
    out = []
    for ph in phs:
        lw = e + states @ (pka - ph)
        w = 10.0 ** (lw - lw.max())
        z = w.sum()
        out.append((states * w[:, None]).sum(0) / z)
    return np.array(out)  # (n_ph, n)
