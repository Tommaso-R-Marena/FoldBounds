"""Monte Carlo baseline: coloured Gibbs sampling, many independent chains vectorised."""
from __future__ import annotations

import math

import numpy as np


def _coloring(n, nbr):
    color = [-1] * n
    for v in sorted(range(n), key=lambda x: -len(nbr[x])):
        used = {color[u] for u in nbr[v] if color[u] >= 0}
        c = 0
        while c in used:
            c += 1
        color[v] = c
    return color


def gibbs_theta(spec, ph: float, seed: int, chains: int = 64, burn: int = 1000, sweeps: int = 4000):
    """Return (theta_mean[n] over sweeps and chains, per-chain theta[n, chains])."""
    n = len(spec["sites"])
    acid = np.array([s["acid"] for s in spec["sites"]], float)
    pka = np.array([float(s["pka"]) for s in spec["sites"]])
    nbr = [dict() for _ in range(n)]
    for key, k in spec["pairs_milli"].items():
        a, b = map(int, key.split(","))
        nbr[a][b] = k
        nbr[b][a] = k
    color = _coloring(n, nbr)
    ln10 = math.log(10.0)
    h = (pka - ph) * ln10
    groups = []
    for c in sorted(set(color)):
        idx = np.array([v for v in range(n) if color[v] == c])
        md = max([len(nbr[v]) for v in idx] + [1])
        J = np.full((len(idx), md), n, dtype=int)
        W = np.zeros((len(idx), md))
        for r, v in enumerate(idx):
            for t, (u, k) in enumerate(nbr[v].items()):
                J[r, t] = u
                W[r, t] = k / 1000.0 * ln10
        groups.append((idx, J, W))
    rng = np.random.default_rng(seed)
    S = rng.integers(0, 2, size=(n, chains)).astype(float)
    qpad = np.zeros((n + 1, chains))
    acc = np.zeros((n, chains))
    for sweep in range(burn + sweeps):
        for idx, J, W in groups:
            qpad[:n] = S - acid[:, None]
            field = h[idx][:, None] - (W[:, :, None] * qpad[J]).sum(1)
            p = 1.0 / (1.0 + np.exp(-np.clip(field, -700, 700)))
            S[idx] = (rng.random(p.shape) < p).astype(float)
        if sweep >= burn:
            acc += S
    per_chain = acc / sweeps
    return per_chain.mean(1), per_chain
