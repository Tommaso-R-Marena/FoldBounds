"""Float (numpy) proposer of separator certificates: branch-and-bound over boundary cubes.

Nothing here is trusted: every claim is re-verified by checker.check in exact rationals.
"""
from __future__ import annotations

import heapq
import itertools
import math
from fractions import Fraction

import numpy as np

from .checker import spec_hash

LN10 = math.log(10.0)
MARGIN = 1e-9


class Prepared:
    """Log-domain tables for site i, interior S (S[0]==i), pH (Fraction)."""

    def __init__(self, spec, nbr, i, S, ph):
        self.spec, self.i, self.S, self.ph = spec, i, list(S), Fraction(ph)
        sites = spec["sites"]
        acid = [s["acid"] for s in sites]
        Sset = set(S)
        bd = sorted({t for s in S for t in nbr[s] if t not in Sset})
        m = len(S) - 1
        R = 1 << m
        X = ((np.arange(R)[:, None] >> np.arange(m)[None, :]) & 1) if m else np.zeros((1, 0), int)
        cfg = {x: np.concatenate([np.full((R, 1), x), X], axis=1) for x in (0, 1)}
        php = float(self.ph)

        def lu(s):
            return (float(Fraction(sites[s]["pka"])) - php) * LN10

        def lphi(a, b, sa, sb):  # arrays sa, sb
            k = nbr[a].get(b, 0)
            if not k:
                return 0.0
            return -(k / 1000.0) * LN10 * (sa - acid[a]) * (sb - acid[b])

        lpsi = {}
        for x in (0, 1):
            c = cfg[x]
            tot = np.zeros(R)
            for a, s in enumerate(S):
                tot = tot + c[:, a] * lu(s)
                for b in range(a + 1, len(S)):
                    tot = tot + lphi(s, S[b], c[:, a], c[:, b])
            lpsi[x] = tot
        B = len(bd)
        H = {x: np.zeros((R, B, 2)) for x in (0, 1)}
        for kt, t in enumerate(bd):
            for tau in (0, 1):
                for x in (0, 1):
                    tot = np.zeros(R)
                    for a, s in enumerate(S):
                        if t in nbr[s]:
                            tot = tot + lphi(s, t, cfg[x][:, a], tau)
                    H[x][:, kt, tau] = tot
        strength = np.max(np.abs(np.concatenate([H[0][:, :, 1] - H[0][:, :, 0],
                                                  H[1][:, :, 1] - H[1][:, :, 0]], axis=0)), axis=0) if B else np.zeros(0)
        perm = list(np.argsort(-strength, kind="stable"))
        self.order = [bd[p] for p in perm]
        self.H0 = H[0][:, perm, :]
        self.H1 = H[1][:, perm, :]
        self.lpsi0, self.lpsi1 = lpsi[0], lpsi[1]
        self.R, self.B = R, B
        rat = self.H1 - self.H0
        def suf(a, fn):
            out = np.zeros((R, B + 1))
            for k in range(B - 1, -1, -1):
                out[:, k] = out[:, k + 1] + fn(a[:, k, 0], a[:, k, 1])
            return out
        self.sLo = suf(self.H0, np.minimum)
        self.sHi = suf(self.H0, np.maximum)
        self.sRlo = suf(rat, np.minimum)
        self.sRhi = suf(rat, np.maximum)

    # T1 bound at a cube
    def bound(self, k, LD, LN, side):
        ld = self.lpsi0 + LD
        lr = (self.lpsi1 + LN) - ld
        lU = ld + self.sHi[:, k]
        s = lU.max()
        Lm = np.exp(ld + self.sLo[:, k] - s)
        Um = np.exp(lU - s)
        if side == "lower":
            r = np.exp(lr + self.sRlo[:, k])
            idx = np.argsort(r, kind="stable")
            r, Lm, Um = r[idx], Lm[idx], Um[idx]
            cu_n = np.concatenate([[0.0], np.cumsum(r * Um)])
            cu_d = np.concatenate([[0.0], np.cumsum(Um)])
            tl_n = np.concatenate([np.cumsum((r * Lm)[::-1])[::-1], [0.0]])
            tl_d = np.concatenate([np.cumsum(Lm[::-1])[::-1], [0.0]])
            return float(np.min((cu_n + tl_n) / (cu_d + tl_d)))
        r = np.exp(lr + self.sRhi[:, k])
        idx = np.argsort(-r, kind="stable")
        r, Lm, Um = r[idx], Lm[idx], Um[idx]
        cu_n = np.concatenate([[0.0], np.cumsum(r * Um)])
        cu_d = np.concatenate([[0.0], np.cumsum(Um)])
        tl_n = np.concatenate([np.cumsum((r * Lm)[::-1])[::-1], [0.0]])
        tl_d = np.concatenate([np.cumsum(Lm[::-1])[::-1], [0.0]])
        return float(np.max((cu_n + tl_n) / (cu_d + tl_d)))


class Side:
    """Branch-and-bound state for one side (lower or upper) of Q."""

    def __init__(self, prep: Prepared, side: str):
        self.p, self.side = prep, side
        self.sign = 1.0 if side == "lower" else -1.0  # heap key = sign*bound (min-heap on worst)
        z = np.zeros(prep.R)
        b = prep.bound(0, z, z, side)
        self.heap = [(self.sign * b, 0, ())]
        self.nleaf = 1
        self.internal = set()
        self.cnt = 1
        self.leafinfo = {(): (z, z, b)}
        self.exact_root = prep.R == 1  # single region => root bound is already exact

    @property
    def bound(self):
        return self.sign * self.heap[0][0]

    def done(self, budget):
        if self.exact_root:
            return True
        key, _, path = self.heap[0]
        return len(path) == self.p.B or self.nleaf + 1 > budget

    def step(self, budget):
        if self.done(budget):
            return False
        key, _, path = heapq.heappop(self.heap)
        LD, LN, _ = self.leafinfo.pop(path)
        k = len(path)
        self.internal.add(path)
        for tau in (0, 1):
            LD2 = LD + self.p.H0[:, k, tau]
            LN2 = LN + self.p.H1[:, k, tau]
            b = self.p.bound(k + 1, LD2, LN2, self.side)
            self.cnt += 1
            np_ = path + (tau,)
            self.leafinfo[np_] = (LD2, LN2, b)
            heapq.heappush(self.heap, (self.sign * b, self.cnt, np_))
        self.nleaf += 1
        return True

    def tree(self):
        leaves = set(self.leafinfo)
        def build(pre):
            if pre in leaves:
                return 0
            return [build(pre + (0,)), build(pre + (1,))]
        return build(())

    def certificate(self, spec, lam_claim: Fraction, site, ph):
        return {"spec_hash": spec_hash(spec), "site": site, "interior": list(self.p.S), "ph": str(ph),
                "side": self.side, "lambda": str(lam_claim), "order": [int(t) for t in self.p.order],
                "tree": self.tree()}

    def claim(self):
        b = self.bound
        f = (1 - MARGIN) if self.side == "lower" else (1 + MARGIN)
        return Fraction(b * f)


def adjacency(spec):
    nbr = [dict() for _ in spec["sites"]]
    for key, v in spec["pairs_milli"].items():
        a, b = map(int, key.split(","))
        nbr[a][b] = v
        nbr[b][a] = v
    return nbr


def interior(spec, i, radius, cap=9):
    xyz = [s["xyz"] for s in spec["sites"]]
    cand = sorted((math.dist(xyz[i], xyz[j]), j) for j in range(len(xyz)) if j != i)
    return [i] + [j for d, j in cand if d <= radius][: cap - 1]


def run_bound(prep: Prepared, side: str, budget: int, target: float | None = None):
    s = Side(prep, side)
    while True:
        if target is not None:
            if (side == "lower" and s.bound > target * (1 + MARGIN)) or (side == "upper" and s.bound < target * (1 - MARGIN)):
                break
        if not s.step(budget):
            break
    return s


def decide(prep: Prepared, budget: int):
    """Certified label at Q=1: 'A' (Q>=1, protonated-majority), 'B' (Q<=1) or 'U'. Returns (label, Side)."""
    lo, up = Side(prep, "lower"), Side(prep, "upper")
    while True:
        if lo.bound > 1 + MARGIN:
            return "A", lo
        if up.bound < 1 - MARGIN:
            return "B", up
        a = lo.step(budget)
        b = up.step(budget)
        if not a and not b:
            return "U", None


def theta_of(q):
    return q / (1.0 + q)
