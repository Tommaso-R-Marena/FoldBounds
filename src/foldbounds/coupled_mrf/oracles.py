"""Exact partition-function oracles, implemented independently of the producer.

* ``brute_force_z``   -- enumerates every joint state (tiny models only).
* ``ve_z``            -- dict-based variable elimination with a min-fill order
                         (moderate treewidth); separate code from tables.py.
* ``grid_transfer_z`` -- row-by-row ("broken row") transfer-matrix DP for
                         binary pairwise models on an L x M grid with
                         arbitrary per-site and per-edge rational tables.

All three return exact Fractions.  Factor tables are scaled to integers so the
inner loops are integer arithmetic; the scale is restored exactly at the end.
"""
from __future__ import annotations

import re
from fractions import Fraction
from itertools import product
from math import lcm, prod


def _int_table(table):
    den = 1
    for t in table:
        den = lcm(den, Fraction(t).denominator)
    return [int(Fraction(t) * den) for t in table], Fraction(1, den)


def brute_force_z(model, which: str) -> Fraction:
    view = model.system(which)
    variables = list(model.variables)
    dims = [view.domains[v] for v in variables]
    if prod(dims) > 2 ** 22:
        raise ValueError("model too large for brute force")
    pos = {v: i for i, v in enumerate(variables)}
    facs = []
    scale = Fraction(1)
    for f in view.factors:
        vals, s = _int_table(f.table)
        scale *= s
        idx = [pos[v] for v in f.scope]
        fd = [view.domains[v] for v in f.scope]
        facs.append((idx, fd, vals))
    total = 0
    for x in product(*(range(d) for d in dims)):
        w = 1
        for idx, fd, vals in facs:
            k = 0
            for i, d in zip(idx, fd):
                k = k * d + x[i]
            w *= vals[k]
            if w == 0:
                break
        total += w
    return scale * total


def _min_fill_order(scopes, variables):
    adj = {v: set() for v in variables}
    for s in scopes:
        for a in s:
            adj[a].update(b for b in s if b != a)
    order = []
    remaining = set(variables)
    while remaining:
        best, best_key = None, None
        for v in remaining:
            nb = list(adj[v])
            fill = sum(1 for i in range(len(nb)) for j in range(i + 1, len(nb))
                       if nb[j] not in adj[nb[i]])
            key = (fill, len(nb), v)
            if best_key is None or key < best_key:
                best, best_key = v, key
        nb = adj.pop(best)
        for w in nb:
            adj[w].discard(best)
            adj[w] |= nb - {w}
        remaining.discard(best)
        order.append(best)
    return order


def ve_z(model, which: str, max_table: int = 2 ** 22) -> Fraction:
    """Exact Z by dict-based variable elimination (min-fill order)."""
    view = model.system(which)
    dom = view.domains
    scale = Fraction(1)
    pool = []
    for f in view.factors:
        vals, s = _int_table(f.table)
        scale *= s
        table = {}
        for k, x in enumerate(product(*(range(dom[v]) for v in f.scope))):
            table[x] = vals[k]
        pool.append((tuple(f.scope), table))
    order = _min_fill_order([s for s, _ in pool], list(model.variables))
    for v in order:
        touching = [p for p in pool if v in p[0]]
        pool = [p for p in pool if v not in p[0]]
        if not touching:
            pool.append(((), {(): dom[v]}))
            continue
        scope = []
        for s, _ in touching:
            for w in s:
                if w != v and w not in scope:
                    scope.append(w)
        if prod(dom[w] for w in scope) * dom[v] > max_table:
            raise ValueError("elimination table too large for the VE oracle")
        out = {}
        full = scope + [v]
        where = [[full.index(w) for w in s] for s, _ in touching]
        tabs = [t for _, t in touching]
        for x in product(*(range(dom[w]) for w in full)):
            w_ = 1
            for wi, t in zip(where, tabs):
                w_ *= t.get(tuple(x[i] for i in wi), 0)
                if w_ == 0:
                    break
            if w_:
                key = x[:-1]
                out[key] = out.get(key, 0) + w_
        pool.append((tuple(scope), out))
    total = 1
    for s, t in pool:
        total *= t.get((), 0)
    return scale * total


_SITE = re.compile(r"^x(\d+)_(\d+)$")


def grid_transfer_z(model, which: str) -> Fraction:
    """Exact Z for a binary pairwise model whose variables are named x{r}_{c}.

    Every variable must be binary and every factor unary or a nearest-
    neighbour pair on the grid.  Sites are processed in row-major order; the
    DP state is the "broken row" of the M most recent spins.  Cost
    O(L * M * 2^M) big-integer operations.
    """
    view = model.system(which)
    coords = {}
    for v in model.variables:
        m = _SITE.match(v)
        if not m or view.domains[v] != 2:
            raise ValueError("grid oracle needs binary variables named x{r}_{c}")
        coords[v] = (int(m.group(1)), int(m.group(2)))
    L = 1 + max(r for r, _ in coords.values())
    M = 1 + max(c for _, c in coords.values())
    if len(coords) != L * M:
        raise ValueError("grid oracle needs a full rectangular grid")
    unary = {rc: [1, 1] for rc in coords.values()}
    up = {}      # (r, c) -> 2x2 integer table indexed [up_spin][spin]
    left = {}    # (r, c) -> 2x2 integer table indexed [left_spin][spin]
    scale = Fraction(1)
    const = 1
    for f in view.factors:
        vals, s = _int_table(f.table)
        scale *= s
        if len(f.scope) == 0:
            const *= vals[0]
        elif len(f.scope) == 1:
            rc = coords[f.scope[0]]
            unary[rc] = [unary[rc][0] * vals[0], unary[rc][1] * vals[1]]
        elif len(f.scope) == 2:
            a, b = coords[f.scope[0]], coords[f.scope[1]]
            t = [[vals[0], vals[1]], [vals[2], vals[3]]]
            if a > b:
                a, b = b, a
                t = [[t[0][0], t[1][0]], [t[0][1], t[1][1]]]
            if b == (a[0] + 1, a[1]):
                d = up
            elif b == (a[0], a[1] + 1):
                d = left
            else:
                raise ValueError("grid oracle: factor is not a nearest-neighbour pair")
            cur = d.get(b, [[1, 1], [1, 1]])
            d[b] = [[cur[i][j] * t[i][j] for j in range(2)] for i in range(2)]
        else:
            raise ValueError("grid oracle supports unary and pairwise factors only")
    n = 1 << M
    f = [0] * n
    f[0] = 1
    for r in range(L):
        for c in range(M):
            bit = 1 << c
            U = unary[(r, c)]
            V = up.get((r, c), [[1, 1], [1, 1]]) if r > 0 else None
            H = left.get((r, c), [[1, 1], [1, 1]]) if c > 0 else None
            new = [0] * n
            lbit = bit >> 1
            for s in range(n):
                t = 1 if s & bit else 0
                base = s & ~bit
                if V is None:
                    acc = f[base]           # row 0: the stored bit is still the initial 0
                else:
                    acc = f[base] * V[0][t] + f[base | bit] * V[1][t]
                if not acc:
                    continue
                acc *= U[t]
                if H is not None:
                    acc *= H[1 if s & lbit else 0][t]
                new[s] = acc
            f = new
    return scale * const * sum(f)
