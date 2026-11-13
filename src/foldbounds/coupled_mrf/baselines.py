"""Separate-bound baseline: exact spanning-tree bounds on Z for pairwise binary MRFs.

For a pairwise model, choose a spanning forest T of the interaction graph.  Every
off-tree factor f (a table on two variables) satisfies  min f <= f <= max f, so

    LB = Z(tree factors, off-tree factors replaced by their minimum constant)
    UB = Z(tree factors, off-tree factors replaced by their maximum constant)

bracket Z exactly (all entries nonnegative).  Both are computed by an exact
tree DP in integers/rationals.  The separate ratio interval is
[LB_B / UB_A, UB_B / LB_A]; it ignores the shared structure between A and B,
which is precisely what the coupled certificate exploits.  (Weitz-type
correlation-decay FPTAS give far better separate bounds in uniqueness
regimes at poly(n) cost; this baseline is the cheap tree bound only.)
"""
from __future__ import annotations

from fractions import Fraction
from math import lcm

from .generators import pairwise_binary


def _int_pair(t):
    den = 1
    for x in t:
        den = lcm(den, Fraction(x).denominator)
    return [int(Fraction(x) * den) for x in t], Fraction(1, den)


def bfs_forest(n, edge_keys, start_order=None):
    adj = [[] for _ in range(n)]
    for (i, j) in edge_keys:
        adj[i].append(j)
        adj[j].append(i)
    parent = [-1] * n
    seen = [False] * n
    order = []           # BFS order, parents before children
    for s in (start_order or range(n)):
        if seen[s]:
            continue
        seen[s] = True
        stack = [s]
        while stack:
            u = stack.pop()
            order.append(u)
            for w in sorted(adj[u]):
                if not seen[w]:
                    seen[w] = True
                    parent[w] = u
                    stack.append(w)
    return parent, order


def tree_bounds(model, which, tree_edges=None):
    """Return (LB, UB) for Z of system ``which`` (exact Fractions).

    ``tree_edges``: optional set of (i, j) index pairs (i<j) forming a forest;
    default is a DFS spanning forest of the A+B interaction graph (same tree for
    both systems so that comparisons are fair).
    """
    variables = list(model.variables)
    const, unary, edges = pairwise_binary(model, which)
    if tree_edges is None:
        ea = pairwise_binary(model, "A")[2].keys()
        eb = pairwise_binary(model, "B")[2].keys()
        keys = sorted(set(ea) | set(eb))
        parent, order = bfs_forest(len(variables), keys)
        tree_edges = {(min(i, p), max(i, p)) for i, p in enumerate(parent) if p >= 0}
    else:
        parent, order = _forest_from_edges(len(variables), tree_edges)
    res = []
    for pick in (min, max):
        scale = Fraction(const)
        for (i, j), t in edges.items():
            if (i, j) not in tree_edges:
                scale *= pick(t)
        # integer tree DP, children before parents
        ui = []
        for u in unary:
            vals, s = _int_pair(u)
            scale *= s
            ui.append(vals)
        msg = [list(u) for u in ui]  # msg[i][x] = weight of subtree of i given x_i = x
        for i in reversed(order):
            p = parent[i]
            if p < 0:
                continue
            key = (min(i, p), max(i, p))
            t = edges.get(key, (Fraction(1),) * 4)
            vals, s = _int_pair(t)
            scale *= s
            # table indexed [x_lo][x_hi]
            def w(xp, xi):
                return vals[xp * 2 + xi] if p < i else vals[xi * 2 + xp]
            new = [sum(w(xp, xi) * msg[i][xi] for xi in (0, 1)) for xp in (0, 1)]
            msg[p] = [msg[p][0] * new[0], msg[p][1] * new[1]]
        # product over forest roots
        z = 1
        for i in range(len(variables)):
            if parent[i] < 0:
                z *= msg[i][0] + msg[i][1]
        res.append(scale * z)
    return res[0], res[1]


def _forest_from_edges(n, tree_edges):
    adj = [[] for _ in range(n)]
    for i, j in tree_edges:
        adj[i].append(j)
        adj[j].append(i)
    parent = [-1] * n
    seen = [False] * n
    order = []
    for s in range(n):
        if seen[s]:
            continue
        seen[s] = True
        stack = [s]
        while stack:
            u = stack.pop()
            order.append(u)
            for w in adj[u]:
                if not seen[w]:
                    seen[w] = True
                    parent[w] = u
                    stack.append(w)
    return parent, order


def separate_ratio_interval(model):
    """[LB_B/UB_A, UB_B/LB_A] from exact spanning-tree bounds (same tree for A and B)."""
    lba, uba = tree_bounds(model, "A")
    lbb, ubb = tree_bounds(model, "B")
    lo = lbb / uba
    hi = ubb / lba if lba > 0 else None
    return lo, hi, (lba, uba, lbb, ubb)
