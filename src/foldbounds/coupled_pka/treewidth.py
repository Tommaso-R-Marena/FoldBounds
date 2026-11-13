"""Greedy min-fill elimination ordering: upper bound on treewidth of the interaction graph."""
from __future__ import annotations


def min_fill_width(n: int, edges) -> tuple[int, list[int]]:
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    alive = set(range(n))
    width, order = 0, []
    while alive:
        best, bestfill = None, None
        for v in alive:
            nb = list(adj[v])
            fill = sum(1 for x in range(len(nb)) for y in range(x + 1, len(nb)) if nb[y] not in adj[nb[x]])
            key = (fill, len(nb), v)
            if bestfill is None or key < bestfill:
                best, bestfill = v, key
        v = best
        nb = list(adj[v])
        width = max(width, len(nb))
        for x in nb:
            for y in nb:
                if x != y:
                    adj[x].add(y)
        for x in nb:
            adj[x].discard(v)
        alive.discard(v)
        adj[v] = set()
        order.append(v)
    return width, order


def spec_width(spec):
    edges = [tuple(map(int, k.split(","))) for k in spec["pairs_milli"]]
    return min_fill_width(len(spec["sites"]), edges)[0]
