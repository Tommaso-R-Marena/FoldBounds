#!/usr/bin/env python3
"""UNTRUSTED certificate producers for ``FOLDBOUNDS-KERNEL 1``.

Nothing here is trusted: a producer only chooses the *shape* of the tree
(which prefixes are split, bounded, enumerated).  The checker recomputes every
number, so a buggy or malicious producer can at worst emit a rejected or a
loose certificate.

Producers (all deterministic, all return the tree as a ``str`` over
``N E X B -``):

``uniform(seq, depth, enum_r=None)``
    split every legal branch until ``depth`` residues are placed; leaves are
    ``B``, or ``X`` when the remaining count ``r <= enum_r``; complete walks
    are ``E``; dead ends are ``N----``.
``hybrid(seq, depth, kappa)``
    uniform depth; a leaf is ``X`` when its prefix already has ``>= kappa``
    H-H contacts, else ``B``.
``adaptive(seq, q, t, budget)``
    best-first refinement of the bound leaf contributing most to the certified
    width ``pHi - pLo`` of the target ``(q, t)``.
``witness_paths(seq, q, t, ...)``
    spines towards a few conformations with ``k >= t`` (``E`` leaves), every
    other legal branch bounded.
``random_tree(seq, rng, ...)``
    random valid shapes for differential testing.

CLI::

    python producer.py --seq HPHPPHHPHH --query 2 4 --query 3/2 2 \\
        --mode uniform|hybrid|adaptive|witness --depth 8 --kappa 1 \\
        --budget 100000 --out FILE
"""

from __future__ import annotations

import argparse
import heapq
import sys
from fractions import Fraction
from pathlib import Path
from typing import Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

from reference_kernel import bound_entry  # noqa: E402  (scoring of leaves only)

FORMAT = "FOLDBOUNDS-KERNEL 1"
DIRS = ((1, 0), (-1, 0), (0, 1), (0, -1))      # 0 = +x, 1 = -x, 2 = +y, 3 = -y
ADAPTIVE_ENUM_R = 8

Site = tuple[int, int]


def _check_seq(seq: str) -> int:
    if not isinstance(seq, str) or len(seq) < 2 or any(ch not in "HP" for ch in seq):
        raise ValueError("sequence must be a string over {H, P} of length >= 2")
    return len(seq)


def _new_contacts(seq: str, occ: dict[Site, int], site: Site, j: int) -> int:
    """Contacts created by placing residue ``j`` on ``site``."""
    if seq[j] != "H":
        return 0
    x, y = site
    count = 0
    for dx, dy in DIRS:
        i = occ.get((x + dx, y + dy))
        if i is not None and i <= j - 2 and seq[i] == "H":
            count += 1
    return count


def _with_recursion(n: int, fn):
    old = sys.getrecursionlimit()
    if old < 4 * n + 200:
        sys.setrecursionlimit(4 * n + 200)
    try:
        return fn()
    finally:
        sys.setrecursionlimit(old)


# ---------------------------------------------------------------------------
# uniform / hybrid
# ---------------------------------------------------------------------------

def _fixed_depth(seq: str, depth: int, leaf_kind) -> str:
    n = _check_seq(seq)
    out: list[str] = []
    occ: dict[Site, int] = {(0, 0): 0, (1, 0): 1}

    def rec(x: int, y: int, m: int, a: int) -> None:
        r = n - m
        if r == 0:
            out.append("E")
            return
        if m >= depth:
            out.append(leaf_kind(r, a))
            return
        out.append("N")
        for dx, dy in DIRS:
            p = (x + dx, y + dy)
            if p in occ:
                out.append("-")
                continue
            gained = _new_contacts(seq, occ, p, m)
            occ[p] = m
            rec(p[0], p[1], m + 1, a + gained)
            del occ[p]

    _with_recursion(n, lambda: rec(1, 0, 2, 0))
    return "".join(out)


def uniform(seq: str, depth: int, enum_r: Optional[int] = None) -> str:
    """Uniform-depth frontier: ``B`` leaves, ``X`` when ``r <= enum_r``."""
    return _fixed_depth(
        seq, depth,
        lambda r, a: "X" if enum_r is not None and r <= enum_r else "B")


def hybrid(seq: str, depth: int, kappa: int) -> str:
    """Uniform depth; leaves with ``>= kappa`` known contacts are enumerated."""
    return _fixed_depth(seq, depth, lambda r, a: "X" if a >= kappa else "B")


# ---------------------------------------------------------------------------
# adaptive
# ---------------------------------------------------------------------------

def _x_cost(r: int) -> int:
    """Estimated conformations under an ``X`` leaf: ``ceil(2.7^r)`` (exact ints)."""
    return -((-(27 ** r)) // (10 ** r))


def _completion_histogram(seq: str, sites: Sequence[Site], a: int) -> dict[int, int]:
    """Contact histogram of all completions of ``sites`` (producer-side only)."""
    n = len(seq)
    occ = {s: i for i, s in enumerate(sites)}
    hist: dict[int, int] = {}

    def rec(x: int, y: int, m: int, k: int) -> None:
        if m == n:
            hist[k] = hist.get(k, 0) + 1
            return
        for dx, dy in DIRS:
            p = (x + dx, y + dy)
            if p in occ:
                continue
            gained = _new_contacts(seq, occ, p, m)
            occ[p] = m
            rec(p[0], p[1], m + 1, k + gained)
            del occ[p]

    rec(sites[-1][0], sites[-1][1], len(sites), a)
    return hist


class _Leaf:
    __slots__ = ("sites", "entry", "r", "parent", "slot", "alive", "seqno")

    def __init__(self, sites, entry, r, parent, slot, seqno):
        self.sites = sites
        self.entry = entry
        self.r = r
        self.parent = parent          # children list of the parent split (or None at the root)
        self.slot = slot
        self.alive = True
        self.seqno = seqno


def adaptive(seq: str, q, t: int, budget: int, enum_r: int = ADAPTIVE_ENUM_R) -> str:
    """Best-first refinement for the target query ``(q, t)``.

    Start from the root ``B``.  Repeatedly take the bound leaf with the largest
    contribution to the certified width ``pHi - pLo`` and replace it by a split
    ``N`` with ``B`` children, by ``X`` when ``1 <= r <= enum_r`` (and the
    enumeration still fits in the budget; otherwise it is split), or by ``E``
    when ``r = 0``.  Stop before the estimated checker work (1 unit per tree
    byte plus ``ceil(2.7^r)`` units per ``X`` leaf) would exceed ``budget``,
    or when no leaf contributes to the width any more.

    Contribution of a leaf: with ``s_e = eHi - eLo`` and ``s_n = nHi - nLo``
    restricted to the leaf, ``s_e * (nLo + nHi) + s_n * (eLo + eHi)``; summed
    over all leaves this is exactly twice the numerator of ``pHi - pLo``
    (``eHi*nHi - eLo*nLo``).  Ties, and the degenerate situation where every
    weighted contribution vanishes although the width is positive, are broken
    by the raw slack ``s_e + s_n`` and then by creation order.  Priorities are
    refreshed periodically because the global sums move.  All arithmetic is
    exact integer arithmetic, so the result is deterministic.
    """
    n = _check_seq(seq)
    q = Fraction(q)
    if q < 1:
        raise ValueError("q must be >= 1")
    if not isinstance(t, int) or t < 0:
        raise ValueError("t must be a natural number")
    top = 2 * n + 2
    p_num, p_den = q.numerator, q.denominator
    weight = [p_num ** k * p_den ** (top - k) for k in range(top + 1)]   # q^k * den^top

    def parts(entry) -> tuple[int, int, int, int]:
        a, cap, cnt, fam = entry
        e_lo = fam * weight[a] if t <= a else 0
        e_hi = cnt * weight[cap] if t <= cap else 0
        n_lo = fam * weight[a] if cap < t else 0
        n_hi = cnt * weight[min(cap, t - 1)] if a < t else 0
        return e_lo, e_hi, n_lo, n_hi

    totals = [0, 0, 0, 0]              # eLo, eHi, nLo, nHi (scaled by den^top)

    def key(leaf: _Leaf) -> tuple[int, int, int, _Leaf]:
        e_lo, e_hi, n_lo, n_hi = parts(leaf.entry)
        s_e = e_hi - e_lo
        s_n = n_hi - n_lo
        weighted = s_e * (totals[2] + totals[3]) + s_n * (totals[0] + totals[1])
        return (-weighted, -(s_e + s_n), leaf.seqno, leaf)

    def width_is_zero() -> bool:
        e_lo, e_hi, n_lo, n_hi = totals
        if e_hi == 0:
            return True                 # pHi = 0, hence pLo = 0
        if e_lo == 0:
            return False                # pLo = 0 < pHi
        if n_lo == 0:
            return n_hi == 0            # pHi = 1; pLo = 1 iff nHi = 0
        return e_lo * (e_hi + n_lo) == e_hi * (e_lo + n_hi)

    counter = 0
    root_sites = ((0, 0), (1, 0))
    root_holder: list = [None]
    root = _Leaf(root_sites, bound_entry(seq, root_sites), n - 2, root_holder, 0, counter)
    root_holder[0] = root
    for idx, v in enumerate(parts(root.entry)):
        totals[idx] += v
    live: list[_Leaf] = [root]
    heap = [key(root)]
    work = 1
    since_refresh = 0

    def refresh() -> None:
        nonlocal live, heap, since_refresh
        live = [lf for lf in live if lf.alive]
        heap = [key(lf) for lf in live]
        heapq.heapify(heap)
        since_refresh = 0

    while heap:
        if since_refresh >= max(16, len(heap) // 8):
            refresh()
        leaf = heapq.heappop(heap)[3]
        if key(leaf)[0] == 0:
            # Possibly stale: recompute every priority (the popped leaf is
            # still alive, so it is re-inserted) and look at the true maximum.
            refresh()
            if heap[0][0] == 0 and (heap[0][1] == 0 or width_is_zero()):
                break                   # nothing contributes to pHi - pLo
            leaf = heapq.heappop(heap)[3]
        since_refresh += 1
        if leaf.r == 0:
            cost = 0
            new_node: object = "E"
        elif leaf.r <= enum_r and work + _x_cost(leaf.r) <= budget:
            cost = _x_cost(leaf.r)
            new_node = "X"
        else:
            # Large subtree, or a small one whose enumeration is not affordable
            # any more: split one level instead.
            cost = 4
            new_node = None
        if work + cost > budget:
            break
        work += cost
        leaf.alive = False
        for idx, v in enumerate(parts(leaf.entry)):
            totals[idx] -= v
        if new_node is not None:
            if leaf.r == 0:
                hist = {leaf.entry[0]: 1}
            else:
                hist = _completion_histogram(seq, leaf.sites, leaf.entry[0])
            for k, count in hist.items():
                mass = count * weight[k]
                if k >= t:
                    totals[0] += mass
                    totals[1] += mass
                else:
                    totals[2] += mass
                    totals[3] += mass
            leaf.parent[leaf.slot] = new_node
            continue
        children: list = [None, None, None, None]
        occupied = set(leaf.sites)
        x, y = leaf.sites[-1]
        for d, (dx, dy) in enumerate(DIRS):
            p = (x + dx, y + dy)
            if p in occupied:
                children[d] = "-"
                continue
            sites = leaf.sites + (p,)
            counter += 1
            child = _Leaf(sites, bound_entry(seq, sites), leaf.r - 1, children, d, counter)
            children[d] = child
            for idx, v in enumerate(parts(child.entry)):
                totals[idx] += v
            live.append(child)
            heapq.heappush(heap, key(child))
        leaf.parent[leaf.slot] = children

    out: list[str] = []
    stack: list = [root_holder[0]]
    while stack:
        node = stack.pop()
        if isinstance(node, str):
            out.append(node)
        elif isinstance(node, _Leaf):
            out.append("B")
        else:
            out.append("N")
            stack.extend(reversed(node))
    return "".join(out)


# ---------------------------------------------------------------------------
# witness paths
# ---------------------------------------------------------------------------

def find_witnesses(seq: str, t: int, max_witnesses: int = 3,
                   node_limit: int = 200_000) -> list[tuple[Site, ...]]:
    """Up to ``max_witnesses`` conformations with ``k >= t`` (depth-first order).

    Falls back to the best conformation seen when the threshold is not reached
    within ``node_limit`` search nodes.  Heuristic and untrusted.
    """
    n = _check_seq(seq)
    occ: dict[Site, int] = {(0, 0): 0, (1, 0): 1}
    walk: list[Site] = [(0, 0), (1, 0)]
    found: list[tuple[Site, ...]] = []
    best: list = [-1, None]
    nodes = [0]

    def rec(a: int) -> bool:
        nodes[0] += 1
        m = len(walk)
        if m == n:
            if a > best[0]:
                best[0], best[1] = a, tuple(walk)
            if a >= t:
                found.append(tuple(walk))
            return len(found) >= max_witnesses or nodes[0] >= node_limit
        if bound_entry(seq, walk)[1] < t and best[1] is not None:
            return nodes[0] >= node_limit
        x, y = walk[-1]
        for dx, dy in DIRS:
            p = (x + dx, y + dy)
            if p in occ:
                continue
            gained = _new_contacts(seq, occ, p, m)
            occ[p] = m
            walk.append(p)
            stop = rec(a + gained)
            walk.pop()
            del occ[p]
            if stop:
                return True
        return nodes[0] >= node_limit

    _with_recursion(n, lambda: rec(0))
    if found:
        return found
    return [best[1]] if best[1] is not None else []


def spine_tree(seq: str, walks: Sequence[Sequence[Site]]) -> str:
    """Tree whose split nodes are the proper prefixes of ``walks``.

    Complete walks end in ``E``; every legal branch off the spines is ``B``.
    """
    n = _check_seq(seq)
    trie: dict = {}
    for walk in walks:
        if len(walk) != n or tuple(walk[:2]) != ((0, 0), (1, 0)):
            raise ValueError("witnesses must be complete anchored walks")
        node = trie
        for site in walk[2:]:
            node = node.setdefault(tuple(site), {})
    out: list[str] = []
    occ = {(0, 0), (1, 0)}

    def rec(x: int, y: int, m: int, node: dict) -> None:
        if m == n:
            out.append("E")
            return
        out.append("N")
        for dx, dy in DIRS:
            p = (x + dx, y + dy)
            if p in occ:
                out.append("-")
            elif p in node:
                occ.add(p)
                rec(p[0], p[1], m + 1, node[p])
                occ.remove(p)
            else:
                out.append("B")

    if not trie and n > 2:
        return "B"
    _with_recursion(n, lambda: rec(1, 0, 2, trie))
    return "".join(out)


def witness_paths(seq: str, q, t: int, max_witnesses: int = 3,
                  node_limit: int = 200_000) -> str:
    """Spines to a few witnesses of ``k >= t``; ``q`` is accepted for symmetry
    with the other producers but the shape does not depend on it."""
    walks = find_witnesses(seq, t, max_witnesses=max_witnesses, node_limit=node_limit)
    if not walks:
        return "B"
    return spine_tree(seq, walks)


# ---------------------------------------------------------------------------
# random shapes (differential testing)
# ---------------------------------------------------------------------------

def random_prefix(n: int, m: int, rng) -> tuple[Site, ...]:
    """A random anchored self-avoiding walk with ``m`` sites (``2 <= m <= n``)."""
    if not 2 <= m <= n:
        raise ValueError("need 2 <= m <= n")
    while True:
        walk: list[Site] = [(0, 0), (1, 0)]
        used = {(0, 0), (1, 0)}
        while len(walk) < m:
            x, y = walk[-1]
            free = [(x + dx, y + dy) for dx, dy in DIRS if (x + dx, y + dy) not in used]
            if not free:
                break
            nxt = rng.choice(free)
            walk.append(nxt)
            used.add(nxt)
        if len(walk) == m:
            return tuple(walk)


def isolate_prefix(seq: str, prefix: Sequence[Site], leaf: str = "B", other: str = "X") -> str:
    """Tree with one ``leaf`` node at the walk ``prefix``; every other legal
    branch off the path to it is ``other``.  With ``other = "X"`` the single
    bound entry of the prefix can be read off the checker's sums."""
    n = _check_seq(seq)
    path = [tuple(s) for s in prefix]
    if len(path) < 2 or len(path) > n or path[0] != (0, 0) or path[1] != (1, 0):
        raise ValueError("prefix must be an anchored walk with 2..n sites")
    if len(set(path)) != len(path) or any(
            abs(a[0] - b[0]) + abs(a[1] - b[1]) != 1 for a, b in zip(path, path[1:])):
        raise ValueError("prefix must be a self-avoiding nearest-neighbour walk")
    # For every split on the path: the children before and after the on-path child.
    around: list[tuple[list[str], list[str]]] = []
    for m in range(2, len(path)):
        occ = set(path[:m])
        x, y = path[m - 1]
        before: list[str] = []
        after: list[str] = []
        seen_path = False
        for dx, dy in DIRS:
            p = (x + dx, y + dy)
            if p == path[m]:
                seen_path = True
                continue
            (after if seen_path else before).append("-" if p in occ else other)
        around.append((before, after))
    out: list[str] = []
    for before, _ in around:
        out.append("N")
        out.extend(before)
    out.append(leaf)
    for _, after in reversed(around):
        out.extend(after)
    return "".join(out)


def random_tree(seq: str, rng, p_split: float = 0.7, p_enum: float = 0.25,
                max_enum_r: int = 6, max_bytes: int = 4000) -> str:
    """A random valid tree: every node independently splits, bounds or
    enumerates; terminals are mostly ``E`` but sometimes ``X`` or ``B``."""
    n = _check_seq(seq)
    out: list[str] = []
    occ = {(0, 0), (1, 0)}

    def rec(x: int, y: int, m: int) -> None:
        r = n - m
        if r == 0:
            u = rng.random()
            out.append("E" if u < 0.8 else ("X" if u < 0.9 else "B"))
            return
        if len(out) < max_bytes and rng.random() < p_split:
            out.append("N")
            for dx, dy in DIRS:
                p = (x + dx, y + dy)
                if p in occ:
                    out.append("-")
                else:
                    occ.add(p)
                    rec(p[0], p[1], m + 1)
                    occ.remove(p)
            return
        out.append("X" if r <= max_enum_r and rng.random() < p_enum else "B")

    _with_recursion(n, lambda: rec(1, 0, 2))
    return "".join(out)


# ---------------------------------------------------------------------------
# text format
# ---------------------------------------------------------------------------

def format_q(q) -> str:
    """Canonical text of a rational weight (``num`` or ``num/den``)."""
    if isinstance(q, str):
        q = Fraction(q)
    q = Fraction(q)
    if q < 0:
        raise ValueError("q must be nonnegative")
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def to_text(seq: str, queries, tree: str) -> bytes:
    """Certificate file content for ``seq``, ``queries = [(q, t), ...]``, ``tree``."""
    _check_seq(seq)
    queries = list(queries)
    if not queries:
        raise ValueError("at least one query is required")
    lines = [FORMAT, f"sequence {seq}"]
    for q, t in queries:
        t = int(t)
        if t < 0:
            raise ValueError("threshold must be a natural number")
        lines.append(f"query {format_q(q)} {t}")
    lines.append(f"tree {tree}")
    lines.append("end")
    return ("\n".join(lines) + "\n").encode("ascii")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Untrusted FOLDBOUNDS-KERNEL 1 certificate producer")
    ap.add_argument("--seq", required=True, help="HP sequence")
    ap.add_argument("--query", nargs=2, action="append", metavar=("Q", "T"), required=True,
                    help="query: exact rational q >= 1 and threshold t (repeatable; "
                         "adaptive/witness modes target the first one)")
    ap.add_argument("--mode", choices=("uniform", "hybrid", "adaptive", "witness"),
                    default="uniform")
    ap.add_argument("--depth", type=int, default=8, help="placed residues at the frontier")
    ap.add_argument("--enum-r", type=int, default=None,
                    help="uniform mode: enumerate (X) leaves with r <= ENUM_R")
    ap.add_argument("--kappa", type=int, default=1, help="hybrid mode: contact threshold for X")
    ap.add_argument("--budget", type=int, default=100000,
                    help="adaptive mode: estimated checker work units")
    ap.add_argument("--witnesses", type=int, default=3, help="witness mode: number of spines")
    ap.add_argument("--out", required=True, help="output file ('-' for stdout)")
    args = ap.parse_args(argv)
    try:
        queries = [(Fraction(q), int(t)) for q, t in args.query]
        for q, t in queries:
            if q < 1 or t < 0:
                raise ValueError("queries need q >= 1 and t >= 0")
        if args.mode == "uniform":
            tree = uniform(args.seq, args.depth, args.enum_r)
        elif args.mode == "hybrid":
            tree = hybrid(args.seq, args.depth, args.kappa)
        elif args.mode == "adaptive":
            tree = adaptive(args.seq, queries[0][0], queries[0][1], args.budget)
        else:
            tree = witness_paths(args.seq, queries[0][0], queries[0][1],
                                 max_witnesses=args.witnesses)
        data = to_text(args.seq, queries, tree)
    except (ValueError, ZeroDivisionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.out == "-":
        sys.stdout.buffer.write(data)
    else:
        Path(args.out).write_bytes(data)
        print(f"wrote {args.out}: n={len(args.seq)} mode={args.mode} "
              f"tree_bytes={len(tree)} file_bytes={len(data)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
