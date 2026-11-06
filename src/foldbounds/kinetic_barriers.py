"""Untrusted proposal of exact finite 2D HP kinetic cut certificates.

The verifier lives in kinetic_barrier_verifier.py and DOES NOT use the
proposal's walk, contact, or move generators. This is a lattice toy model.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

Point = tuple[int, int]
Path = tuple[Point, ...]
STEPS = ((1, 0), (0, 1), (-1, 0), (0, -1))
ROOT: Path = ((0, 0), (1, 0))


@dataclass(frozen=True)
class BasinCertificate:
    sequence: str
    threshold: int
    seed: Path
    members: tuple[Path, ...]
    mode: str = "superlevel"  # alternative: arbitrary, exact-mincut
    flow_q: Fraction | None = None
    flow_edges: tuple[tuple[Path, Path, int], ...] = ()


class ProposalLimit(ValueError):
    """The producer did not enumerate the full fixed-root ensemble."""


def _walks(n: int, cap: int) -> list[Path]:
    found: list[Path] = []

    def visit(p: Path, occupied: set[Point]) -> None:
        if len(p) == n:
            found.append(p)
            if len(found) > cap:
                raise ProposalLimit("walk budget exceeded; no partial certificate")
            return
        x, y = p[-1]
        for dx, dy in STEPS:
            u = (x + dx, y + dy)
            if u not in occupied:
                visit(p + (u,), occupied | {u})

    visit(ROOT, set(ROOT))
    return found


def _hh_score(p: Path, hydrophobic: frozenset[int]) -> int:
    occupied = {pos: i for i, pos in enumerate(p)}
    out = 0
    for i in hydrophobic:
        x, y = p[i]
        for dx, dy in STEPS:
            j = occupied.get((x + dx, y + dy), -1)
            if j >= i + 3 and j in hydrophobic:
                out += 1
    return out


def _local_moves(p: Path):
    occupied = set(p)
    for i in range(2, len(p)):
        x, y = p[i - 1]
        for dx, dy in STEPS:
            v = (x + dx, y + dy)
            if v in occupied or v == p[i]:
                continue
            if i + 1 < len(p):
                w = p[i + 1]
                if abs(v[0] - w[0]) + abs(v[1] - w[1]) != 1:
                    continue
            yield p[:i] + (v,) + p[i + 1:]


def propose_kinetic_basin(
    sequence: str,
    q: Fraction,
    *,
    allow_closed: bool = False,
    min_barrier: int = 1,
    max_states: int = 40000,
) -> BasinCertificate:
    """Enumerate finite cut components and propose the strongest q-weighted cut.

    Membership is the complete connected component in {score >= threshold}.
    It may be metastable (positive boundary) or a permanently closed class.
    This is an UNTRUSTED exhaustive producer. It never changes FoldBounds' core.
    """
    if not isinstance(sequence, str) or len(sequence) < 4 or set(sequence) - {"H", "P"}:
        raise ValueError("use a length >= 4 H/P sequence")
    if not isinstance(q, Fraction) or q < 1 or type(min_barrier) is not int or min_barrier < 1:
        raise ValueError("q must be an exact Fraction >= 1; min_barrier positive integer")
    if type(max_states) is not int or max_states < 1:
        raise ValueError("invalid state budget")
    walks = _walks(len(sequence), max_states)
    ix = {p: i for i, p in enumerate(walks)}
    hydro = frozenset(i for i, a in enumerate(sequence) if a == "H")
    score = [_hh_score(p, hydro) for p in walks]
    target = max(score)
    neighbors = [tuple(ix[m] for m in _local_moves(p)) for p in walks]
    proposal_slots = 4 * (len(sequence) - 2)
    best: tuple[tuple, BasinCertificate] | None = None
    for threshold in range(1, target):
        visited: set[int] = set()
        for start in range(len(walks)):
            if score[start] < threshold or start in visited:
                continue
            visited.add(start)
            part = [start]
            for u in part:
                for v in neighbors[u]:
                    if v not in visited and score[v] >= threshold:
                        visited.add(v)
                        part.append(v)
            high = max(score[i] for i in part)
            if high == target:
                continue
            members = set(part)
            boundary = [(u, v) for u in part for v in neighbors[u] if v not in members]
            if not boundary and not allow_closed:
                continue
            cut_score = max((min(score[u], score[v]) for u, v in boundary), default=-1)
            exponent = high - cut_score if boundary else high + 1
            if exponent < min_barrier:
                continue
            z_part = sum((q ** score[u] for u in part), Fraction(0))
            flux = sum((q ** min(score[u], score[v]) for u, v in boundary), Fraction(0))
            leak = flux / (proposal_slots * z_part)
            # Closed classes are ranked first if explicitly allowed; otherwise
            # maximize the rigorous all-q exponent, then minimize exact leakage.
            ranking = (int(not boundary), exponent, -leak, -len(part))
            cert = BasinCertificate(sequence, threshold, walks[min(part)],
                                    tuple(sorted((walks[i] for i in part))))
            if best is None or ranking > best[0]:
                best = (ranking, cert)
    if best is None:
        raise ProposalLimit("no suitable nonnative kinetic cut found within limit")
    return best[1]


def propose_optimal_kinetic_cut(
    sequence: str,
    q: Fraction,
    *,
    source: Path | None = None,
    max_states: int = 40000,
    max_augmentations: int = 2000,
) -> BasinCertificate:
    """Propose a rational exact s-to-native minimum cut with a flow witness.

    Uses a pure-Python integer residual-network augmenting-path algorithm.
    Denominators are cleared with denominator(q)**max_contact_score. The
    returned net-edge flow is independently audited by the verifier against
    conservation, exact capacities, and cut/flow equality. No trust is placed
    in this producer's claimed minimum-cut optimization.
    """
    if not isinstance(sequence, str) or len(sequence) < 4 or set(sequence) - {"H", "P"}:
        raise ValueError("invalid HP sequence")
    if not isinstance(q, Fraction) or q < 1:
        raise ValueError("exact rational q>=1 required")
    if type(max_states) is not int or max_states < 1 or type(max_augmentations) is not int or max_augmentations < 1:
        raise ValueError("invalid work budget")
    walks = _walks(len(sequence), max_states)
    index = {p: i for i, p in enumerate(walks)}
    h = frozenset(i for i, r in enumerate(sequence) if r == "H")
    scores = [_hh_score(p, h) for p in walks]
    target = max(scores)
    if target == 0:
        raise ProposalLimit("no nonnative score class exists")
    if source is None:
        # Select a high-score kinetic trap as the cut's fixed deterministic
        # starting state. Min-cut then optimizes the boundary around it.
        baseline = propose_kinetic_basin(sequence, q, max_states=max_states)
        source = min((p for p in baseline.members
                      if _hh_score(p, h) == max(_hh_score(z, h) for z in baseline.members)))
    if source not in index:
        raise ValueError("starting path must belong to canonical full ensemble")
    source_i = index[source]
    if scores[source_i] == target:
        raise ValueError("source is already a native maximum-contact state")

    scale = q.denominator ** target
    N = len(walks)
    sink = N
    residual: list[dict[int, int]] = [{} for _ in range(N + 1)]
    original_edges: list[tuple[int, int, int]] = []
    total_c = 0
    for u, p in enumerate(walks):
        for new in _local_moves(p):
            v = index[new]
            if u >= v:
                continue
            cap = scale * q**min(scores[u], scores[v])
            if cap.denominator != 1:
                raise AssertionError("integer capacity scaling failed")
            c = int(cap)
            original_edges.append((u, v, c))
            residual[u][v] = c
            residual[v][u] = c
            total_c += c
    inf = total_c + 1
    for i, score in enumerate(scores):
        if score == target:
            residual[i][sink] = inf
            residual[sink][i] = 0
    from collections import deque
    for attempt in range(max_augmentations + 1):
        parent = {source_i: -1}
        queue = deque([source_i])
        while queue and sink not in parent:
            u = queue.popleft()
            for v, cap in residual[u].items():
                if cap > 0 and v not in parent:
                    parent[v] = u
                    queue.append(v)
        if sink not in parent:
            break
        if attempt == max_augmentations:
            raise ProposalLimit("max-flow augment budget exceeded; no partial certificate")
        v = sink
        bottleneck = None
        while v != source_i:
            u = parent[v]
            cap = residual[u][v]
            bottleneck = cap if bottleneck is None else min(bottleneck, cap)
            v = u
        assert bottleneck is not None and bottleneck > 0
        v = sink
        while v != source_i:
            u = parent[v]
            residual[u][v] -= bottleneck
            residual[v][u] += bottleneck
            v = u
    else:
        raise ProposalLimit("failed to certify residual termination")

    seen = {source_i}
    q_walk = deque([source_i])
    while q_walk:
        u = q_walk.popleft()
        for v, cap in residual[u].items():
            if cap > 0 and v not in seen:
                seen.add(v)
                q_walk.append(v)
    if sink in seen or any(scores[i] == target for i in seen):
        raise AssertionError("maxflow residual mincut has reached sink")
    flow_edges = tuple((walks[u], walks[v], c-residual[u][v])
                       for u, v, c in original_edges if c != residual[u][v])
    return BasinCertificate(
        sequence=sequence, threshold=0, seed=source,
        members=tuple(sorted(walks[i] for i in seen)),
        mode="exact-mincut", flow_q=q, flow_edges=flow_edges,
    )
