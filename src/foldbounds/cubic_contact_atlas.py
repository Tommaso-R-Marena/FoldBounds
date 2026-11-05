"""3D cubic HP contact maps with independent geometry and polynomial audits.

The first bond is anchored +x. A D4 transverse symmetry selects the first
sideways turn +y and first out-of-plane turn +z; orbit weights 1/4/8 recover
EVERY anchored 3D walk, including reflections. A separate full-walk auditor
uses no such reduction and computes contacts from complete vertex paths.

This distinct cubic atlas type cannot be mistaken for an audited 2D atlas.
Exact finite toy-model mathematics; not Lean-refined protein science.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import comb

from .contact_map_atlas import contact_pairs
from .universal_landscape import UniversalLandscape

Point3 = tuple[int, int, int]
_DIRECTIONS: tuple[Point3, ...] = (
    (1, 0, 0), (-1, 0, 0), (0, 1, 0),
    (0, -1, 0), (0, 0, 1), (0, 0, -1),
)
_INDEPENDENT_DIRECTIONS: tuple[Point3, ...] = (
    (0, 0, -1), (0, -1, 0), (-1, 0, 0),
    (0, 0, 1), (0, 1, 0), (1, 0, 0),
)


def _hmask(sequence: str, length: int, pairs: tuple[tuple[int, int], ...]) -> int:
    if (not isinstance(sequence, str) or len(sequence) != length
            or set(sequence) - {"H", "P"}):
        raise ValueError("sequence must have the atlas length and only H/P residues")
    return sum(1 << b for b, (i, j) in enumerate(pairs)
               if sequence[i] == sequence[j] == "H")


@dataclass(frozen=True)
class CubicContactAtlas:
    """Cubic geometry. Structural sanity checks are NOT a completeness proof."""

    length: int
    multiplicities: tuple[tuple[int, int], ...]
    conformations: int

    def __post_init__(self) -> None:
        pairs = contact_pairs(self.length)
        if not isinstance(self.multiplicities, tuple) or not self.multiplicities:
            raise ValueError("a nonempty sorted tuple of contact masks is required")
        last = -1
        total = 0
        for entry in self.multiplicities:
            if not isinstance(entry, tuple) or len(entry) != 2:
                raise ValueError("contact-map entries must be pairs")
            mask, count = entry
            if (type(mask) is not int or type(count) is not int
                    or count <= 0 or not last < mask < (1 << len(pairs))):
                raise ValueError("invalid, duplicate, or unsorted contact map")
            # Exactly one straight walk; all other transverse D4 orbits
            # have size four (planar) or eight (nonplanar).
            if count % 4 != (1 if mask == 0 else 0):
                raise ValueError("contact-map multiplicity violates D4 congruence")
            last = mask
            total += count
        if (self.multiplicities[0][0] != 0 or type(self.conformations) is not int
                or self.conformations != total):
            raise ValueError("missing straight walk or incorrect conformation total")

    @property
    def pairs(self) -> tuple[tuple[int, int], ...]:
        return contact_pairs(self.length)

    def spectrum(self, sequence: str) -> dict[int, int]:
        selected = _hmask(sequence, self.length, self.pairs)
        histogram: Counter[int] = Counter()
        for mask, count in self.multiplicities:
            histogram[(mask & selected).bit_count()] += count
        return dict(sorted(histogram.items()))

    def partition(self, sequence: str, q: Fraction) -> Fraction:
        if type(q) is not Fraction or q < 1:
            raise ValueError("q must be an exact Fraction >=1")
        return sum((v * q**k for k, v in self.spectrum(sequence).items()), Fraction(0))

    def event_probability(self, sequence: str, q: Fraction, threshold: int) -> Fraction:
        if type(q) is not Fraction or q < 1 or type(threshold) is not int:
            raise ValueError("invalid rational weight or event threshold")
        spectrum = self.spectrum(sequence)
        z = sum((v * q**k for k, v in spectrum.items()), Fraction(0))
        event = sum((v * q**k for k, v in spectrum.items() if k >= threshold), Fraction(0))
        return event / z


def build_cubic_atlas(length: int, *, max_nodes: int = 3_000_000) -> CubicContactAtlas:
    """Enumerate transverse D4 orbit representatives; fail on node budget.

    A straight walk has weight 1. First lateral +y stands for 4 possible
    directions. Among these, the first z turn +z stands for its mirror -z,
    giving 8 for each genuinely nonplanar labelled path.
    """
    bits = {pair: 1 << k for k, pair in enumerate(contact_pairs(length))}
    if type(max_nodes) is not int or max_nodes < 1:
        raise ValueError("max_nodes must be a positive integer")
    occupied: dict[Point3, int] = {(0, 0, 0): 0, (1, 0, 0): 1}
    counts: Counter[int] = Counter()
    visited_nodes = 0

    def visit(x: int, y: int, z: int, sideways: bool, out_of_plane: bool, mask: int) -> None:
        nonlocal visited_nodes
        visited_nodes += 1
        if visited_nodes > max_nodes:
            raise ValueError("cubic enumeration exceeded max_nodes; no partial result")
        index = len(occupied)
        if index == length:
            counts[mask] += 8 if out_of_plane else 4 if sideways else 1
            return
        for dx, dy, dz in _DIRECTIONS:
            if (not sideways and (dy < 0 or dz != 0)) or (sideways and not out_of_plane and dz < 0):
                continue
            next_point = (x + dx, y + dy, z + dz)
            if next_point in occupied:
                continue
            new_mask = mask
            for ddx, ddy, ddz in _DIRECTIONS:
                previous = occupied.get((
                    next_point[0] + ddx, next_point[1] + ddy, next_point[2] + ddz
                ))
                if previous is not None and previous < index - 1:
                    new_mask |= bits.get((previous, index), 0)
            occupied[next_point] = index
            visit(*next_point, sideways or dy != 0,
                  out_of_plane or dz != 0, new_mask)
            del occupied[next_point]

    visit(1, 0, 0, False, False, 0)
    return CubicContactAtlas(
        length, tuple(sorted(counts.items())), sum(counts.values())
    )


def audit_cubic_atlas(atlas: CubicContactAtlas, *,
                      max_length: int = 11, max_leaves: int = 2_000_000) -> bool:
    """Full unsymmetrized 3D walk enumeration and independent leaf contact scoring."""
    if not isinstance(atlas, CubicContactAtlas):
        raise TypeError("expected CubicContactAtlas")
    if (type(max_length) is not int or max_length < 2 or atlas.length > max_length
            or type(max_leaves) is not int or max_leaves < 1):
        raise ValueError("independent audit resource limit exceeded")
    path: list[Point3] = [(0, 0, 0), (1, 0, 0)]
    occupied = set(path)
    counts: Counter[int] = Counter()
    pairs = atlas.pairs
    leaves = 0

    def walk() -> None:
        nonlocal leaves
        if len(path) == atlas.length:
            leaves += 1
            if leaves > max_leaves:
                raise ValueError("independent audit leaf budget exceeded")
            mask = 0
            for b, (i, j) in enumerate(pairs):
                if sum(abs(a - c) for a, c in zip(path[i], path[j])) == 1:
                    mask |= 1 << b
            counts[mask] += 1
            return
        x, y, z = path[-1]
        for dx, dy, dz in _INDEPENDENT_DIRECTIONS:
            point = (x + dx, y + dy, z + dz)
            if point in occupied:
                continue
            occupied.add(point)
            path.append(point)
            walk()
            path.pop()
            occupied.remove(point)

    walk()
    return tuple(sorted(counts.items())) == atlas.multiplicities


def build_cubic_landscape(
    atlas: CubicContactAtlas, *, require_audit: bool = True,
    max_length: int = 11, max_edge_subsets: int = 1_000_000,
    max_audit_leaves: int = 2_000_000,
) -> UniversalLandscape:
    """Reuse exact graph-subset/zeta partition algebra with a 3D provenance gate."""
    if not isinstance(atlas, CubicContactAtlas):
        raise TypeError("a distinct CubicContactAtlas is required")
    if type(require_audit) is not bool or type(max_length) is not int or not 2 <= max_length <= 20:
        raise ValueError("invalid audit or length limit")
    if atlas.length > max_length or type(max_edge_subsets) is not int or max_edge_subsets < 1:
        raise ValueError("cubic expansion resource limit exceeded")
    if require_audit and not audit_cubic_atlas(
        atlas, max_length=max_length, max_leaves=max_audit_leaves
    ):
        raise ValueError("cubic geometry failed independent exhaustive audit")

    pairs = atlas.pairs
    endpoints = tuple((1 << i) | (1 << j) for i, j in pairs)
    degree = max(mask.bit_count() for mask, _ in atlas.multiplicities)
    size = 1 << atlas.length
    interactions = [[0] * (degree + 1) for _ in range(size)]
    edges_processed = 0
    for mask, weight in atlas.multiplicities:
        edges = [j for j in range(len(pairs)) if mask & (1 << j)]
        edges_processed += 1 << len(edges)
        if edges_processed > max_edge_subsets:
            raise ValueError("edge-subset enumeration exceeded budget; no partial result")
        states = {(0, 0): 1}
        for edge in edges:
            nxt = dict(states)
            for (vertices, count), multiplicity in states.items():
                key = (vertices | endpoints[edge], count + 1)
                nxt[key] = nxt.get(key, 0) + multiplicity
            states = nxt
        for (vertices, count), multiplicity in states.items():
            interactions[vertices][count] += weight * multiplicity

    partitions = [row.copy() for row in interactions]
    for i in range(atlas.length):
        bit = 1 << i
        for mask in range(size):
            if mask & bit:
                parent = partitions[mask ^ bit]
                child = partitions[mask]
                for j in range(degree + 1):
                    child[j] += parent[j]
    return UniversalLandscape(
        atlas.length, atlas.conformations,
        tuple(tuple(row) for row in interactions),
        tuple(tuple(row) for row in partitions),
    )


def audit_cubic_landscape(
    landscape: UniversalLandscape, atlas: CubicContactAtlas, *,
    max_length: int = 11, max_leaves: int = 2_000_000,
    max_comparisons: int = 4_000_000,
) -> bool:
    """Direct per-map contact scoring checks EVERY sequence polynomial.

    Also independently enumerates all 3D walks. This does NOT use the
    builder's graph-subset expansion or Boolean zeta transform.
    """
    if not isinstance(landscape, UniversalLandscape) or not isinstance(atlas, CubicContactAtlas):
        raise TypeError("expected polynomial table and cubic atlas")
    if (type(max_comparisons) is not int or max_comparisons < 1
            or (1 << atlas.length) * len(atlas.multiplicities) > max_comparisons):
        raise ValueError("complete polynomial audit comparison budget exceeded")
    if landscape.length != atlas.length or landscape.conformations != atlas.conformations:
        return False
    if not audit_cubic_atlas(atlas, max_length=max_length, max_leaves=max_leaves):
        return False
    degree = len(landscape.partitions[0])
    for residues in range(1 << atlas.length):
        contact_mask = sum(1 << b for b, (i, j) in enumerate(atlas.pairs)
                           if residues & (1 << i) and residues & (1 << j))
        expected = [0] * degree
        for mask, weight in atlas.multiplicities:
            k = (mask & contact_mask).bit_count()
            for j in range(k + 1):
                expected[j] += weight * comb(k, j)
        if tuple(expected) != landscape.partitions[residues]:
            return False
    return True
