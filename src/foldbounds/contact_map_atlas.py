"""Reusable exact contact-map density of states for the anchored 2D HP model.

Geometry is enumerated *once*, independently of the H/P sequence.  A lattice
reflection halves the search, but its multiplicity is restored exactly.

This is an exponential reference oracle, not a Lean-certified implementation
or a claim about physical proteins.  Serialized atlases are untrusted until
independently re-enumerated with :func:`audit_atlas`.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

Point = tuple[int, int]
DIRECTIONS: tuple[Point, ...] = ((1, 0), (0, 1), (-1, 0), (0, -1))


def contact_pairs(length: int) -> tuple[tuple[int, int], ...]:
    """Potential nonconsecutive neighbor pairs (checkerboard parity exact)."""
    if type(length) is not int or length < 2:
        raise ValueError("length must be an integer >= 2")
    return tuple((i, j) for j in range(length) for i in range(j - 2)
                 if (j - i) % 2 == 1)


def _sequence_mask(sequence: str, length: int, pairs: tuple[tuple[int, int], ...]) -> int:
    if not isinstance(sequence, str) or len(sequence) != length or set(sequence) - {"H", "P"}:
        raise ValueError("sequence must be H/P and have atlas length")
    return sum(1 << bit for bit, (i, j) in enumerate(pairs)
               if sequence[i] == "H" and sequence[j] == "H")


def _q(q: Fraction) -> Fraction:
    if type(q) is not Fraction or q < 1:
        raise ValueError("q must be an exact Fraction >= 1")
    return q


@dataclass(frozen=True)
class ContactMapAtlas:
    """A geometry-only histogram of contact bitmasks with exact multiplicities.

    Structural consistency checks do *not* establish that a claimed histogram
    enumerates every walk.  Use ``audit_atlas`` for a second full enumeration.
    """

    length: int
    multiplicities: tuple[tuple[int, int], ...]
    conformations: int

    def __post_init__(self) -> None:
        pairs = contact_pairs(self.length)
        if not isinstance(self.multiplicities, tuple) or not self.multiplicities:
            raise ValueError("multiplicities must be a nonempty tuple")
        last = -1
        total = 0
        for entry in self.multiplicities:
            if not isinstance(entry, tuple) or len(entry) != 2:
                raise ValueError("each entry must be (mask, positive multiplicity)")
            mask, count = entry
            if (type(mask) is not int or type(count) is not int
                    or not last < mask < (1 << len(pairs)) or count < 1):
                raise ValueError("invalid, duplicate, or unsorted contact mask")
            # Exactly one fixed-by-reflection walk (the straight one) exists.
            if (mask == 0 and count % 2 != 1) or (mask != 0 and count % 2 != 0):
                raise ValueError("reflection multiplicities are inconsistent")
            last = mask
            total += count
        if self.multiplicities[0][0] != 0 or type(self.conformations) is not int or total != self.conformations:
            raise ValueError("missing straight walk or incorrect total")

    @property
    def pairs(self) -> tuple[tuple[int, int], ...]:
        return contact_pairs(self.length)

    def spectrum(self, sequence: str) -> dict[int, int]:
        """Exact score histogram *assuming the atlas was correctly built*."""
        hydrophobic = _sequence_mask(sequence, self.length, self.pairs)
        out: Counter[int] = Counter()
        for mask, count in self.multiplicities:
            out[(mask & hydrophobic).bit_count()] += count
        return dict(sorted(out.items()))

    def partition(self, sequence: str, q: Fraction) -> Fraction:
        """Exact Z(q) in the finite anchored lattice model."""
        q = _q(q)
        return sum((count * q ** k for k, count in self.spectrum(sequence).items()), Fraction(0))

    def event_probability(self, sequence: str, q: Fraction, threshold: int) -> Fraction:
        """Exact probability of at least ``threshold`` H-H contacts."""
        q = _q(q)
        if type(threshold) is not int:
            raise ValueError("threshold must be an integer")
        spectrum = self.spectrum(sequence)
        partition = sum((count * q ** k for k, count in spectrum.items()), Fraction(0))
        numerator = sum((count * q ** k for k, count in spectrum.items() if k >= threshold), Fraction(0))
        return numerator / partition

    def contact_moments(self, sequence: str, q: Fraction) -> tuple[Fraction, Fraction]:
        """Return exact mean and variance of K under the Boltzmann ensemble.

        In particular, d E[K]/d(log q) = Var(K) >= 0; this routine evaluates
        the finite sums, not approximate derivatives or physical temperatures.
        """
        q = _q(q)
        spectrum = self.spectrum(sequence)
        z = sum((count * q ** k for k, count in spectrum.items()), Fraction(0))
        first = sum((k * count * q ** k for k, count in spectrum.items()), Fraction(0)) / z
        second = sum((k * k * count * q ** k for k, count in spectrum.items()), Fraction(0)) / z
        return first, second - first * first

    def score_covariance(self, reference: str, comparison: str, q: Fraction) -> Fraction:
        """Cov(K_ref, K_comp) under the reference sequence's q-weighted ensemble.

        Both score functions are evaluated on *the same geometries*. This is
        not a cross-protein covariance or a model of mutation kinetics.
        """
        q = _q(q)
        first = _sequence_mask(reference, self.length, self.pairs)
        second = _sequence_mask(comparison, self.length, self.pairs)
        z = ra = rb = rab = 0
        for contact_map, count in self.multiplicities:
            a = (contact_map & first).bit_count()
            b = (contact_map & second).bit_count()
            w = count * q ** a
            z += w
            ra += a * w
            rb += b * w
            rab += a * b * w
        return Fraction(rab, z) - Fraction(ra, z) * Fraction(rb, z)


def build_atlas(length: int, *, max_nodes: int = 2_000_000) -> ContactMapAtlas:
    """Enumerate canonical SAWs with reflection pairing and an explicit budget.

    Every nonstraight anchored walk has a unique mirror with the opposite
    *first vertical turn*. Only the positive-turn representative is explored;
    its contact mask has multiplicity two. The straight walk has weight one.
    Exceeding the search budget raises, never returns a partial atlas.
    """
    pairs = contact_pairs(length)
    if type(max_nodes) is not int or max_nodes < 1:
        raise ValueError("max_nodes must be a positive integer")
    bit_for_pair = {pair: 1 << bit for bit, pair in enumerate(pairs)}
    occupied: dict[Point, int] = {(0, 0): 0, (1, 0): 1}
    counts: Counter[int] = Counter()
    visited_nodes = 0

    def visit(x: int, y: int, turned: bool, mask: int) -> None:
        nonlocal visited_nodes
        visited_nodes += 1
        if visited_nodes > max_nodes:
            raise ValueError("atlas enumeration exceeded max_nodes; no partial result")
        index = len(occupied)
        if index == length:
            counts[mask] += 2 if turned else 1
            return
        for dx, dy in DIRECTIONS:
            candidate = (x + dx, y + dy)
            if candidate in occupied or (not turned and dy < 0):
                continue
            new_mask = mask
            for nx, ny in ((candidate[0] + ddx, candidate[1] + ddy) for ddx, ddy in DIRECTIONS):
                previous = occupied.get((nx, ny))
                if previous is not None and previous < index - 1:
                    new_mask |= bit_for_pair.get((previous, index), 0)
            occupied[candidate] = index
            visit(candidate[0], candidate[1], turned or dy != 0, new_mask)
            del occupied[candidate]

    visit(1, 0, False, 0)
    histogram = tuple(sorted(counts.items()))
    return ContactMapAtlas(length, histogram, sum(counts.values()))


def audit_atlas(atlas: ContactMapAtlas, *, max_length: int = 12) -> bool:
    """Independently recompute *all* walks without reflection reduction.

    This deliberately uses a separate full-walk enumeration and leaf-pair
    scorer, not ``build_atlas`` or its incremental mask update. It still runs
    in ordinary Python and is not a formal proof. Refuses large audits.
    """
    if not isinstance(atlas, ContactMapAtlas):
        raise TypeError("expected ContactMapAtlas")
    if type(max_length) is not int or max_length < 2 or atlas.length > max_length:
        raise ValueError("refusing full independent audit beyond max_length")
    pairs = atlas.pairs
    seen = {(0, 0), (1, 0)}
    walk: list[Point] = [(0, 0), (1, 0)]
    independent: Counter[int] = Counter()

    def visit() -> None:
        if len(walk) == atlas.length:
            bits = 0
            for bit, (i, j) in enumerate(pairs):
                if abs(walk[i][0] - walk[j][0]) + abs(walk[i][1] - walk[j][1]) == 1:
                    bits |= 1 << bit
            independent[bits] += 1
            return
        x, y = walk[-1]
        for dx, dy in ((0, -1), (-1, 0), (0, 1), (1, 0)):
            nxt = (x + dx, y + dy)
            if nxt not in seen:
                seen.add(nxt)
                walk.append(nxt)
                visit()
                walk.pop()
                seen.remove(nxt)

    visit()
    return tuple(sorted(independent.items())) == atlas.multiplicities
