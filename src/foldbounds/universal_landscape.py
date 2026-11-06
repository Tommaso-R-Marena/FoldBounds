"""Exact nonnegative all-sequence HP partition polynomials from a contact atlas.

For a conformation contact graph G and hydrophobic indicator x_i in {0,1},
q**K_G(x) = product_{(i,j) in G} (1 + (q-1)*x_i*x_j).
Expanding edge subsets yields nonnegative integer coefficients c[T,k] such that
Z_S(q) = sum_{T subset S} sum_k c[T,k]*(q-1)**k.

A Boolean subset zeta transform computes all 2**n sequence polynomials at once.
This is exact finite-model combinatorics, not a Lean proof or a real-protein model.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb

from .contact_map_atlas import ContactMapAtlas, audit_atlas, contact_pairs


def _sequence_to_bits(sequence: str, length: int) -> int:
    if not isinstance(sequence, str) or len(sequence) != length or set(sequence) - {"H", "P"}:
        raise ValueError("sequence must contain exactly atlas-length H/P residues")
    return sum((1 << i) for i, r in enumerate(sequence) if r == "H")


def _rational_q(q: Fraction) -> Fraction:
    if type(q) is not Fraction or q < 1:
        raise ValueError("q must be an exact Fraction >= 1")
    return q


def _evaluate(coefficients: tuple[int, ...], t: Fraction) -> Fraction:
    """Evaluate integer polynomial via Horner, without floating point."""
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * t + coefficient
    return result


def _trim(poly: tuple[int, ...]) -> tuple[int, ...]:
    last = len(poly)
    while last > 1 and poly[last - 1] == 0:
        last -= 1
    return poly[:last]


@dataclass(frozen=True)
class UniversalLandscape:
    """Immutable exact polynomial table indexed by the H-residue bitmask.

    interactions[T] contains the coefficients c[T,k] of (q-1)^k.
    partitions[S] contains the coefficients of Z_S(q), and equals the
    subset-zeta transform of interactions. Both arrays are dense and are
    therefore intended for modest n only.

    This type is not proof-carrying if instantiated with arbitrary data.
    Prefer build_landscape(require_audit=True) from an independently
    cross-checked atlas and audit_landscape for independent polynomial checks.
    """

    length: int
    conformations: int
    interactions: tuple[tuple[int, ...], ...]
    partitions: tuple[tuple[int, ...], ...]

    def __post_init__(self) -> None:
        if type(self.length) is not int or not 2 <= self.length <= 20:
            raise ValueError("invalid landscape length")
        size = 1 << self.length
        if (type(self.conformations) is not int or self.conformations < 1
                or len(self.interactions) != size or len(self.partitions) != size):
            raise ValueError("incorrect landscape dimensions")
        degree = len(self.interactions[0])
        if not 1 <= degree <= len(contact_pairs(self.length)) + 1:
            raise ValueError("invalid polynomial degree")
        for family in (self.interactions, self.partitions):
            for row in family:
                if (not isinstance(row, tuple) or len(row) != degree
                        or any(type(x) is not int or x < 0 for x in row)):
                    raise ValueError("polynomial coefficients must be nonnegative integers")
        if self.interactions[0][0] != self.conformations:
            raise ValueError("the empty edge-subgraph coefficient must count all walks")
        if any(self.partitions[s][0] != self.conformations for s in range(size)):
            raise ValueError("all sequences must have Z(q=1) = number of walks")

    def polynomial(self, sequence: str) -> tuple[int, ...]:
        """Integer coefficients of Z(q) in t=q-1."""
        return _trim(self.partitions[_sequence_to_bits(sequence, self.length)])

    def partition(self, sequence: str, q: Fraction) -> Fraction:
        """Exact partition function at rational q >= 1."""
        return _evaluate(self.polynomial(sequence), _rational_q(q) - 1)

    def all_partitions(self, q: Fraction) -> tuple[Fraction, ...]:
        """All 2**n exact partition functions in bitmask index order."""
        t = _rational_q(q) - 1
        return tuple(_evaluate(row, t) for row in self.partitions)

    def interaction(self, residues: str, q: Fraction) -> Fraction:
        """The nonnegative Moebius coefficient c_T(q), with T = H positions."""
        mask = _sequence_to_bits(residues, self.length)
        return _evaluate(self.interactions[mask], _rational_q(q) - 1)

    def mixed_difference(self, baseline: str, additions: tuple[int, ...], q: Fraction) -> Fraction:
        """Exact mixed addition difference of Z, proven >= 0 for q >= 1.

        Residue indices are zero-based, distinct and initially P in baseline.
        This computes a positive sum of interaction coefficients, not a
        potentially ill-conditioned alternating finite difference.
        """
        t = _rational_q(q) - 1
        base = _sequence_to_bits(baseline, self.length)
        if not isinstance(additions, tuple) or any(type(i) is not int or not 0 <= i < self.length
                                                   for i in additions):
            raise ValueError("additions must be a tuple of valid zero-based indices")
        add = sum(1 << i for i in additions)
        if add.bit_count() != len(additions) or add & base:
            raise ValueError("addition residues must be distinct P positions")
        result = Fraction(0)
        subset = base
        while True:
            result += _evaluate(self.interactions[add | subset], t)
            if subset == 0:
                break
            subset = (subset - 1) & base
        return result

    def contact_spectrum(self, sequence: str) -> dict[int, int]:
        """Recover exact contact-count density from the t-basis polynomial.

        If Z(q) = sum_k N_k q^k = sum_j a_j (q-1)^j, then
        N_k = sum_{j>=k} (-1)^(j-k) binom(j,k) a_j.
        """
        a = self.partitions[_sequence_to_bits(sequence, self.length)]
        counts = [sum(((-1) ** (j - k)) * comb(j, k) * a[j]
                      for j in range(k, len(a))) for k in range(len(a))]
        if any(v < 0 for v in counts) or sum(counts) != self.conformations:
            raise AssertionError("internal spectrum identity failed")
        return {k: value for k, value in enumerate(counts) if value}

    def contact_moments(self, sequence: str, q: Fraction) -> tuple[Fraction, Fraction]:
        """Exact (mean, variance) of contact score via partition derivatives."""
        q = _rational_q(q)
        poly = self.polynomial(sequence)
        t = q - 1
        z = _evaluate(poly, t)
        first = _evaluate(tuple(k * poly[k] for k in range(1, len(poly))), t)
        second = _evaluate(tuple(k * (k-1) * poly[k] for k in range(2, len(poly))), t)
        mean = q * first / z
        variance = mean + (q * q * second / z) - mean * mean
        if variance < 0:
            raise AssertionError("nonnegative variance identity violated")
        return mean, variance


def build_landscape(atlas: ContactMapAtlas, *, require_audit: bool = True,
                    max_length: int = 16, max_edge_subsets: int = 1_000_000) -> UniversalLandscape:
    """Build the universal nonnegative partition polynomial by edge expansion.

    max_edge_subsets bounds the total edge-subgraph subsets processed across
    contact maps. A budget failure raises rather than returning partial data.
    require_audit=True independently re-enumerates the entire atlas and
    refuses lengths above the atlas audit limit (12).
    """
    if not isinstance(atlas, ContactMapAtlas):
        raise TypeError("expected a ContactMapAtlas")
    if type(require_audit) is not bool:
        raise ValueError("require_audit must be bool")
    if type(max_length) is not int or not 2 <= max_length <= 20 or atlas.length > max_length:
        raise ValueError("landscape exceeds requested max_length")
    if type(max_edge_subsets) is not int or max_edge_subsets < 1:
        raise ValueError("max_edge_subsets must be positive int")
    if require_audit and not audit_atlas(atlas):
        raise ValueError("contact-map atlas failed independent exhaustive audit")
    pairs = atlas.pairs
    edge_vertices = tuple((1 << i) | (1 << j) for i, j in pairs)
    maximum_degree = max(mask.bit_count() for mask, _ in atlas.multiplicities)
    size = 1 << atlas.length
    interactions = [[0] * (maximum_degree + 1) for _ in range(size)]
    total_edge_subsets = 0
    for mask, multiplicity in atlas.multiplicities:
        edge_indices = [j for j in range(len(pairs)) if mask & (1 << j)]
        total_edge_subsets += 1 << len(edge_indices)
        if total_edge_subsets > max_edge_subsets:
            raise ValueError("contact-subgraph expansion exceeded budget; no partial landscape")
        states: dict[tuple[int, int], int] = {(0, 0): 1}
        for index in edge_indices:
            endpoints = edge_vertices[index]
            new = states.copy()
            for (vertices, count), ways in states.items():
                key = (vertices | endpoints, count + 1)
                new[key] = new.get(key, 0) + ways
            states = new
        for (vertices, count), ways in states.items():
            interactions[vertices][count] += multiplicity * ways

    partitions = [row.copy() for row in interactions]
    for i in range(atlas.length):
        bit = 1 << i
        for subset in range(size):
            if subset & bit:
                parent = partitions[subset ^ bit]
                child = partitions[subset]
                for degree in range(maximum_degree + 1):
                    child[degree] += parent[degree]
    return UniversalLandscape(atlas.length, atlas.conformations,
                              tuple(tuple(row) for row in interactions),
                              tuple(tuple(row) for row in partitions))


def audit_landscape(landscape: UniversalLandscape, atlas: ContactMapAtlas, *,
                    max_length: int = 12, max_mask_checks: int = 6_000_000) -> bool:
    """Independently check every sequence partition polynomial.

    For each HP sequence, score every contact map directly and reconstruct
    coefficients of (q-1)^k as sum_G m_G * binom(K_G,k).
    Then independently invert the subset-zeta transform to check the stored
    interaction coefficients. The atlas itself is checked by separate full
    SAW enumeration. Resource limits fail closed. This is a Python audit,
    not a Lean theorem or cryptographic authentication.
    """
    if not isinstance(landscape, UniversalLandscape) or not isinstance(atlas, ContactMapAtlas):
        raise TypeError("expected UniversalLandscape and ContactMapAtlas")
    if type(max_length) is not int or not 2 <= max_length <= 20:
        raise ValueError("invalid independent audit length limit")
    if type(max_mask_checks) is not int or max_mask_checks < 1:
        raise ValueError("invalid mask-check budget")
    if atlas.length > max_length or (1 << atlas.length) * len(atlas.multiplicities) > max_mask_checks:
        raise ValueError("full independent landscape audit exceeds resource budget")
    if landscape.length != atlas.length or landscape.conformations != atlas.conformations:
        return False
    if not audit_atlas(atlas, max_length=max_length):
        return False
    degree = len(landscape.partitions[0])
    pair_residue_bits = tuple((1 << i) | (1 << j) for i, j in atlas.pairs)
    pair_masks = []
    for sequence_mask in range(1 << atlas.length):
        mask = sum((1 << i) for i, endpoints in enumerate(pair_residue_bits)
                   if sequence_mask & endpoints == endpoints)
        pair_masks.append(mask)
    for sequence_mask, pair_mask in enumerate(pair_masks):
        expected = [0] * degree
        for graph, multiplicity in atlas.multiplicities:
            k = (graph & pair_mask).bit_count()
            if k >= degree:
                return False
            for j in range(k + 1):
                expected[j] += multiplicity * comb(k, j)
        if tuple(expected) != landscape.partitions[sequence_mask]:
            return False
    inverse = [list(row) for row in landscape.partitions]
    for i in range(atlas.length):
        bit = 1 << i
        for mask in range(1 << atlas.length):
            if mask & bit:
                source, destination = inverse[mask ^ bit], inverse[mask]
                for k in range(degree):
                    destination[k] -= source[k]
    return tuple(tuple(row) for row in inverse) == landscape.interactions
