"""A small, independently checked 2D square-lattice HP model.

This module implements *model-relative* mathematical guarantees; it does not
certify real protein physics, machine floating point, or new algorithmic novelty.

Canonical conformations are square-lattice self-avoiding walks starting at
(0, 0) with the first bond to (1, 0). Reflection is NOT quotiented out.
E = -epsilon * (number of nonconsecutive H-H lattice-neighbor pairs).
Set q = exp(beta*epsilon) >= 1, supplied as an exact rational parameter.
The Boltzmann weight is therefore q**contacts; no approximate exponentials.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

Point = tuple[int, int]
Path = tuple[Point, ...]
Direction = tuple[int, int]
DIRECTIONS: tuple[Direction, ...] = ((1, 0), (0, 1), (-1, 0), (0, -1))
ROOT: Path = ((0, 0), (1, 0))


class CertificateError(ValueError):
    """Malformed, invalid, or incomplete conformational frontier."""


def _check_sequence(sequence: str) -> None:
    if not isinstance(sequence, str) or len(sequence) < 2 or set(sequence) - {"H", "P"}:
        raise ValueError("sequence must contain at least two H/P residues")


def _check_q(q: Fraction) -> Fraction:
    if not isinstance(q, Fraction) or q < 1:
        raise ValueError("q must be an exact Fraction >= 1")
    return q


def _neighbors(path: Path) -> tuple[Path, ...]:
    occupied = set(path)
    x, y = path[-1]
    return tuple(path + ((x + dx, y + dy),) for dx, dy in DIRECTIONS
                 if (x + dx, y + dy) not in occupied)


def _validate_path(path: Path, length: int) -> None:
    if not isinstance(path, tuple) or len(path) < 2 or len(path) > length or path[:2] != ROOT:
        raise CertificateError("frontier node must extend the canonical two-site root")
    if any(not isinstance(p, tuple) or len(p) != 2 or
           any(type(v) is not int for v in p) for p in path):
        raise CertificateError("coordinates must be pairs of integers")
    if len(set(path)) != len(path):
        raise CertificateError("self-intersecting path")
    if any(abs(x1-x0) + abs(y1-y0) != 1 for (x0, y0), (x1, y1) in zip(path, path[1:])):
        raise CertificateError("non-adjacent lattice steps")


def _possible_hh_pairs(sequence: str) -> tuple[tuple[int, int], ...]:
    # Two vertices separated by an even chain distance have the same checkerboard
    # parity and cannot be adjacent on a square lattice.
    return tuple((i, j) for i in range(len(sequence)) for j in range(i + 3, len(sequence))
                 if (j - i) % 2 and sequence[i] == "H" and sequence[j] == "H")


def _contacts(path: Path, pairs: tuple[tuple[int, int], ...]) -> int:
    n = len(path)
    return sum(abs(path[i][0]-path[j][0]) + abs(path[i][1]-path[j][1]) == 1
               for i, j in pairs if j < n)


def _score_range(path: Path, pairs: tuple[tuple[int, int], ...]) -> tuple[int, int]:
    known = _contacts(path, pairs)
    remaining_potential = sum(j >= len(path) for _, j in pairs)
    return known, known + remaining_potential


def exact_spectrum(sequence: str) -> dict[int, int]:
    """Return the exact contact-count density of states (exponential runtime)."""
    _check_sequence(sequence)
    pairs = _possible_hh_pairs(sequence)
    spectrum: Counter[int] = Counter()

    def visit(path: Path) -> None:
        if len(path) == len(sequence):
            spectrum[_contacts(path, pairs)] += 1
            return
        for child in _neighbors(path):
            visit(child)

    visit(ROOT)
    return dict(sorted(spectrum.items()))


@dataclass(frozen=True)
class Node:
    path: Path
    kind: Literal["exact", "bound"]


@dataclass(frozen=True)
class Certificate:
    sequence: str
    q: Fraction
    event_at_least: int
    frontier: tuple[Node, ...]


@dataclass(frozen=True)
class VerifiedBounds:
    partition_lower: Fraction
    partition_upper: Fraction
    event_lower: Fraction
    event_upper: Fraction
    complement_lower: Fraction
    complement_upper: Fraction
    probability_lower: Fraction
    probability_upper: Fraction
    exact_paths: int
    bounded_subtrees: int

    @property
    def exact(self) -> bool:
        return self.partition_lower == self.partition_upper


def build_frontier(sequence: str, q: Fraction, event_at_least: int,
                   depth: int) -> Certificate:
    """Construct a valid prefix frontier, exact at terminal nodes.

    Depth is the number of sequence residues fixed in each frontier prefix.
    """
    _check_sequence(sequence)
    _check_q(q)
    if type(event_at_least) is not int or type(depth) is not int or not 2 <= depth <= len(sequence):
        raise ValueError("integer event threshold and depth within [2, sequence length] required")
    nodes: list[Node] = []

    def visit(path: Path) -> None:
        if len(path) == len(sequence):
            nodes.append(Node(path, "exact"))
        elif len(path) == depth:
            nodes.append(Node(path, "bound"))
        else:
            children = _neighbors(path)
            for child in children:
                visit(child)

    visit(ROOT)
    return Certificate(sequence, q, event_at_least, tuple(nodes))


def refine_frontier(cert: Certificate, path: Path) -> Certificate:
    """Replace one bounded node with all its legal immediate successors.

    This is a deterministic refinement operation; verification is independent.
    """
    _check_sequence(cert.sequence)
    indices = [i for i, n in enumerate(cert.frontier) if n.path == path and n.kind == "bound"]
    if len(indices) != 1:
        raise CertificateError("refinement requires exactly one bounded frontier node")
    i = indices[0]
    if len(path) >= len(cert.sequence):
        raise CertificateError("cannot refine completed path")
    successors = tuple(Node(child, "exact" if len(child) == len(cert.sequence) else "bound")
                       for child in _neighbors(path))
    return Certificate(cert.sequence, cert.q, cert.event_at_least,
                       cert.frontier[:i] + successors + cert.frontier[i+1:])


def verify_certificate(cert: Certificate) -> VerifiedBounds:
    """Check full canonical SAW-tree coverage and return rigorous rational bounds.

    IMPORTANT: the verifier does not trust frontier completeness or node validity.
    Each subtree upper bound is 3**remaining * q**max_possible_contacts:
    a SAW has at most 3 choices after each step, and all future contacts must
    belong to nonconsecutive parity-compatible H-H pairs. Missing branches,
    duplicate nodes, invalid paths, and prefix-overlapping nodes are rejected.
    """
    _check_sequence(cert.sequence)
    q = _check_q(cert.q)
    if type(cert.event_at_least) is not int:
        raise CertificateError("event_at_least must be an integer")
    length = len(cert.sequence)
    pairs = _possible_hh_pairs(cert.sequence)
    nodes: dict[Path, str] = {}
    for node in cert.frontier:
        if node.kind not in ("exact", "bound"):
            raise CertificateError("unknown node kind")
        _validate_path(node.path, length)
        if node.kind == "exact" and len(node.path) != length:
            raise CertificateError("exact node must be terminal")
        if node.path in nodes:
            raise CertificateError("duplicate frontier path")
        nodes[node.path] = node.kind

    # For efficient validation of frontier paths: follow only prefixes that
    # lead to submitted nodes, rejecting hidden extras after traversal.
    known_prefixes: set[Path] = set()
    for path in nodes:
        for k in range(2, len(path) + 1):
            known_prefixes.add(path[:k])

    zlo = zhi = elo = ehi = nlo = nhi = Fraction(0)
    exact_count = bounded_count = 0
    encountered: set[Path] = set()

    def visit(path: Path) -> None:
        nonlocal zlo, zhi, elo, ehi, nlo, nhi, exact_count, bounded_count
        if path in nodes:
            encountered.add(path)
            kind = nodes[path]
            if kind == "exact" or len(path) == length:
                contacts = _contacts(path, pairs)
                weight = q ** contacts
                zlo += weight
                zhi += weight
                if contacts >= cert.event_at_least:
                    elo += weight
                    ehi += weight
                else:
                    nlo += weight
                    nhi += weight
                exact_count += 1
                return
            minimum, maximum = _score_range(path, pairs)
            count_upper = 3 ** (length - len(path))
            upper = count_upper * q ** maximum
            # Conservatively omit any lower mass contribution for this subtree.
            zhi += upper
            if minimum >= cert.event_at_least:
                ehi += upper
            elif maximum < cert.event_at_least:
                nhi += upper
            else:
                ehi += upper
                nhi += upper
            bounded_count += 1
            return
        if len(path) == length:
            raise CertificateError("certificate omits a complete conformation")
        for child in _neighbors(path):
            if child not in known_prefixes:
                raise CertificateError("certificate omits a nonempty subtree")
            visit(child)

    if ROOT not in known_prefixes:
        raise CertificateError("certificate omits the root")
    visit(ROOT)
    if encountered != set(nodes):
        raise CertificateError("frontier has nodes hidden by a covering ancestor")
    # Every sequence has a straight-line conformation, so true Z > 0.
    # Unknown-subtree lower bounds are zero. Treat limiting ratios safely.
    p_lo = elo / (elo + nhi) if elo else Fraction(0)
    p_hi = Fraction(0) if not ehi else (ehi / (ehi + nlo) if nlo else Fraction(1))
    if not (zlo <= zhi and elo <= ehi and nlo <= nhi and 0 <= p_lo <= p_hi <= 1):
        raise AssertionError("internal bounds inconsistency")
    return VerifiedBounds(zlo, zhi, elo, ehi, nlo, nhi, p_lo, p_hi,
                          exact_count, bounded_count)
