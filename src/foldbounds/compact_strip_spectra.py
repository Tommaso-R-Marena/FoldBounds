"""Exact HP contact spectra of a guaranteed, compact, exponential SAW subset.

Construct rooted walks by concatenating h-by-4 Hamiltonian *tiles*. Every
full-width tile uses all 4h grid sites exactly once, starts at the upper-left
corner and finishes at the upper-right corner; the first tile must begin East
so that the total walk uses FoldBounds' canonical anchored first edge.
Consecutive tiles are connected by one East backbone edge. Optional leftover
residues extend as a straight East ray. This is a subset of the infinite 2D
rooted square-lattice SAW ensemble, not an approximation to its geometry.

The intertile contact score depends only on H occupancy of the two exposed
boundary columns. A small transfer dynamic program computes an *exact integer
contact-score histogram for the restricted subset*, without generating the
exponentially many concatenated paths. All statements are relative to the HP
toy model and Python; no Lean or real-protein claim is implied.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

Point = tuple[int, int]
Path = tuple[Point, ...]
STEPS = ((1, 0), (0, 1), (-1, 0), (0, -1))


def _check_sequence(sequence: str) -> None:
    if not isinstance(sequence, str) or len(sequence) < 2 or set(sequence) - {"H", "P"}:
        raise ValueError("expected H/P sequence with at least two residues")


@lru_cache(maxsize=8)
def tile_paths(height: int = 8) -> tuple[Path, ...]:
    """All finite Hamiltonian height-by-4 paths from (0,0) to (w-1,0).

    This is intentionally a small, explicit, independently recheckable witness
    generator. No unverified asymptotic Hamiltonian-path count is assumed.
    """
    if type(height) is not int or not 3 <= height <= 8:
        raise ValueError("tile height must be an integer in [3,8]")
    width = 4
    end = (3, 0)
    length = height * width
    used = {(0, 0)}
    result: list[Path] = []
    walk: list[Point] = [(0, 0)]

    def visit() -> None:
        if len(walk) == length:
            if walk[-1] == end:
                result.append(tuple(walk))
            return
        x, y = walk[-1]
        for dx, dy in STEPS:
            p = (x + dx, y + dy)
            if not (0 <= p[0] < width and 0 <= p[1] < height):
                continue
            if p in used or (p == end and len(walk) + 1 < length):
                continue
            used.add(p)
            walk.append(p)
            visit()
            walk.pop()
            used.remove(p)

    visit()
    return tuple(result)


def validate_tile_library(height: int, paths: tuple[Path, ...]) -> bool:
    """Check each supplied tile rather than trust Hamiltonian-path generator."""
    if type(height) is not int or not 3 <= height <= 8 or not isinstance(paths, tuple):
        raise ValueError("invalid tile height or path library")
    width = 4
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("tile library empty or has duplicates")
    cells = {(x, y) for x in range(width) for y in range(height)}
    for path in paths:
        if (not isinstance(path, tuple) or len(path) != len(cells) or
                path[0] != (0, 0) or path[-1] != (width - 1, 0) or
                any(not isinstance(p, tuple) or len(p) != 2 or
                    any(type(v) is not int for v in p) for p in path)):
            raise ValueError("malformed tile")
        if set(path) != cells or len(set(path)) != len(cells):
            raise ValueError("tile does not visit each site exactly once")
        if any(abs(a[0]-b[0]) + abs(a[1]-b[1]) != 1
               for a, b in zip(path, path[1:])):
            raise ValueError("tile has nonadjacent step")
    return True


def _details(sequence: str, block_index: int, height: int,
             paths: tuple[Path, ...]) -> tuple[tuple[int, int, int], ...]:
    """For each tile: (within-tile HH score, exposed left H mask, right H mask)."""
    width=4
    left_offset = block_index * (height * width)
    out = []
    for path in paths:
        indices = {p: i for i, p in enumerate(path)}
        block = sequence[left_offset:left_offset + height * width]
        hsites = {p for p in path if block[indices[p]] == "H"}
        contact = 0
        for x, y in hsites:
            for v in ((x + 1, y), (x, y + 1)):
                if v in hsites and abs(indices[(x, y)] - indices[v]) != 1:
                    contact += 1
        first = sum(1 << y for y in range(height) if (0, y) in hsites)
        last = sum(1 << y for y in range(height) if (width-1, y) in hsites)
        out.append((contact, first, last))
    return tuple(out)


def _direction_bits(path: Path) -> int:
    """W=1, N=2, S=4; the seam and straight tail introduce E only."""
    out = 0
    for a, b in zip(path, path[1:]):
        dx, dy = b[0] - a[0], b[1] - a[1]
        if dx == -1:
            out |= 1
        elif dy == 1:
            out |= 2
        elif dy == -1:
            out |= 4
    return out


def _spectrum_from_library(sequence: str, height: int,
                           paths: tuple[Path, ...], *,
                           disjoint_from_partially_directed: bool = False
                           ) -> tuple[tuple[int, int], ...]:
    blocks = len(sequence) // (height * 4)
    if blocks == 0:
        return () if disjoint_from_partially_directed else ((0, 1),)
    if not disjoint_from_partially_directed and len(set(sequence)) == 1 and sequence[0] == "H":
        # All occupied sites in each tile are H. The entire completed set
        # forms a filled rectangle and the straight tail adds no contacts.
        # Perimeter identity gives (height-1)*(4*blocks-1) contacts exactly.
        total = sum(path[1] == (1, 0) for path in paths)
        count = total * len(paths) ** (blocks - 1)
        return (((height - 1) * (4 * blocks - 1), count),)
    if not disjoint_from_partially_directed and len(set(sequence)) == 1 and sequence[0] == "P":
        total = sum(path[1] == (1, 0) for path in paths)
        return ((0, total * len(paths) ** (blocks - 1)),)
    states: dict[tuple[int, int, int], int] = {}
    directions = tuple(_direction_bits(path) for path in paths)
    for block_index in range(blocks):
        nxt: dict[tuple[int, int, int], int] = defaultdict(int)
        # The score and exposed H masks depend on the local H/P substring,
        # not on which other tile variants were chosen.
        for j, (within, leftmask, rightmask) in enumerate(
                _details(sequence, block_index, height, paths)):
            bits = directions[j]
            if block_index == 0:
                if paths[j][1] != (1, 0):
                    continue
                nxt[(rightmask, bits, within)] += 1
            else:
                for (prev_right, prev_bits, score), count in states.items():
                    # There is exactly one seam backbone edge at y=0.
                    seam = ((prev_right & leftmask) >> 1).bit_count()
                    nxt[(rightmask, prev_bits | bits, score + within + seam)] += count
        states = nxt
    totals: dict[int, int] = defaultdict(int)
    for (_, bits, score), count in states.items():
        if disjoint_from_partially_directed and bits != 7:
            continue
        totals[score] += count
    return tuple(sorted(totals.items()))


@lru_cache(maxsize=128)
def compact_strip_spectrum(sequence: str, *, height: int = 8,
                           disjoint_from_partially_directed: bool = False,
                           include_reflection: bool = False
                           ) -> tuple[tuple[int, int], ...]:
    """Exact count-by-contact histogram of a *genuine* rooted HP SAW subset.

    For n not divisible by 4h the remaining vertices form a straight East
    tail outside the rectangle. The first tail step is a backbone edge, with
    no additional nonconsecutive contact to any earlier tile vertex.
    """
    _check_sequence(sequence)
    if type(disjoint_from_partially_directed) is not bool or type(include_reflection) is not bool:
        raise ValueError("disjoint and reflection flags must be boolean")
    paths = tile_paths(height)
    validate_tile_library(height, paths)
    spectrum = _spectrum_from_library(sequence, height, paths,
        disjoint_from_partially_directed=disjoint_from_partially_directed)
    if include_reflection and len(sequence) >= height*4:
        # Reflect each lattice site over y=0. The first East edge remains
        # canonical, scores are identical and the mirrored walk uses y<0.
        # Every tiled path occupies (0,1), whereas every reflected path
        # occupies (0,-1), proving the two families are disjoint.
        return tuple((score,2*count) for score,count in spectrum)
    return spectrum


@dataclass(frozen=True)
class CompactStripCertificate:
    sequence: str
    height: int
    spectrum: tuple[tuple[int, int], ...]
    disjoint_from_partially_directed: bool = False
    include_reflection: bool = False


@dataclass(frozen=True)
class CompactStripMass:
    partition_lower: Fraction
    event_lower: Fraction
    complement_lower: Fraction
    mean_numerator_lower: Fraction
    second_numerator_lower: Fraction
    witnessed_paths: int
    contact_min: int | None
    contact_max: int | None


def build_compact_strip_certificate(sequence: str, *, height: int = 8,
                                    disjoint_from_partially_directed: bool = False,
                                    include_reflection: bool = False
                                    ) -> CompactStripCertificate:
    return CompactStripCertificate(sequence, height, compact_strip_spectrum(
        sequence, height=height,
        disjoint_from_partially_directed=disjoint_from_partially_directed,
        include_reflection=include_reflection),
        disjoint_from_partially_directed, include_reflection)


def verify_compact_strip_certificate(cert: CompactStripCertificate, q: Fraction,
                                     event_at_least: int = 1) -> CompactStripMass:
    """Fail-closed rational checker: validate tile geometry and replay spectrum.

    Does not trust certificate histogram. Full rerun is independent of the
    certificate producer's claimed values, but shares the Python recurrence.
    """
    if not isinstance(cert, CompactStripCertificate):
        raise ValueError("certificate must be a CompactStripCertificate")
    _check_sequence(cert.sequence)
    if type(cert.height) is not int or not 3 <= cert.height <= 8:
        raise ValueError("invalid tile height")
    if not isinstance(q, Fraction) or q < 1 or type(event_at_least) is not int:
        raise ValueError("expected exact Fraction q >= 1 and integer threshold")
    if type(cert.disjoint_from_partially_directed) is not bool or type(cert.include_reflection) is not bool:
        raise ValueError("disjoint and reflection flags must be boolean")
    paths = tile_paths(cert.height)
    validate_tile_library(cert.height, paths)
    expected = _spectrum_from_library(cert.sequence, cert.height, paths,
        disjoint_from_partially_directed=cert.disjoint_from_partially_directed)
    if cert.include_reflection and len(cert.sequence) >= cert.height*4:
        expected = tuple((score,2*count) for score,count in expected)
    if not isinstance(cert.spectrum, tuple) or expected != cert.spectrum:
        raise ValueError("forged compact strip histogram")
    z = e = complement = m1 = m2 = Fraction(0)
    count = 0
    for k, multiplicity in expected:
        if type(k) is not int or type(multiplicity) is not int or multiplicity < 1:
            raise ValueError("invalid histogram coefficient")
        mass = multiplicity * q ** k
        z += mass
        m1 += k * mass
        m2 += k * k * mass
        count += multiplicity
        if k >= event_at_least:
            e += mass
        else:
            complement += mass
    return CompactStripMass(z, e, complement, m1, m2, count,
                            expected[0][0] if expected else None,
                            expected[-1][0] if expected else None)


def best_compact_strip_mass(sequence: str, q: Fraction, event_at_least: int = 1,
                            heights: tuple[int, ...] = (4, 5, 6, 7, 8)) -> tuple[int, CompactStripMass]:
    """Select best *individual* lower bound (no unjustified family addition)."""
    _check_sequence(sequence)
    if not isinstance(heights, tuple) or not heights:
        raise ValueError("nonempty tuple of tile heights required")
    results = [(height, verify_compact_strip_certificate(
        build_compact_strip_certificate(sequence, height=height), q, event_at_least))
        for height in heights]
    return max(results, key=lambda entry: entry[1].partition_lower)


def compact_plus_directed_mass(sequence: str, q: Fraction, *, height: int = 8,
                               event_at_least: int = 1,
                               include_reflection: bool = True) -> CompactStripMass:
    """Certified (software) sum of *disjoint* PDSAW and compact-tile masses.

    PDSAW is the existing PR #56 union of E-progressive, N-progressive
    with a West step, and S-progressive with a West step. Any path with
    W, N and S steps fails every language, so the compact subset selected
    by the three-direction mask cannot overlap it. The two spectra may
    therefore be ADDED rather than only maximized. This function
    recomputes the PDSAW histogram through the existing research module.
    """
    from .partial_directed import partially_directed_spectrum
    _check_sequence(sequence)
    tile_cert = build_compact_strip_certificate(
        sequence, height=height, disjoint_from_partially_directed=True,
        include_reflection=include_reflection)
    verify_compact_strip_certificate(tile_cert, q, event_at_least)
    combined: dict[int, int] = defaultdict(int)
    for k, c in partially_directed_spectrum(sequence):
        combined[k] += c
    for k, c in tile_cert.spectrum:
        combined[k] += c
    z = e = complement = m1 = m2 = Fraction(0)
    count = 0
    for k, c in combined.items():
        w = c * q**k
        z += w
        m1 += k * w
        m2 += k*k*w
        count += c
        if k >= event_at_least:
            e += w
        else:
            complement += w
    return CompactStripMass(z, e, complement, m1, m2, count,
                            min(combined), max(combined))


@dataclass(frozen=True)
class TailDominanceCertificate:
    """Continuum lower-bound improvement factor for every q >= q_min >= 1.

    The separate disjoint compact spectrum must have no scores below the
    highest score in the existing PDSAW spectrum. Under this restricted
    hypothesis compact/PDSAW is nondecreasing in q, since all compact
    exponents dominate all PDSAW exponents. The finite-sum cross-product
    proves the ordering with nonnegative coefficients.
    """
    sequence: str
    height: int
    q_min: Fraction
    ratio_lower: Fraction
    include_reflection: bool = True


def build_tail_dominance_certificate(sequence: str, *, height: int = 8,
                                     q_min: Fraction = Fraction(10),
                                     include_reflection: bool = True
                                     ) -> TailDominanceCertificate:
    """Propose a rational continuum factor, explicitly refusing invalid cases."""
    from .partial_directed import partially_directed_spectrum
    _check_sequence(sequence)
    if not isinstance(q_min, Fraction) or q_min < 1:
        raise ValueError("q_min must be Fraction >= 1")
    compact = compact_strip_spectrum(sequence, height=height,
        disjoint_from_partially_directed=True,
        include_reflection=include_reflection)
    directed = partially_directed_spectrum(sequence)
    if not compact or not directed or min(k for k,c in compact) < max(k for k,c in directed):
        raise ValueError("contact-score ordering does not establish tail dominance")
    extra = sum((c*q_min**k for k,c in compact), Fraction(0))
    base = sum((c*q_min**k for k,c in directed), Fraction(0))
    return TailDominanceCertificate(sequence, height, q_min, 1 + extra/base,
                                    include_reflection)


def verify_tail_dominance_certificate(cert: TailDominanceCertificate) -> bool:
    """Independently check both spectra, support ordering, and exact ratio."""
    if not isinstance(cert, TailDominanceCertificate):
        raise ValueError("invalid tail dominance certificate")
    if not isinstance(cert.ratio_lower,Fraction) or cert.ratio_lower <= 1:
        raise ValueError("invalid rational ratio claim")
    expected = build_tail_dominance_certificate(cert.sequence,
        height=cert.height, q_min=cert.q_min,
        include_reflection=cert.include_reflection)
    if cert != expected:
        raise ValueError("tail-dominance certificate does not match checked spectra")
    return True
