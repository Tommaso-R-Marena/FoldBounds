"""Auditable exponential witness families for 3D cubic HP lattice ensembles.

Each 2x2x2 tile admits finitely many Hamiltonian walks, enumerated by a
self-contained small DFS. A Hamiltonian *macro-tile* route and an integer
transfer calculation count exponentially many distinct complete walks
without enumerating that full family. Every family walk is a genuine
anchored, self-avoiding walk in Z^3.

A separate 3D self-avoiding-walk block count bounds ALL competing walks.
Together with contact lower bounds and arbitrary positive microstate
reweighting, this gives exact-rational population certificates.

NOT a Lean kernel checked algorithm, real-protein prediction, or historical
claim of novel transfer-matrix or HP model mathematics.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache

Point = tuple[int, int, int]
Path = tuple[Point, ...]
OFFSETS: tuple[Point, ...] = (
    (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)
)
CUBE: tuple[Point, ...] = tuple(
    (x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)
)
ROOT: Path = ((0, 0, 0), (1, 0, 0))


@lru_cache(maxsize=1)
def tile_paths() -> tuple[Path, ...]:
    """Independently enumerate ALL ordered 8-site cube Hamiltonian walks."""
    adjacent = {
        p: tuple(
            (p[0] + dx, p[1] + dy, p[2] + dz)
            for dx, dy, dz in OFFSETS
            if (p[0] + dx, p[1] + dy, p[2] + dz) in CUBE
        ) for p in CUBE
    }
    output: list[Path] = []

    def visit(path: Path, seen: set[Point]) -> None:
        if len(path) == 8:
            output.append(path)
            return
        for candidate in adjacent[path[-1]]:
            if candidate not in seen:
                seen.add(candidate)
                visit(path + (candidate,), seen)
                seen.remove(candidate)

    for start in CUBE:
        visit((start,), {start})
    if not output or len(set(output)) != len(output):
        raise AssertionError("invalid or duplicated 8-site tile paths")
    return tuple(output)


def macro_snake(nx: int, ny: int, nz: int) -> tuple[Point, ...]:
    """Hamiltonian traversal of a rectangular macro-cube grid."""
    if any(type(v) is not int or v < 1 for v in (nx, ny, nz)):
        raise ValueError("positive integer macro dimensions required")
    output: list[Point] = []
    for z in range(nz):
        ys = range(ny) if z % 2 == 0 else range(ny - 1, -1, -1)
        for y in ys:
            xs = range(nx) if (len(output) // nx) % 2 == 0 else range(nx - 1, -1, -1)
            for x in xs:
                output.append((x, y, z))
    if output[0] != (0, 0, 0) or len(output) != nx * ny * nz:
        raise AssertionError("invalid macro route")
    if len(set(output)) != len(output):
        raise AssertionError("duplicate macro tile")
    if any(sum(abs(a - b) for a, b in zip(u, v)) != 1
           for u, v in zip(output, output[1:])):
        raise AssertionError("disconnected macro route")
    return tuple(output)


def _local_contacts(block: str, path: Path) -> int:
    lookup = {p: i for i, p in enumerate(path)}
    contacts = 0
    for i, (x, y, z) in enumerate(path):
        if block[i] != "H":
            continue
        for dx, dy, dz in OFFSETS:
            j = lookup.get((x + dx, y + dy, z + dz))
            if j is not None and j > i + 1 and block[j] == "H":
                contacts += 1
    return contacts


@lru_cache(maxsize=1024)
def _tile_options(block: str, first: bool = False) -> tuple[tuple[Point, Point, int, int], ...]:
    """(start, end, local_HH, multiplicity) over validated local paths."""
    if len(block) != 8 or set(block) - {"H", "P"}:
        raise ValueError("each tile must be eight H/P residues")
    grouped: dict[tuple[Point, Point, int], int] = defaultdict(int)
    for path in tile_paths():
        if first and path[:2] != ROOT:
            continue
        grouped[(path[0], path[-1], _local_contacts(block, path))] += 1
    return tuple(
        (start, end, k, multiplicity)
        for (start, end, k), multiplicity in sorted(grouped.items())
    )


def _bridge(u: Point, v: Point, before: Point, after: Point) -> bool:
    """Unit bridge between adjacent macro tiles; their vertex sets are disjoint."""
    return sum(
        abs(2 * a + x - 2 * b - y)
        for a, x, b, y in zip(before, u, after, v)
    ) == 1


@lru_cache(maxsize=9)
def _free_saw_counts_3d(block_steps: int = 8) -> tuple[int, ...]:
    """Exact small fixed-bond free 3D self-avoiding-walk counts (not a trusted table)."""
    if type(block_steps) is not int or not 1 <= block_steps <= 8:
        raise ValueError("block_steps must be an integer in [1,8]")
    occupied = {(0, 0, 0), (1, 0, 0)}
    counts = [0] * (block_steps + 1)

    def dfs(steps: int, point: Point) -> None:
        counts[steps] += 1
        if steps == block_steps:
            return
        for dx, dy, dz in OFFSETS:
            candidate = (point[0] + dx, point[1] + dy, point[2] + dz)
            if candidate not in occupied:
                occupied.add(candidate)
                dfs(steps + 1, candidate)
                occupied.remove(candidate)

    dfs(0, (1, 0, 0))
    return tuple(counts)


@lru_cache(maxsize=2048)
def all_saw_count_upper_3d(length: int, block_steps: int = 8) -> int:
    """Exact integer cap from optimal small-block decomposition of remaining steps.

    For any SAW prefix, old occupied vertices only eliminate possible
    extensions. After lattice translation/rotation, at most c[k] free
    anchored k-step choices remain. Products over consecutive blocks
    therefore upper-bound the complete walk tree.
    """
    if type(length) is not int or length < 2:
        raise ValueError("length must be an integer >=2")
    counts = _free_saw_counts_3d(block_steps)
    best = [1]
    for remaining in range(1, length - 1):
        best.append(
            min(counts[k] * best[remaining - k]
                for k in range(1, min(remaining, block_steps) + 1))
        )
    return best[-1]


@dataclass(frozen=True)
class TiledCertificate:
    macro_dimensions: Point
    sequence: str
    q: Fraction
    event_at_least: int
    correction_lower: Fraction = Fraction(1)
    correction_upper: Fraction = Fraction(1)


@dataclass(frozen=True)
class TiledResult:
    residues: int
    macro_tiles: int
    family_paths: int
    family_min_contacts: int
    coordination_lower_contacts: int
    local_lower_contacts: int
    total_all_h_contacts: int
    total_walks_upper: int
    event_mass_lower: Fraction
    complement_mass_upper: Fraction
    probability_lower: Fraction

    @property
    def guaranteed_event_majority(self) -> bool:
        return self.event_mass_lower > self.complement_mass_upper


def check_tiled_family(cert: TiledCertificate) -> TiledResult:
    """Independently validate exact local motifs, transfer family, and mass bound.

    All complete paths counted by the transfer are distinct and globally
    self-avoiding. Each tile visits its disjoint 8-site set exactly once,
    and each consecutive intertile bridge is one lattice step. For a
    heterogeneous H/P sequence, the sum of intratile H/H contacts is a
    LOWER bound (all omitted intertile contacts are nonnegative).
    Independently, a completely filled rectangular box has a fixed
    total contact edge count; every internal P removes at most four
    noncovalent H/H slots and each endpoint P removes at most five.

    A non-event conformation has K <= t-1. Its total mass is bounded
    by the universal walk count times q**(t-1), and the independently
    permitted correction upper bound. No untrusted score or path count
    is supplied to the checker.
    """
    if not isinstance(cert, TiledCertificate):
        raise ValueError("expected TiledCertificate")
    dims = cert.macro_dimensions
    if not isinstance(dims, tuple) or len(dims) != 3:
        raise ValueError("macro_dimensions must be three positive integers")
    nx, ny, nz = dims
    route = macro_snake(nx, ny, nz)
    length = 8 * len(route)
    if (not isinstance(cert.sequence, str) or len(cert.sequence) != length
            or set(cert.sequence) - {"H", "P"}):
        raise ValueError("sequence must consist of exactly n H/P residues")
    if not isinstance(cert.q, Fraction) or cert.q < 1:
        raise ValueError("q must be exact Fraction >=1")
    if type(cert.event_at_least) is not int or cert.event_at_least < 1:
        raise ValueError("event threshold must be a positive integer")
    if (not isinstance(cert.correction_lower, Fraction)
            or not isinstance(cert.correction_upper, Fraction)
            or not 0 < cert.correction_lower <= cert.correction_upper):
        raise ValueError("correction factors must be ordered positive exact Fractions")

    # Two semirings are propagated over the SAME verified DAG:
    #  1) integer (+,*) for unique path counts,
    #  2) min-plus for smallest intratile HH score among those paths.
    count: dict[Point, int] = {}
    minscore: dict[Point, int] = {}
    for index, macro in enumerate(route):
        block = cert.sequence[index * 8: (index + 1) * 8]
        candidates = _tile_options(block, first=(index == 0))
        new_count: dict[Point, int] = defaultdict(int)
        new_min: dict[Point, int] = {}
        for start, end, local_k, multiplicity in candidates:
            if index == 0:
                precursors = ((1, 0),) if start == (0, 0, 0) else ()
            else:
                precursors = (
                    (ways, minscore[previous_end])
                    for previous_end, ways in count.items()
                    if _bridge(previous_end, start, route[index - 1], macro)
                )
            for previous_count, previous_score in precursors:
                new_count[end] += previous_count * multiplicity
                score = previous_score + local_k
                if end not in new_min or score < new_min[end]:
                    new_min[end] = score
        if not new_count:
            raise ValueError("no canonical tile-compatible Hamiltonian family")
        count, minscore = dict(new_count), new_min

    family = sum(count.values())
    local_lower = min(minscore.values())
    sx, sy, sz = (2 * nx, 2 * ny, 2 * nz)
    grid_edges = (sx - 1) * sy * sz + sx * (sy - 1) * sz + sx * sy * (sz - 1)
    all_h_contacts = grid_edges - (length - 1)
    polar_count = cert.sequence.count("P")
    coordination_lower = max(
        0, all_h_contacts - 4 * polar_count
        - int(cert.sequence[0] == "P") - int(cert.sequence[-1] == "P")
    )
    minimum_contacts = max(local_lower, coordination_lower)
    upper_count = all_saw_count_upper_3d(length)
    if family > upper_count or minimum_contacts > all_h_contacts:
        raise AssertionError("inconsistent verified witness family")

    event_lower = (
        cert.correction_lower * family * cert.q ** minimum_contacts
        if minimum_contacts >= cert.event_at_least else Fraction(0)
    )
    other_upper = (
        cert.correction_upper * upper_count * cert.q ** (cert.event_at_least - 1)
    )
    probability_lower = (
        event_lower / (event_lower + other_upper)
        if event_lower else Fraction(0)
    )
    return TiledResult(
        length, len(route), family, minimum_contacts, coordination_lower,
        local_lower, all_h_contacts, upper_count,
        event_lower, other_upper, probability_lower
    )


def minimum_certified_integer_q(cert: TiledCertificate, max_q: int = 100) -> int | None:
    """Find first sufficient integer q; NOT the physical folding transition.

    The lower event mass grows as q**Kmin, and the upper competing
    non-event mass grows as q**(t-1), so if Kmin>=t the ratio is
    nondecreasing for real q>=1. Binary search cannot skip success.
    """
    if type(max_q) is not int or max_q < 1:
        raise ValueError("max_q must be a positive integer")

    def success(q: int) -> bool:
        return check_tiled_family(replace(cert, q=Fraction(q))).guaranteed_event_majority

    if not success(max_q):
        return None
    lo, hi = 1, max_q
    while lo < hi:
        mid = (lo + hi) // 2
        if success(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def certify_majority_halfline(cert: TiledCertificate) -> TiledResult:
    """Sufficient strict majority certificate for EVERY real q>=cert.q.

    This is valid even when corrections vary arbitrarily by state and q,
    provided all remain within the globally stated positive interval.
    """
    result = check_tiled_family(cert)
    if not result.guaranteed_event_majority:
        raise ValueError("no certified majority at the stated contact weight")
    return result
