"""Scalable, fail-closed strong-coupling certificates for the anchored 2D HP model.

A small verified witness set and an independently computed *upper bound* on the
number of anchored self-avoiding walks imply rigorous population bounds for
arbitrarily large chain lengths. All arithmetic is exact; no exhaustive
enumeration of large conformational ensembles is needed.

This is a Python-checked mathematical claim about a simplified lattice
model, NOT a Lean-kernel proof or a prediction about real proteins.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache

Point = tuple[int, int]
Path = tuple[Point, ...]
ROOT: Path = ((0, 0), (1, 0))
DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1))


@dataclass(frozen=True)
class StrongCouplingCertificate:
    sequence: str
    q: Fraction
    event_at_least: int
    witnesses: tuple[Path, ...]
    correction_lower: Fraction = Fraction(1)
    correction_upper: Fraction = Fraction(1)


@dataclass(frozen=True)
class StrongCouplingResult:
    sequence_length: int
    q: Fraction
    event_at_least: int
    witness_count: int
    event_witness_count: int
    max_possible_contacts: int
    all_walk_count_upper: int
    unseen_walk_count_upper: int
    event_lower: Fraction
    event_upper: Fraction
    complement_lower: Fraction
    complement_upper: Fraction
    partition_lower: Fraction
    partition_upper: Fraction
    probability_lower: Fraction
    probability_upper: Fraction

    @property
    def guaranteed_event_majority(self) -> bool:
        return self.probability_lower > Fraction(1, 2)

    @property
    def guaranteed_complement_majority(self) -> bool:
        return self.probability_upper < Fraction(1, 2)


def _valid_fraction(value: Fraction, *, positive: bool = False) -> bool:
    return isinstance(value, Fraction) and (value > 0 if positive else value >= 1)


def _validate_certificate(cert: StrongCouplingCertificate) -> None:
    if not isinstance(cert, StrongCouplingCertificate):
        raise ValueError("expected a StrongCouplingCertificate")
    if (not isinstance(cert.sequence, str) or len(cert.sequence) < 2
            or set(cert.sequence) - {"H", "P"}):
        raise ValueError("sequence must have at least two H/P residues")
    if not _valid_fraction(cert.q):
        raise ValueError("q must be an exact Fraction >= 1")
    if type(cert.event_at_least) is not int:
        raise ValueError("event threshold must be an integer")
    if (not _valid_fraction(cert.correction_lower, positive=True)
            or not _valid_fraction(cert.correction_upper, positive=True)
            or cert.correction_lower > cert.correction_upper):
        raise ValueError("require exact positive correction bounds lower <= upper")
    if not isinstance(cert.witnesses, tuple):
        raise ValueError("witness collection must be a tuple")
    seen: set[Path] = set()
    for path in cert.witnesses:
        if (not isinstance(path, tuple) or len(path) != len(cert.sequence)
                or path[:2] != ROOT):
            raise ValueError("witness must be a complete canonical-root walk")
        if any(not isinstance(point, tuple) or len(point) != 2
               or any(type(x) is not int for x in point) for point in path):
            raise ValueError("witness coordinates must be integer pairs")
        if len(set(path)) != len(path):
            raise ValueError("self-intersecting witness")
        if any(abs(x - a) + abs(y - b) != 1
               for (a, b), (x, y) in zip(path, path[1:])):
            raise ValueError("non-adjacent steps in witness")
        if path in seen:
            raise ValueError("duplicate witness")
        seen.add(path)


def _contact_count(sequence: str, path: Path) -> int:
    return sum(sequence[i] == "H" and sequence[j] == "H"
               and abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) == 1
               for i in range(len(sequence))
               for j in range(i + 2, len(sequence)))


def _contact_cap(sequence: str) -> int:
    """Every noncovalent 2D lattice contact has odd chain distance >=3."""
    return sum(sequence[i] == "H" and sequence[j] == "H"
               for i in range(len(sequence))
               for j in range(i + 3, len(sequence), 2))


@lru_cache(maxsize=8)
def _free_block_counts(max_steps: int) -> tuple[int, ...]:
    """EXACT small free-space counts for k steps after a fixed first bond.

    Starting from ROOT, enumerate up to eight further non-repeating moves.
    With any preceding path, legal extensions embed in this free-space set
    by translation and rotation; old obstacles only remove extensions.
    """
    frontier = [ROOT]
    counts = [1]
    for _ in range(max_steps):
        nxt: list[Path] = []
        for path in frontier:
            x, y = path[-1]
            occupied = set(path)
            for dx, dy in DIRECTIONS:
                p = (x + dx, y + dy)
                if p not in occupied:
                    nxt.append(path + (p,))
        frontier = nxt
        counts.append(len(frontier))
    return tuple(counts)


@lru_cache(maxsize=2048)
def anchored_walk_count_upper(length: int, *, block_steps: int = 8) -> int:
    """Exact upper bound by optimal product of free-space block counts.

    After any SAW prefix ending in a bond, its next k legal steps are a
    subset of the anchored free-space k-step walks. Their number is at
    most c[k]. Therefore any composition r=k1+...+km of remaining steps
    gives an upper bound prod(c[ki]) for the full walk tree.

    Dynamic programming minimizes over all block partitions with ki<=8.
    This is a universal *count bound*, not a claim to enumerate all paths.
    """
    if type(length) is not int or length < 2:
        raise ValueError("length must be an integer >= 2")
    if type(block_steps) is not int or not 1 <= block_steps <= 8:
        raise ValueError("block_steps must be an integer in [1,8]")
    counts = _free_block_counts(block_steps)
    remaining = length - 2
    best = [1]
    for r in range(1, remaining + 1):
        best.append(min(counts[k] * best[r - k]
                        for k in range(1, min(block_steps, r) + 1)))
    return best[-1]


def check_strong_coupling(cert: StrongCouplingCertificate, *,
                          block_steps: int = 8) -> StrongCouplingResult:
    """Verify witnesses, then bound all missing conformations adversarially.

    A full ensemble state has nominal weight q**K; its actual weight may be
    independently multiplied by ANY factor in [correction_lower,
    correction_upper], including fully correlated/adversarial choices.
    A known event witness contributes at least lower*q**K. At most
    (universal walk cap - number of witnesses) other states remain.
    Every non-event K < event_at_least, so their unobserved total mass
    is at most unseen_count*upper*q**min(event_at_least-1, M).
    This yields a rigorous probability lower bound without exhaustive search.
    """
    _validate_certificate(cert)
    walk_cap = anchored_walk_count_upper(len(cert.sequence), block_steps=block_steps)
    if len(cert.witnesses) > walk_cap:
        raise AssertionError("checked witnesses exceed universal walk cap")
    unseen = walk_cap - len(cert.witnesses)
    maximum = _contact_cap(cert.sequence)
    known_event = known_complement = Fraction(0)
    event_witness_count = 0
    for path in cert.witnesses:
        score = _contact_count(cert.sequence, path)
        if score > maximum:
            raise AssertionError("parity contact cap is inconsistent")
        weight = cert.q ** score
        if score >= cert.event_at_least:
            known_event += weight
            event_witness_count += 1
        else:
            known_complement += weight

    event_possible = cert.event_at_least <= maximum
    complement_possible = cert.event_at_least > 0
    event_missing_upper = unseen * cert.q**maximum if event_possible else Fraction(0)
    non_event_max = min(maximum, cert.event_at_least - 1)
    complement_missing_upper = (unseen * cert.q**non_event_max
                                if complement_possible else Fraction(0))

    l, u = cert.correction_lower, cert.correction_upper
    e_lo = l * known_event
    e_hi = u * (known_event + event_missing_upper)
    n_lo = l * known_complement
    n_hi = u * (known_complement + complement_missing_upper)
    z_lo = e_lo + n_lo
    z_hi = u * (known_event + known_complement + unseen * cert.q**maximum)

    p_lo = (Fraction(1) if n_hi == 0 else
            e_lo / (e_lo + n_hi) if e_lo else Fraction(0))
    p_hi = (Fraction(0) if e_hi == 0 else
            e_hi / (e_hi + n_lo) if n_lo else Fraction(1))
    if not (0 <= e_lo <= e_hi and 0 <= n_lo <= n_hi
            and 0 <= z_lo <= z_hi and 0 <= p_lo <= p_hi <= 1):
        raise AssertionError("inconsistent exact-rational strong-coupling envelope")
    return StrongCouplingResult(
        len(cert.sequence), cert.q, cert.event_at_least,
        len(cert.witnesses), event_witness_count, maximum, walk_cap, unseen,
        e_lo, e_hi, n_lo, n_hi, z_lo, z_hi, p_lo, p_hi)


def minimum_integer_majority_q(cert: StrongCouplingCertificate, *,
                               max_q: int = 1_000_000,
                               block_steps: int = 8) -> int | None:
    """Find smallest integer q>=1 for which this checker certifies event majority.

    No claim that this is the true phase-transition point: it is the threshold
    of the conservative witness/count certificate. If no event witness exists
    or the threshold exceeds max_q, return None.

    Monotonic predicate proof: each known event score K>=t, and each term
    in the complement upper bound has K<t. Thus E_lower/N_upper
    increases (or is constant for an empty complement) as q increases.
    Binary search therefore never skips an earlier certifiable integer q.
    """
    _validate_certificate(cert)
    if type(max_q) is not int or max_q < 1:
        raise ValueError("max_q must be a positive integer")
    def succeeds(q: int) -> bool:
        updated = replace(cert, q=Fraction(q))
        return check_strong_coupling(updated, block_steps=block_steps).guaranteed_event_majority
    if not succeeds(max_q):
        return None
    lo, hi = 1, max_q
    while lo < hi:
        mid = (lo + hi) // 2
        if succeeds(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def certify_majority_halfline(cert: StrongCouplingCertificate, *,
                              block_steps: int = 8) -> StrongCouplingResult:
    """Prove event probability >1/2 for EVERY real q >= cert.q.

    The uncertainty multipliers may independently vary with q and with
    the microstate. All known event powers have K>=t, and all possible
    complement powers have K<t, so the certified lower event/upper
    complement ratio cannot decrease when q increases.
    """
    result = check_strong_coupling(cert, block_steps=block_steps)
    if not result.guaranteed_event_majority:
        raise ValueError("not certified: event lower mass <= complement upper mass")
    return result


def serpentine_witness(width: int, height: int) -> Path:
    """Untrusted candidate constructor for an all-H rectangular grid walk.

    Its output is never accepted without the same independent path validation
    as any other submitted witness.
    """
    if type(width) is not int or type(height) is not int or width < 2 or height < 1:
        raise ValueError("width >= 2 and height >= 1 required")
    return tuple((x, y) for y in range(height)
                 for x in (range(width) if y % 2 == 0 else range(width - 1, -1, -1)))
