"""Exact finite-size and infinite-subsequence density certificates for a 2D HP homopolymer.

This is a *model-relative* Python-rational checker, not a Lean kernel proof.
For square m x m Hamiltonian walks of n=m**2 all-H residues:
  K* = (m-1)**2 contact edges.
Using the exact 16-step SAW count c16, the total number of anchored walks
of length n is bounded above by c16**(n//16) when m is divisible by four.
Combining a single square witness and the entire non-event walk bound gives
an exact, adversarial-reweighting-robust population lower bound.

The arithmetic additionally certifies a *whole infinite subsequence* of
square sizes with a finite set of checks. This is a Peierls-type elementary
entropy-vs-energy argument, not a claim of new polymer mathematics.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import lcm

# Unlike a lookup or a user-supplied bound, the small count is recalculated.
DIRECTIONS = ((1, 0), (-1, 0), (0, 1), (0, -1))


@lru_cache(maxsize=1)
def _exact_anchored_saw_counts16() -> tuple[int, ...]:
    """Recompute all 0..16-step anchored square-lattice SAW counts by exact DFS.

    Only a 16-step free-space tree is enumerated, never the long-chain model.
    The root fixes two adjacent sites. Earlier occupied sites are recorded
    explicitly and cannot be revisited. No floating point is involved.
    """
    occupied = {(0, 0), (1, 0)}
    counts = [0] * 17

    def visit(steps: int, x: int, y: int) -> None:
        counts[steps] += 1
        if steps == 16:
            return
        for dx, dy in DIRECTIONS:
            other = (x + dx, y + dy)
            if other not in occupied:
                occupied.add(other)
                visit(steps + 1, *other)
                occupied.remove(other)

    visit(0, 1, 0)
    return tuple(counts)


def exact_anchored_saw16() -> int:
    """Return the recomputed number of exact 16-step free anchored walks."""
    return _exact_anchored_saw_counts16()[-1]


def square_serpentine(side: int) -> tuple[tuple[int, int], ...]:
    """Independently construct a canonical Hamiltonian square witness."""
    if type(side) is not int or side < 2:
        raise ValueError("side must be an integer >= 2")
    return tuple((x, y) for y in range(side)
                 for x in (range(side) if y % 2 == 0 else
                           range(side - 1, -1, -1)))


def check_square_witness(side: int) -> int:
    """Verify all vertices, steps, and exact contact count of the square path.

    For n=m**2 sites filling a square, the number of unit lattice adjacency
    edges is 2*m*(m-1), of which n-1 are chain bonds. All residues are H.
    Thus every full Hamiltonian square walk has exactly (m-1)**2 contacts.
    """
    path = square_serpentine(side)
    if path[:2] != ((0, 0), (1, 0)) or len(path) != side**2:
        raise AssertionError("bad canonical square-witness generator")
    if len(set(path)) != len(path):
        raise AssertionError("self-intersecting square witness")
    if any(abs(x - a) + abs(y - b) != 1
           for (a, b), (x, y) in zip(path, path[1:])):
        raise AssertionError("nonunit step in square witness")
    occupied = set(path)
    edges = sum((x + 1, y) in occupied for x, y in path)
    edges += sum((x, y + 1) in occupied for x, y in path)
    contacts = edges - (len(path) - 1)
    if contacts != (side - 1) ** 2:
        raise AssertionError("square contact identity failed")
    return contacts


@dataclass(frozen=True)
class SquareTailCertificate:
    """Witness-based, uniform result over every side m=m0+k*step (k>=0)."""
    side_min: int
    side_step: int
    q: int
    density: Fraction
    correction_lower: Fraction
    correction_upper: Fraction
    goal_probability: Fraction
    c16: int
    event_threshold: int
    witness_contacts: int
    walk_count_upper: int
    failure_odds_upper: Fraction
    probability_lower: Fraction
    next_side_ratio: Fraction
    next_ratio_multiplier: Fraction

    @property
    def success(self) -> bool:
        return self.probability_lower >= self.goal_probability


def _valid_inputs(side_min: int, q: int, density: Fraction,
                  correction_lower: Fraction, correction_upper: Fraction,
                  goal_probability: Fraction) -> int:
    if type(q) is not int or q < 2:
        raise ValueError("q must be an integer >=2 for this sufficient-gate checker")
    if not isinstance(density, Fraction) or not 0 < density < 1:
        raise ValueError("density must be an exact Fraction in (0,1)")
    if (not isinstance(correction_lower, Fraction)
            or not isinstance(correction_upper, Fraction)
            or not 0 < correction_lower <= correction_upper):
        raise ValueError("correction factors must be exact positive ordered Fractions")
    if not isinstance(goal_probability, Fraction) or not 0 < goal_probability < 1:
        raise ValueError("goal probability must be an exact Fraction in (0,1)")
    step = lcm(4, density.denominator)
    if type(side_min) is not int or side_min < step or side_min % step:
        raise ValueError("side_min must be a positive multiple of lcm(4,density denominator)")
    return step


def _majorant(side: int, q: int, density: Fraction,
              correction_lower: Fraction, correction_upper: Fraction,
              c16: int) -> tuple[int, int, int, Fraction]:
    n = side * side
    threshold = density * n
    if threshold.denominator != 1 or n % 16:
        raise AssertionError("incompatible square size and exact density")
    t = threshold.numerator
    contacts = (side - 1)**2
    if contacts < t:
        raise ValueError("constructive square walk does not satisfy the event")
    # For m divisible by 4, n=16*b and there are n-2 steps after the root.
    # Split into b-1 blocks of 16 steps, plus one block of 14 steps.
    # Exact short-block enumeration directly checks c14 <= c16; it does
    # not assume every shorter SAW can be extended (trapped tips exist).
    # Therefore number of full walks <= c16**b, even with no shape constraint.
    cap = c16**(n // 16)
    exponent = t - 1 - contacts
    # An unseen non-event has K<=t-1, and there are <= cap-1 such states.
    # We overbound further using cap (not cap-1) for the simple tail proof.
    odds = (correction_upper / correction_lower) * cap * Fraction(q)**exponent
    return t, contacts, cap, odds


def certify_square_density_tail(
    side_min: int, *, q: int = 4, density: Fraction = Fraction(1, 4),
    correction_lower: Fraction = Fraction(1, 2),
    correction_upper: Fraction = Fraction(2),
    goal_probability: Fraction = Fraction(99, 100),
) -> SquareTailCertificate:
    """Sound population certificate on an *infinite arithmetic progression*.

    For all sides m=side_min+k*step, k>=0, where
    step=lcm(4,density.denominator), the event
        K >= density * (m*m)
    has population at least probability_lower for every choice of per-state
    weight multipliers in [lower,upper] at that m.

    The proof is elementary: 1) existence of a Hamiltonian serpentine
    witness with K*=(m-1)^2; 2) short-block SAW submultiplicative cap;
    3) a bound on every unobserved non-event score K<t; 4) an exact
    rational majorant F(m) of non-event/event odds; 5) the ratio
    F(m+step)/F(m) is positive and decreases geometrically when the
    checked next_ratio_multiplier is <1. It suffices to require the
    initial ratio <=1. No finite grid sampling of larger sizes is used.
    """
    step = _valid_inputs(side_min, q, density, correction_lower,
                         correction_upper, goal_probability)
    c16 = exact_anchored_saw16()
    counts = _exact_anchored_saw_counts16()
    if c16 <= 0 or counts[14] > c16:
        raise AssertionError("invalid exact free-space block upper count")
    # Sufficient leading-coefficient gate: c16^(1/16) < q^(1-density).
    a, b = density.numerator, density.denominator
    if c16**b >= q**(16*(b-a)):
        raise ValueError("leading entropy/energy inequality is not certified")
    # Anchor the constructive identity with direct coordinate validation.
    contact_count = check_square_witness(side_min)
    t, k, walk_cap, f0 = _majorant(
        side_min, q, density, correction_lower, correction_upper, c16)
    if k != contact_count:
        raise AssertionError("source witness does not match symbolic contact count")
    _, _, _, f1 = _majorant(
        side_min + step, q, density, correction_lower, correction_upper, c16)
    _, _, _, f2 = _majorant(
        side_min + 2*step, q, density, correction_lower, correction_upper, c16)
    ratio = f1 / f0
    ratio_progression = (f2 / f1) / ratio
    if not (0 < ratio <= 1 and 0 < ratio_progression < 1):
        raise ValueError("whole-tail monotonicity arithmetic gate failed")
    p = 1 / (1 + f0)
    return SquareTailCertificate(
        side_min, step, q, density, correction_lower, correction_upper,
        goal_probability, c16, t, k, walk_cap, f0, p,
        ratio, ratio_progression)


def find_first_square_tail(*, q: int = 4,
                           density: Fraction = Fraction(1, 4),
                           goal_probability: Fraction = Fraction(99, 100),
                           correction_lower: Fraction = Fraction(1, 2),
                           correction_upper: Fraction = Fraction(2),
                           max_side: int = 256) -> SquareTailCertificate | None:
    """Find the earliest *certifiable* side in the chosen congruence class.

    This is NOT a sharp physical onset or phase-transition threshold.
    An unsuccessful trial is not proof that the event majority is absent.
    """
    if type(max_side) is not int or max_side < 4:
        raise ValueError("max_side must be an integer >=4")
    if not isinstance(density, Fraction) or not 0 < density < 1:
        raise ValueError("density must be an exact Fraction in (0,1)")
    step = _valid_inputs(lcm(4, density.denominator), q,
                         density, correction_lower, correction_upper,
                         goal_probability)
    for m in range(step, max_side + 1, step):
        try:
            cert = certify_square_density_tail(
                m, q=q, density=density,
                correction_lower=correction_lower,
                correction_upper=correction_upper,
                goal_probability=goal_probability)
        except ValueError as error:
            if "leading entropy" in str(error):
                return None
            continue
        if cert.success:
            return cert
    return None
