"""Infinite-subsequence 3D heterogeneous HP concentration, exact-rational gate.

A finite exhaustive eight-vertex cube enumeration supplies >=4 Hamiltonian
walks for EVERY opposite-parity endpoint pair, >=2 rooted walks for a
fixed first exit, and >=16 terminal walks for any fixed even start.
The fixed intertile bridge parity allows composition of an exponential
family, so on B>=2 macro tiles C >= 2*4**B. This is a LOWER bound on
distinct complete rooted 3D SAWs, not a claim to count the whole ensemble.

For macro cubes (m,m,m), m divisible by 11, n=8*m**3, and sequence
P,H*10 repeated, P count is exactly n/11. The fully occupied cubic
rectangle supplies HH >= (144/11)*m**3-12*m**2 for every family path.
For A={K_HH >= n/2}, t=4*m**3, an exact finite free-step count gives
all SAWs <= c8**(m**3) and therefore a rational failure-odds bound.

One finite integer inequality provides a *uniform tail* guarantee for all
m=m0+11*k, k>=0. No large-m ensemble or family is enumerated.

This is a Python-checked combinatorial proof for a synthetic HP toy
model, NOT Lean-verified or an experimental protein folding result.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

from .cubic_tiling import CUBE, ROOT, _free_saw_counts_3d, tile_paths


@dataclass(frozen=True)
class CubicDensityTail:
    side_min: int
    step: int
    q: int
    residues: int
    polar_residues: int
    family_count_lower: int
    all_walk_count_upper: int
    hh_contacts_lower: int
    event_threshold: int
    failure_odds_upper: Fraction
    probability_lower: Fraction
    c8: int

    @property
    def certified(self) -> bool:
        return self.probability_lower > Fraction(1, 2)


def _local_multiplicity_gate() -> None:
    """Recount every local opposite-parity endpoint pair independently."""
    motifs = tile_paths()
    counts = Counter((path[0], path[-1]) for path in motifs)
    even = tuple(p for p in CUBE if sum(p) % 2 == 0)
    odd = tuple(p for p in CUBE if sum(p) % 2 == 1)
    if not all(counts[a, b] >= 4 for a in even for b in odd):
        raise AssertionError("local Hamiltonian endpoint multiplicity gate failed")
    # The first macro tile always exits through its +x face for m>=2.
    if sum(path[:2] == ROOT and path[-1] == (1, 1, 1)
           for path in motifs) < 2:
        raise AssertionError("rooted-first-tile multiplicity gate failed")
    if any(sum(counts[a, b] for b in odd) < 16 for a in even):
        raise AssertionError("last-tile free-end multiplicity gate failed")


def certify_periodic_cubic_density_tail(
    side_min: int, *, q: int,
    correction_lower: Fraction = Fraction(1, 2),
    correction_upper: Fraction = Fraction(2),
    goal_probability: Fraction = Fraction(99, 100),
) -> CubicDensityTail:
    """Check one finite gate and certify all cubic sides m=side_min+11*k.

    Every model uses a 3D cubic-lattice H/P sequence with P at precisely
    positions 0,11,22,... and all other residues H. The coarse event
    is K_HH >= n/2, i.e. >=4*m**3 noncovalent H-H nearest-neighbor
    contacts, NOT a claim that half of residues participate in contacts.

    The analytical upper bound on non-event/event mass is

      F(m) = (u/(2*l)) * (c8/4)^(m**3)
                    * q^[-100*m**3/11 + 12*m**2 - 1].

    m is divisible by 11 so every exponent is integer. For every
    m>=m0, the ratio Delta(m**2)/Delta(m**3) <= 1/m0.
    Therefore the *exact integer* gate

      (c8/4)^(11*m0) <= q^(100*m0-132)

    implies F(m+11)<=F(m) for ALL eligible m>=m0.
    The same strict asymptotic gap forces F->0 along the progression.
    All microstates may be adversarially reweighted independently
    within the stated positive exact-rational correction box.
    """
    if type(side_min) is not int or side_min < 11 or side_min % 11:
        raise ValueError("side_min must be an integer multiple of 11, >=11")
    if type(q) is not int or q < 2:
        raise ValueError("q must be an integer >= 2")
    if (not isinstance(correction_lower, Fraction)
            or not isinstance(correction_upper, Fraction)
            or not 0 < correction_lower <= correction_upper):
        raise ValueError("correction bounds must be ordered positive exact Fractions")
    if not isinstance(goal_probability, Fraction) or not 0 < goal_probability < 1:
        raise ValueError("goal_probability must be an exact Fraction in (0,1)")

    _local_multiplicity_gate()
    counts = _free_saw_counts_3d(8)
    c8 = counts[-1]
    if not 0 < counts[6] <= c8:
        raise AssertionError("exact 6-step SAW cap invalid")
    tiles = side_min**3
    residues = 8 * tiles
    if residues % 11 != 0:
        raise AssertionError("polar periodicity arithmetic failed")
    polar = residues // 11
    threshold = 4 * tiles
    # Fully occupied (2m)^3 box contains exactly 16m^3-12m^2+1
    # nonbonded contacts. Every P internal to the chain removes at
    # most four HH contact slots; the initial P removes at most five,
    # the final residue is H because n is a multiple of 11.
    contact_lower = 16 * tiles - 12 * side_min**2 - 4 * polar
    if contact_lower != 144 * tiles // 11 - 12 * side_min**2:
        raise AssertionError("3D geometric HH contact floor inconsistent")
    if contact_lower < threshold:
        raise ValueError("constructive family is not guaranteed inside the event")
    # Constructively selected one bridge per macro interface:
    # first tile >=2 valid paths, internal tiles >=4 per pair,
    # terminal tile >=16 free-end paths. B>=2 yields >=2*4**B.
    family_lower = 2 * 4**tiles
    # n-2=8(B-1)+6 free steps. c6 <= c8, so W(n)<=c8**B.
    total_upper = c8**tiles
    if family_lower > total_upper:
        raise AssertionError("family count exceeds universal count cap")

    # Sufficient monotonically decreasing-tail gate. Raising a
    # rational inequality to 11*m0 clears ALL exponents.
    if c8**(11 * side_min) > (
        4**(11 * side_min) * q**(100 * side_min - 132)
    ):
        raise ValueError("infinite-tail monotonicity gate not certified")

    odds = (
        (correction_upper / correction_lower)
        * Fraction(total_upper, family_lower)
        * Fraction(q)**(threshold - 1 - contact_lower)
    )
    probability = 1 / (1 + odds)
    if probability < goal_probability:
        raise ValueError("target population lower bound not certified")
    return CubicDensityTail(
        side_min, 11, q, residues, polar, family_lower, total_upper,
        contact_lower, threshold, odds, probability, c8
    )
