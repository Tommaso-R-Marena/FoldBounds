"""Parametric infinite 3D heterogeneous HP contact-density certificates.

A rigorous sufficient checker for a PERIODIC H/P sequence with one P every
period positions, and the event K_HH>=rho*n (rho exact rational). No full
long-chain enumeration, floating-point estimates or statistical assumptions.

This extends the fixed 1-in-11, rho=1/2 case; the congruence restriction
here is sufficient but not optimized (side multiple of lcm(period,rho.den)).
Both models remain synthetic 3D HP polymers, not natural proteins.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import lcm

from .cubic_density_tail import _local_multiplicity_gate
from .cubic_tiling import _free_saw_counts_3d


@dataclass(frozen=True)
class GeneralCubicDensityTail:
    side_min: int
    side_step: int
    polar_period: int
    density: Fraction
    q: int
    residues: int
    polar_residues: int
    witness_family_lower: int
    all_walks_upper: int
    hh_contacts_lower: int
    event_threshold: int
    failure_odds_upper: Fraction
    probability_lower: Fraction

    @property
    def guaranteed_majority(self) -> bool:
        return self.probability_lower > Fraction(1, 2)


def certify_general_cubic_density_tail(
    side_min: int, *, polar_period: int, density: Fraction, q: int,
    correction_lower: Fraction = Fraction(1, 2),
    correction_upper: Fraction = Fraction(2),
    goal_probability: Fraction = Fraction(99, 100),
) -> GeneralCubicDensityTail:
    """Certify ALL sides m=side_min+k*side_step, k>=0, by one exact gate.

    The periodic sequence has P at i=0,p,2p,... and H otherwise.
    Requiring macro side m to be a multiple of step=lcm(p,rho.denominator)
    guarantees both P count n/p and event threshold rho*n are integers.

    A fixed full-box contact budget gives, for every constructive family walk,

        K_HH >= (16 - 32/p)*m^3 - 12*m^2.

    The event K_HH >= rho*n has t=rho*8*m^3. All non-events have K<=t-1.
    Thus with c8 the exact 8-step cubic SAW count and C>=2*4**(m^3):

        N/E <= (u/(2*l))*(c8/4)^(m^3) *
               q^[-8*(2-rho-4/p)*m^3+12*m^2-1].

    For any fixed step, Delta(m^2)/Delta(m^3)<=1/m<=1/m0.
    Let gain=8*(2-rho-4/p)-12/m0. The exact integer inequality

        c8**gain.denominator <=
           4**gain.denominator * q**gain.numerator

    implies this odds upper bound never increases at later eligible sides.
    The same positive gain implies decay to zero at infinite size.

    These are *sufficient mathematical inequalities* only. Failure to pass
    the gate says nothing about the true thermodynamic event probability.
    """
    if type(polar_period) is not int or polar_period < 2:
        raise ValueError("polar_period must be an integer >=2")
    if not isinstance(density, Fraction) or not 0 < density < 2:
        raise ValueError("density must be an exact Fraction in (0,2)")
    step = lcm(polar_period, density.denominator)
    if type(side_min) is not int or side_min < step or side_min % step:
        raise ValueError("side_min must be a positive multiple of the required step")
    if type(q) is not int or q < 2:
        raise ValueError("contact weight q must be an integer >=2")
    if (not isinstance(correction_lower, Fraction)
            or not isinstance(correction_upper, Fraction)
            or not 0 < correction_lower <= correction_upper):
        raise ValueError("corrections must be ordered positive exact Fractions")
    if not isinstance(goal_probability, Fraction) or not 0 < goal_probability < 1:
        raise ValueError("goal_probability must be a Fraction in (0,1)")
    _local_multiplicity_gate()
    counts = _free_saw_counts_3d(8)
    c8 = counts[8]
    if not 0 < counts[6] <= c8:
        raise AssertionError("invalid exact 3D free-walk block count")

    tiles = side_min**3
    n = 8 * tiles
    polar = n // polar_period
    if n % polar_period:
        raise AssertionError("sequence period not exactly aligned")
    threshold_q = density * n
    if threshold_q.denominator != 1:
        raise AssertionError("contact threshold is not an exact integer")
    threshold = threshold_q.numerator
    # Every complete Hamiltonian box walk has K_all=16B-12m^2+1,
    # and P sites can destroy no more than four contacts each, except
    # the first P endpoint destroys at most one extra contact.
    kmin = 16 * tiles - 12 * side_min**2 - 4 * polar
    if kmin < threshold:
        raise ValueError("no guaranteed event member from the contact budget")

    coefficient = 2 - density - Fraction(4, polar_period)
    gain = 8 * coefficient - Fraction(12, side_min)
    if gain <= 0:
        raise ValueError("no positive asymptotic energy/entropy margin")
    if c8**gain.denominator > (
        4**gain.denominator * q**gain.numerator
    ):
        raise ValueError("infinite-side monotonicity integer gate not certified")

    # Exactly regenerated local motif counts imply a constructive distinct
    # family of 2*4**B rooted Hamiltonian paths for every B>=2.
    family_lower = 2 * 4**tiles
    upper_walks = c8**tiles
    if family_lower > upper_walks:
        raise AssertionError("family count violates universal walk cap")
    odds = (
        correction_upper / correction_lower
        * Fraction(upper_walks, family_lower)
        * Fraction(q)**(threshold - 1 - kmin)
    )
    probability = 1 / (1 + odds)
    if probability < goal_probability:
        raise ValueError("requested probability threshold not certified")
    return GeneralCubicDensityTail(
        side_min, step, polar_period, density, q, n, polar,
        family_lower, upper_walks, kmin, threshold, odds, probability
    )
