"""Exact-rational response enclosures for the anchored 2D HP lattice model.

Purely additive: this consumer checks an existing finite frontier using the
reference verifier, then derives bounds on first/second contact moments, the
contact-count variance, and logarithmic Boltzmann-parameter derivatives.

The continuous optimization of an *unknown* subtree's mass is deliberately
relaxed to [0, U]. Thus these are conservative software/model-relative bounds,
not a Lean-verified refinement or real-protein thermodynamics.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING

from .variance_dual import (VarianceDualCertificate, maximize_relaxed_variance,
                            verify_variance_dual)
from .monotone_spectral import (monotone_spectral_histogram,
                                is_monotone_cone_suffix)

if TYPE_CHECKING:
    from .lattice import Certificate


@dataclass(frozen=True)
class ThermalResponseBounds:
    """Enclosures at the certificate's exact rational q >= 1.

    Contact variance is the susceptibility d<E[k]>/d(log q).
    Event susceptibility is dP(k >= t)/d(log q).
    Both derivative identities follow from the finite Gibbs sum and are
    not numerical finite differences. A bound is an interval, not a sampled
    estimate; exact may also mean validated witnesses exhausted upper caps.
    """

    mean_lower: Fraction
    mean_upper: Fraction
    second_moment_lower: Fraction
    second_moment_upper: Fraction
    variance_lower: Fraction
    variance_upper: Fraction
    event_susceptibility_lower: Fraction
    event_susceptibility_upper: Fraction
    partition_lower: Fraction
    partition_upper: Fraction
    support_lower: int
    support_upper: int
    exact_paths: int
    witness_paths: int
    bounded_subtrees: int
    probability_lower: Fraction
    probability_upper: Fraction
    exact: bool
    variance_dual: VarianceDualCertificate
    monotone_paths: int
    monotone_subtrees: int


@dataclass(frozen=True)
class _Unknown:
    min_contacts: int
    max_contacts: int
    mass_upper: Fraction


def _contact_range(sequence: str, path: tuple[tuple[int, int], ...]) -> tuple[int, int]:
    """Independent parity-aware score enclosure for any legal HP prefix.

    On the square lattice, contact partners have odd sequence separation.
    A non-backbone pair must therefore be separated by >= 3 indices. All
    pairs touching an as-yet-unplaced vertex are conservatively possible.
    """
    fixed = len(path)
    minimum = 0
    pending = 0
    for i, residue in enumerate(sequence):
        if residue != "H":
            continue
        for j in range(i + 3, len(sequence), 2):
            if sequence[j] != "H":
                continue
            if j >= fixed:
                pending += 1
            elif abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) == 1:
                minimum += 1
    return minimum, minimum + pending


def _ratio_extreme(known_mass: Fraction, known_weighted_value: Fraction,
                   candidates: tuple[tuple[int, Fraction], ...], *,
                   minimize: bool) -> Fraction:
    """Exact extremum of (M+sum x_i a_i)/(Z+sum x_i), 0<=x_i<=U_i.

    When Z>0, each partial derivative has the sign of (a_i - average).
    Therefore a minimizer fills a prefix of candidates sorted by ascending
    a_i, and a maximizer fills a prefix in descending order. Checking every
    prefix attains the exact relaxed fractional-programming optimum.

    When Z=0 the actual model still has positive partition mass, so the
    conservative extreme is the least/greatest candidate coefficient.
    """
    if known_mass == 0:
        if not candidates:
            raise AssertionError("verified HP frontier has no possible mass")
        values = (value for value, _ in candidates)
        return min(values) if minimize else max(values)
    if not candidates:
        return known_weighted_value / known_mass

    sorted_candidates = sorted(candidates, key=lambda entry: entry[0],
                               reverse=not minimize)
    mass = known_mass
    weighted = known_weighted_value
    best = weighted / mass
    for value, mass_upper in sorted_candidates:
        mass += mass_upper
        weighted += mass_upper * value
        candidate = weighted / mass
        if (candidate < best) if minimize else (candidate > best):
            best = candidate
    return best


def _maximum_bernoulli_variance(lo: Fraction, hi: Fraction) -> Fraction:
    """max p(1-p) for p in [lo, hi], exactly in Q."""
    half = Fraction(1, 2)
    p = min(hi, max(lo, half))
    return p * (1 - p)


def _validate_witness(path: tuple[tuple[int, int], ...], n: int) -> None:
    """Independently validate a terminal SAW; never trust witness coordinates."""
    if (not isinstance(path, tuple) or len(path) != n or
            path[:2] != ((0, 0), (1, 0))):
        raise ValueError("witness must be a full anchored self-avoiding walk")
    if any(not isinstance(site, tuple) or len(site) != 2 or
           any(type(coord) is not int for coord in site) for site in path):
        raise ValueError("witness coordinates must be integer lattice pairs")
    if len(set(path)) != len(path):
        raise ValueError("self-intersecting witness")
    if any(abs(x2-x1) + abs(y2-y1) != 1
           for (x1, y1), (x2, y2) in zip(path, path[1:])):
        raise ValueError("witness contains a non-lattice step")


def thermal_response_bounds(
        cert: Certificate, *,
        witnesses: tuple[tuple[tuple[int, int], ...], ...] = (),
        use_monotone_spectra: bool = False
) -> ThermalResponseBounds:
    """Check an existing frontier, then enclose Gibbs response observables.

    Optional complete-conformation witnesses are independently checked,
    assigned uniquely to a bounded subtree and excluded from its residual
    upper mass, so they increase the known moments without double counting.

    The base verifier independently checks SAW paths, frontier completeness,
    missing children, dead nodes and node status before data are inspected.
    Unknown subtree mass is relaxed to [0, 3**r * q**kmax]; conditioned on
    positive mass, its mean score lies in [kmin,kmax]. This bounds each
    moment by a rational fractional-programming relaxation, rather than by
    dividing uncoupled, potentially loose numerator/denominator intervals.

    Variance lower bound also retains all known-known pairwise differences:
       Z**2 Var(k) = sum_{i<j} w_i w_j (k_i-k_j)**2.
    The event susceptibility uses event/complement cross-pairs. Because
    every event score exceeds every complement score, omitted pairs cannot
    reduce either lower bound. Upper bounds use E[k^2]-E[k]^2 and the
    Bhatia-Davis support/mean inequality (plus a Bernoulli covariance cap).
    """
    from .lattice import verify_certificate

    if type(use_monotone_spectra) is not bool:
        raise ValueError('use_monotone_spectra must be a boolean')
    verified = verify_certificate(cert)
    n = len(cert.sequence)
    q = cert.q
    t = cert.event_at_least
    z_known = m1_known = m2_known = Fraction(0)
    event_mass = event_moment = comp_mass = comp_moment = Fraction(0)
    known_count = 0
    unknown: list[_Unknown] = []
    unknown_paths: list[tuple[tuple[int, int], ...]] = []
    known_scores: list[int] = []
    known_atoms: list[tuple[int, Fraction]] = []

    for node in cert.frontier:
        if node.kind == "dead":
            continue
        low, high = _contact_range(cert.sequence, node.path)
        if len(node.path) == n or node.kind == "exact":
            if low != high:
                raise AssertionError("completed conformation score not exact")
            w = q ** low
            z_known += w
            m1_known += w * low
            m2_known += w * low * low
            if low >= t:
                event_mass += w
                event_moment += w * low
            else:
                comp_mass += w
                comp_moment += w * low
            known_count += 1
            known_scores.append(low)
            known_atoms.append((low, w))
        else:
            upper = Fraction(3 ** (n - len(node.path))) * q ** high
            unknown.append(_Unknown(low, high, upper))
            unknown_paths.append(node.path)

    raw_z_upper = z_known + sum((x.mass_upper for x in unknown), Fraction(0))
    if (z_known != verified.partition_lower or
            raw_z_upper != verified.partition_upper or
            known_count != verified.exact_paths or
            len(unknown) != verified.bounded_subtrees):
        raise AssertionError("thermal-response tally disagrees with base verifier")

    extra_mass = [Fraction(0) for _ in unknown]
    seen_witnesses: set[tuple[tuple[int, int], ...]] = set()
    cone_witness_scores: list[dict[int, int]] = [{} for _ in unknown]
    for path in witnesses:
        _validate_witness(path, n)
        if path in seen_witnesses:
            raise ValueError("duplicate complete-conformation witness")
        seen_witnesses.add(path)
        owning = [i for i, prefix in enumerate(unknown_paths)
                  if path[:len(prefix)] == prefix]
        if len(owning) != 1:
            raise ValueError("witness is not in exactly one bounded subtree")
        index = owning[0]
        k_low, k_high = _contact_range(cert.sequence, path)
        if k_low != k_high:
            raise AssertionError("validated terminal witness has uncertain score")
        w = q ** k_low
        extra_mass[index] += w
        if use_monotone_spectra and is_monotone_cone_suffix(unknown_paths[index], path):
            scores = cone_witness_scores[index]
            scores[k_low] = scores.get(k_low, 0) + 1
        z_known += w
        m1_known += w * k_low
        m2_known += w * k_low ** 2
        known_scores.append(k_low)
        known_atoms.append((k_low, w))
        if k_low >= t:
            event_mass += w
            event_moment += w * k_low
        else:
            comp_mass += w
            comp_moment += w * k_low

    # Count exponentially many genuine monotone extensions in polynomial-state
    # exact-rational DP; exclude all externally credited witnesses in the union.
    monotone_paths = monotone_subtrees = 0
    if use_monotone_spectra:
        for index, prefix in enumerate(unknown_paths):
            spectrum = dict(monotone_spectral_histogram(cert.sequence, prefix))
            if spectrum:
                monotone_subtrees += 1
            for score, count in cone_witness_scores[index].items():
                if count > spectrum.get(score, 0):
                    raise AssertionError('cone overlap exceeds checked spectrum')
            for score, total in spectrum.items():
                if not unknown[index].min_contacts <= score <= unknown[index].max_contacts:
                    raise AssertionError('cone score outside universally checked range')
                fresh = total-cone_witness_scores[index].get(score, 0)
                if fresh < 0:
                    raise AssertionError('negative disjoint cone multiplicity')
                if not fresh:
                    continue
                mass = fresh*q**score
                monotone_paths += fresh
                extra_mass[index] += mass
                z_known += mass
                m1_known += mass*score
                m2_known += mass*score*score
                known_atoms.append((score, mass))
                known_scores.append(score)
                if score >= t:
                    event_mass += mass
                    event_moment += mass*score
                else:
                    comp_mass += mass
                    comp_moment += mass*score

    residual: list[_Unknown] = []
    for node, credited in zip(unknown, extra_mass):
        if credited > node.mass_upper:
            raise AssertionError("validated witness mass exceeds universal cap")
        if credited < node.mass_upper:
            residual.append(_Unknown(node.min_contacts, node.max_contacts,
                                     node.mass_upper - credited))
    unknown = residual
    z_upper = z_known + sum((x.mass_upper for x in unknown), Fraction(0))
    if z_upper != raw_z_upper:
        raise AssertionError("witness accounting changed the partition cap")
    support_min = known_scores + [x.min_contacts for x in unknown]
    support_max = known_scores + [x.max_contacts for x in unknown]
    if not support_min:
        raise AssertionError("every anchored HP model has a straight SAW")

    event_upper, complement_upper = event_mass, comp_mass
    for x in unknown:
        if x.max_contacts >= t:
            event_upper += x.mass_upper
        if x.min_contacts < t:
            complement_upper += x.mass_upper
    prob_lo = (event_mass / (event_mass + complement_upper)
               if event_mass else Fraction(0))
    prob_hi = (event_upper / (event_upper + comp_mass) if comp_mass
               else Fraction(1) if event_upper else Fraction(0))

    candidates_1_lo = tuple((x.min_contacts, x.mass_upper) for x in unknown)
    candidates_1_hi = tuple((x.max_contacts, x.mass_upper) for x in unknown)
    candidates_2_lo = tuple((x.min_contacts ** 2, x.mass_upper) for x in unknown)
    candidates_2_hi = tuple((x.max_contacts ** 2, x.mass_upper) for x in unknown)
    mu_lo = _ratio_extreme(z_known, m1_known, candidates_1_lo, minimize=True)
    mu_hi = _ratio_extreme(z_known, m1_known, candidates_1_hi, minimize=False)
    m2_lo = _ratio_extreme(z_known, m2_known, candidates_2_lo, minimize=True)
    m2_hi = _ratio_extreme(z_known, m2_known, candidates_2_hi, minimize=False)

    variance_known_pair = z_known * m2_known - m1_known ** 2
    event_known_pair = event_moment * comp_mass - event_mass * comp_moment
    if variance_known_pair < 0 or event_known_pair < 0:
        raise AssertionError("known score pairs must have nonnegative gaps")
    var_lo = max(Fraction(0), variance_known_pair / z_upper ** 2,
                 m2_lo - mu_hi ** 2)
    support_lo = min(support_min)
    support_hi = max(support_max)
    middle = Fraction(support_lo + support_hi, 2)
    mu_star = min(mu_hi, max(mu_lo, middle))
    upper_bhatia_davis = (support_hi - mu_star) * (mu_star - support_lo)
    var_hi = min(m2_hi - mu_lo ** 2, upper_bhatia_davis)

    # Exact maximization over the relaxed remaining score/mass polytope.
    # Both rational primal attainment and the universal centered-quadratic
    # dual inequality are checked independently before accepting the cap.
    dual_caps = tuple((x.min_contacts, x.max_contacts, x.mass_upper)
                      for x in unknown)
    dual = maximize_relaxed_variance(tuple(known_atoms), dual_caps)
    if not verify_variance_dual(tuple(known_atoms), dual_caps, dual):
        raise AssertionError('independent primal-dual verifier rejected certificate')
    if dual.upper > var_hi:
        raise AssertionError('sharp relaxed variance exceeds general variance cap')
    var_hi = dual.upper

    susceptibility_lo = event_known_pair / z_upper ** 2
    susceptibility_hi = (support_hi - support_lo) * _maximum_bernoulli_variance(
        prob_lo, prob_hi)

    if not unknown:
        var_lo = var_hi = variance_known_pair / z_known ** 2
        susceptibility_lo = susceptibility_hi = event_known_pair / z_known ** 2

    if not (0 <= prob_lo <= prob_hi <= 1 and
            0 <= support_lo <= mu_lo <= mu_hi <= support_hi and
            0 <= m2_lo <= m2_hi <= support_hi ** 2 and
            0 <= var_lo <= var_hi and
            0 <= susceptibility_lo <= susceptibility_hi):
        raise AssertionError("internal response enclosure inconsistency")
    return ThermalResponseBounds(
        mu_lo, mu_hi, m2_lo, m2_hi, var_lo, var_hi,
        susceptibility_lo, susceptibility_hi, z_known,
        z_upper, support_lo, support_hi,
        known_count, len(seen_witnesses), len(unknown), prob_lo, prob_hi,
        not unknown, dual, monotone_paths, monotone_subtrees)
