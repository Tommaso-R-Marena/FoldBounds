"""Exact-rational transport of *checked* HP certificate bounds across contact weights.

For a finite anchored 2D HP ensemble with w_q(s) = q**K(s), q >= 1,
transporting from q0 to q1 multiplies each weight by (q1/q0)**K(s).
An event K >= threshold and its complement have known integer contact ranges.
These ranges give sound, if potentially loose, bounds without enumerating the
conformations at q1.  The original Python checker, not this module, validates
frontier completeness; this is NOT a Lean-kernel-checked certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

from .lattice import Certificate, VerifiedBounds
from .lookahead import verify_lookahead


@dataclass(frozen=True)
class TransportBounds:
    q: Fraction
    partition_lower: Fraction
    partition_upper: Fraction
    event_lower: Fraction
    event_upper: Fraction
    complement_lower: Fraction
    complement_upper: Fraction
    probability_lower: Fraction
    probability_upper: Fraction

    @property
    def guaranteed_majority(self) -> str:
        """'event', 'complement', or 'undetermined' (strict > 1/2)."""
        if self.probability_lower > Fraction(1, 2):
            return "event"
        if self.probability_upper < Fraction(1, 2):
            return "complement"
        return "undetermined"


@dataclass(frozen=True)
class WindowBounds:
    """Uniform enclosure for *every real* q in [q_lower, q_upper]."""
    q_lower: Fraction
    q_upper: Fraction
    partition_lower: Fraction
    partition_upper: Fraction
    probability_lower: Fraction
    probability_upper: Fraction

    @property
    def guaranteed_majority(self) -> str:
        if self.probability_lower > Fraction(1, 2):
            return "event"
        if self.probability_upper < Fraction(1, 2):
            return "complement"
        return "undetermined"


def _check_target(q: Fraction) -> None:
    if not isinstance(q, Fraction) or q < 1:
        raise ValueError("target q must be an exact Fraction >= 1")


def _contact_cap(sequence: str) -> int:
    """Upper bound on K, computed independently of the checker score routine.

    Any square-lattice contact has odd chain separation >= 3, and two H
    residues. Counting all such pairs is conservative even if geometrically
    incompatible, so this does not assume realizability of the upper cap.
    """
    return sum(sequence[i] == "H" and sequence[j] == "H"
               for i in range(len(sequence))
               for j in range(i + 3, len(sequence), 2))


def _factors(ratio: Fraction, k_min: int, k_max: int) -> tuple[Fraction, Fraction]:
    lo, hi = ratio ** k_min, ratio ** k_max
    return (lo, hi) if ratio >= 1 else (hi, lo)


def _ratio_lower(event_lower: Fraction, other_upper: Fraction) -> Fraction:
    if other_upper == 0:
        return Fraction(1)
    return event_lower / (event_lower + other_upper) if event_lower else Fraction(0)


def _ratio_upper(event_upper: Fraction, other_lower: Fraction) -> Fraction:
    if event_upper == 0:
        return Fraction(0)
    return event_upper / (event_upper + other_lower) if other_lower else Fraction(1)


def _from_checked(sequence: str, q_from: Fraction, threshold: int,
                  bounds: VerifiedBounds, q_target: Fraction) -> TransportBounds:
    """Internal step: caller MUST provide a freshly checked source bound."""
    ratio = q_target / q_from
    max_k = _contact_cap(sequence)
    if threshold <= 0:
        e_range, n_range = (0, max_k), (0, 0)  # no complement states
    elif threshold > max_k:
        e_range, n_range = (0, 0), (0, max_k)  # no event states
    else:
        e_range, n_range = (threshold, max_k), (0, threshold - 1)
    el, eh = _factors(ratio, *e_range)
    nl, nh = _factors(ratio, *n_range)
    zl, zh = _factors(ratio, 0, max_k)
    e_lo, e_hi = bounds.event_lower * el, bounds.event_upper * eh
    n_lo, n_hi = bounds.complement_lower * nl, bounds.complement_upper * nh
    z_lo = max(bounds.partition_lower * zl, e_lo + n_lo)
    z_hi = min(bounds.partition_upper * zh, e_hi + n_hi)
    p_lo, p_hi = _ratio_lower(e_lo, n_hi), _ratio_upper(e_hi, n_lo)
    # Event K >= threshold is monotone in q: every event-state contact
    # count exceeds every complement-state contact count. This can tighten
    # transferred bounds even when mass envelopes are broad.
    if q_target >= q_from:
        p_lo = max(p_lo, bounds.probability_lower)
    if q_target <= q_from:
        p_hi = min(p_hi, bounds.probability_upper)
    if not (0 <= e_lo <= e_hi and 0 <= n_lo <= n_hi
            and 0 <= z_lo <= z_hi and 0 <= p_lo <= p_hi <= 1):
        raise AssertionError("inconsistent transported checked bounds")
    return TransportBounds(q_target, z_lo, z_hi, e_lo, e_hi,
                           n_lo, n_hi, p_lo, p_hi)


def transport_certificate(cert: Certificate, q_target: Fraction, *,
                          lookahead: int = 0, **verify_options: object) -> TransportBounds:
    """Recheck a complete source certificate, then enclose the ensemble at q_target.

    All arithmetic in this operation is exact Fraction arithmetic. The source
    `verify_lookahead` is executed every call: no user-supplied mass numbers,
    contact caps, or claimed refinements are trusted. Optional verifier flags
    (matching/transfer/witnesses etc.) pass through unchanged.
    """
    _check_target(q_target)
    checked = verify_lookahead(cert, lookahead=lookahead, **verify_options)
    return _from_checked(cert.sequence, cert.q, cert.event_at_least, checked, q_target)


def transport_certificates(certificates: Iterable[Certificate], q_target: Fraction, *,
                           lookahead: int = 0, **verify_options: object) -> TransportBounds:
    """Intersect independently checked anchors (same sequence and event only).

    Ancillary transport errors remain correlated; taking max of lower bounds
    and min of upper bounds is sound without an independence assumption.
    """
    _check_target(q_target)
    anchors = tuple(certificates)
    if not anchors:
        raise ValueError("at least one source certificate is required")
    identity = (anchors[0].sequence, anchors[0].event_at_least)
    if any((c.sequence, c.event_at_least) != identity for c in anchors):
        raise ValueError("transport anchors must describe the same HP event")
    parts = [transport_certificate(c, q_target, lookahead=lookahead,
                                   **verify_options) for c in anchors]
    e_lo, e_hi = max(p.event_lower for p in parts), min(p.event_upper for p in parts)
    n_lo, n_hi = max(p.complement_lower for p in parts), min(p.complement_upper for p in parts)
    z_lo = max(*(p.partition_lower for p in parts), e_lo + n_lo)
    z_hi = min(*(p.partition_upper for p in parts), e_hi + n_hi)
    p_lo = max(*(p.probability_lower for p in parts), _ratio_lower(e_lo, n_hi))
    p_hi = min(*(p.probability_upper for p in parts), _ratio_upper(e_hi, n_lo))
    if not (0 <= e_lo <= e_hi and 0 <= n_lo <= n_hi
            and 0 <= z_lo <= z_hi and 0 <= p_lo <= p_hi <= 1):
        raise AssertionError("inconsistent independently checked anchor bounds")
    return TransportBounds(q_target, z_lo, z_hi, e_lo, e_hi,
                           n_lo, n_hi, p_lo, p_hi)


def certify_temperature_window(certificates: Iterable[Certificate],
                               q_lower: Fraction, q_upper: Fraction, *,
                               lookahead: int = 0,
                               **verify_options: object) -> WindowBounds:
    """Bound population and partition mass for all real q in a closed interval.

    q is a *contact-weight* parameter, not temperature in kelvins.
    Since K>=0, Z(q) is nondecreasing; since the event is K>=threshold,
    its probability is nondecreasing as well. Endpoint enclosures therefore
    apply uniformly to the infinite interval, not only sampled q values.
    """
    _check_target(q_lower)
    _check_target(q_upper)
    if q_lower > q_upper:
        raise ValueError("q_lower must not exceed q_upper")
    anchors = tuple(certificates)
    low = transport_certificates(anchors, q_lower, lookahead=lookahead,
                                 **verify_options)
    high = transport_certificates(anchors, q_upper, lookahead=lookahead,
                                  **verify_options)
    result = WindowBounds(q_lower, q_upper, low.partition_lower, high.partition_upper,
                          low.probability_lower, high.probability_upper)
    if not (0 <= result.partition_lower <= result.partition_upper
            and 0 <= result.probability_lower <= result.probability_upper <= 1):
        raise AssertionError("inconsistent whole-interval bounds")
    return result
