"""Exact global optimization of HP contact-event probabilities for all q>1.

For A_t={K>=threshold}, every conformation has an exact rational q-weight.
P(A_t)=E(q)/Z(q), with Z(q)>0. A comparison of two sequences is the sign
of the integer polynomial E_A*Z_B-E_B*Z_A in x=q-1.

Coefficientwise-positive cross products provide sufficient *all-q*
domination; exact Sturm root isolation certifies the remaining order
across the entire real continuum. Degenerate, overlapping, and resource
limit cases explicitly fail closed rather than guess phase order.

Python reference research only; not a Lean-refined executable or a physical
protein folding predictor.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Literal

from .exact_phase import (
    _diff, _eval, _product, _trim,
    certify_sequence_comparison, q_coefficients_to_t,
)
from .fixed_composition_dominance import sequence_from_mask


@dataclass(frozen=True)
class EventPhaseTransition:
    t_lower: Fraction
    t_upper: Fraction
    first_class: int
    second_class: int
    kind: Literal["global_switch", "global_tangent", "irrelevant"]
    before: tuple[str, ...]
    after: tuple[str, ...]


@dataclass(frozen=True)
class EventProbabilityEnvelope:
    length: int
    h_count: int
    event_at_least: int
    population: int
    coefficient_survivors: int
    skyline_classes: int
    final_classes: tuple[tuple[str, ...], ...]
    root_free_winners: tuple[tuple[str, ...], ...]
    events: tuple[EventPhaseTransition, ...]

    @property
    def switches(self) -> tuple[EventPhaseTransition, ...]:
        return tuple(e for e in self.events if e.kind == "global_switch")


def certify_event_probability_envelope(
    landscape, h_count: int, event_at_least: int, *,
    width: Fraction = Fraction(1, 1_000_000),
    max_sequences: int = 20_000,
    max_classes: int = 128,
    max_events: int = 512,
    max_refinements: int = 4,
) -> EventProbabilityEnvelope:
    """Certify the GLOBAL maximizers of P(K>=threshold) for ALL real q>1.

    Unlike optimizing Z, one cannot infer probability dominance from Z
    coefficient dominance. This routine instead uses exact cross-products
    of the event and partition polynomials and proves their signs.

    It certifies every original sequence is dominated by the final
    cross-coefficient skyline, then every discarded skyline class has
    a strict all-q final dominator. The complete pairwise Sturm root
    inventory of the remaining classes gives all global switches and
    tangent ties. Unseparated roots cause a bounded, fail-closed abstention.
    """
    n = landscape.length
    if (type(n) is not int or not 2 <= n <= 20
            or type(h_count) is not int or not 0 <= h_count <= n
            or type(event_at_least) is not int or event_at_least < 0):
        raise ValueError("invalid chain length, composition or event threshold")
    if type(width) is not Fraction or width <= 0:
        raise ValueError("width must be a positive exact Fraction")
    if any(type(v) is not int or v < 1
           for v in (max_sequences, max_classes, max_events, max_refinements)):
        raise ValueError("work budgets must be positive integers")

    masks = tuple(sum(1 << i for i in indices)
                  for indices in combinations(range(n), h_count))
    if len(masks) > max_sequences:
        raise ValueError("composition enumeration exceeds budget")
    cache: dict[int, tuple[tuple[int, ...], tuple[int, ...]]] = {}
    conformations = None

    def polynomials(mask: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
        nonlocal conformations
        if mask not in cache:
            sequence = sequence_from_mask(mask, n)
            z = _trim(tuple(landscape.polynomial(sequence)))
            e = q_coefficients_to_t(
                landscape.contact_spectrum(sequence), event_at_least
            )
            if (not z or not e or any(type(x) is not int or x < 0 for x in z + e)
                    or z[0] < 1 or any(
                        (e[i] if i < len(e) else 0) >
                        (z[i] if i < len(z) else 0)
                        for i in range(max(len(z), len(e)))
                    )):
                raise ValueError("inconsistent exact event/partition coefficients")
            if conformations is None:
                conformations = z[0]
            if z[0] != conformations:
                raise ValueError("sequences disagree on q=1 conformation count")
            cache[mask] = (e, z)
        return cache[mask]

    def difference(first: int, second: int):
        ea, za = polynomials(first)
        eb, zb = polynomials(second)
        return _diff(_product(ea, zb), _product(eb, za))

    def weak_dominates(first: int, second: int) -> bool:
        # Sufficiency only: nonnegative coefficients in t=q-1 ensure
        # P_first >= P_second for the ENTIRE real q>=1 continuum.
        return all(c >= 0 for c in difference(first, second))

    skyline: list[int] = []
    for mask in masks:
        if any(weak_dominates(old, mask)
               and difference(old, mask) != (0,) for old in skyline):
            continue
        skyline = [
            old for old in skyline
            if not (weak_dominates(mask, old)
                    and difference(mask, old) != (0,))
        ]
        skyline.append(mask)
    # Audit every originally considered sequence against the *final* skyline.
    if not all(any(weak_dominates(old, mask) for old in skyline)
               for mask in masks):
        raise ArithmeticError("final event skyline fails whole-composition coverage")

    # Group sequences with the same RATIONAL event probability, even if
    # their event and partition polynomials differ individually.
    classes: list[tuple[int, list[str]]] = []
    for mask in skyline:
        for representative, names in classes:
            if difference(representative, mask) == (0,):
                names.append(sequence_from_mask(mask, n))
                break
        else:
            classes.append((mask, [sequence_from_mask(mask, n)]))
    if len(classes) > max_classes:
        raise ValueError("too many event-probability skyline classes")

    comparisons = {}
    discarded = set()
    for i in range(len(classes)):
        for j in range(i + 1, len(classes)):
            cert = certify_sequence_comparison(
                landscape, classes[i][1][0], classes[j][1][0],
                event_at_least=event_at_least, width=width
            )
            comparisons[i, j] = cert
            if cert.positive_root_count == 0:
                if cert.cell_signs == (1,):
                    discarded.add(j)
                elif cert.cell_signs == (-1,):
                    discarded.add(i)
                else:
                    raise ArithmeticError("distinct event probability classes compare identically")
    survivors = tuple(i for i in range(len(classes)) if i not in discarded)
    if not survivors:
        raise ArithmeticError("no event probability finalists")

    def strict_greater(first: int, second: int) -> bool:
        cert = comparisons[min(first, second), max(first, second)]
        return (cert.positive_root_count == 0
                and cert.cell_signs == ((1,) if first < second else (-1,)))

    if any(not any(strict_greater(i, j) for i in survivors)
           for j in discarded):
        raise ArithmeticError("final event classes do not dominate every discarded class")

    finals = tuple(tuple(classes[i][1]) for i in survivors)
    finalists = tuple(classes[i][0] for i in survivors)

    def regenerate(current_width):
        roots = []
        for i in range(len(survivors)):
            for j in range(i + 1, len(survivors)):
                old_i, old_j = survivors[i], survivors[j]
                cert = (comparisons[old_i, old_j] if current_width == width
                        else certify_sequence_comparison(
                            landscape, finals[i][0], finals[j][0],
                            event_at_least=event_at_least, width=current_width
                        ))
                for root in cert.roots:
                    roots.append((
                        root.lo, root.hi, i, j,
                        root.sign_left, root.sign_right
                    ))
        roots.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        if len(roots) > max_events:
            raise ValueError("too many event probability pair-root events")
        separated = ((not roots or roots[0][0] > 0)
                     and all(roots[i - 1][1] < roots[i][0]
                             for i in range(1, len(roots))))
        return roots, separated

    for attempt in range(max_refinements):
        roots, separated = regenerate(width / 2**attempt)
        if separated:
            break
    else:
        raise ValueError("overlapping or coincident event roots: exact ordering unresolved")

    def probability(mask, t: Fraction) -> Fraction:
        e, z = polynomials(mask)
        denominator = _eval(z, t)
        if denominator <= 0:
            raise ArithmeticError("event denominator is not positive")
        return _eval(e, t) / denominator

    samples = []
    previous = Fraction(0)
    for lo, hi, *_ in roots:
        if previous >= lo:
            raise ArithmeticError("unavailable root-free probability interval")
        samples.append((previous + lo) / 2)
        previous = hi
    samples.append(previous + 1)

    owners = []
    for t in samples:
        values = tuple(probability(mask, t) for mask in finalists)
        maximum = max(values)
        winners = tuple(i for i, v in enumerate(values) if v == maximum)
        if len(winners) != 1:
            raise ArithmeticError("unisolated event-probability tie in root-free cell")
        owners.append(winners[0])

    events = []
    for index, (lo, hi, i, j, before_sign, after_sign) in enumerate(roots):
        before, after = owners[index], owners[index + 1]
        if before != after:
            if {before, after} != {i, j} or before_sign == after_sign:
                raise ArithmeticError("event global switch inconsistent with its unique pair root")
            kind = "global_switch"
        elif before in (i, j):
            if before_sign != after_sign:
                raise ArithmeticError("event sign-changing root cannot preserve winning class")
            kind = "global_tangent"
        else:
            kind = "irrelevant"
        events.append(EventPhaseTransition(
            lo, hi, i, j, kind, finals[before], finals[after]
        ))

    return EventProbabilityEnvelope(
        n, h_count, event_at_least, len(masks), len(skyline),
        len(classes), finals, tuple(finals[i] for i in owners),
        tuple(events)
    )
