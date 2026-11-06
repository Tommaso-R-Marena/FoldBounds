"""Exact global HP sequence optimization with arbitrarily many finalist classes.

For each fixed H count, coefficientwise partition-polynomial dominance and
complete Sturm root isolation certify all optimal sequences over the *entire*
real parameter continuum q>1. The complete pair-root inventory divides the
line into root-free cells with constant order, not merely sampled q values.

Overlapping algebraic roots that this checker cannot order are rejected.
Conditional on valid atlas/polynomial inputs; Python is not Lean-refined.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Literal

from .fixed_composition_dominance import dominates, sequence_from_mask
from .exact_phase import certify_sequence_comparison, _eval, _trim


@dataclass(frozen=True)
class CriticalEvent:
    """One isolated positive root of a finalist pair difference in t=q-1."""

    t_lower: Fraction
    t_upper: Fraction
    first_class: int
    second_class: int
    kind: Literal["global_switch", "global_tangent", "irrelevant"]
    before: tuple[str, ...]
    after: tuple[str, ...]


@dataclass(frozen=True)
class MultiwayEnvelope:
    length: int
    h_count: int
    population: int
    skyline_sequences: int
    skyline_classes: int
    final_classes: tuple[tuple[str, ...], ...]
    root_free_winners: tuple[tuple[str, ...], ...]
    events: tuple[CriticalEvent, ...]

    @property
    def switches(self) -> tuple[CriticalEvent, ...]:
        return tuple(e for e in self.events if e.kind == "global_switch")


def certify_multiway_envelope(
    landscape, h_count: int, *,
    width: Fraction = Fraction(1, 1_000_000),
    max_sequences: int = 20_000,
    max_classes: int = 32,
    max_events: int = 512,
    max_refinements: int = 4,
) -> MultiwayEnvelope:
    """Certify the global maximum of Z_S(q) for all real q>1.

    Every input sequence is directly covered by a final coefficientwise
    skyline member. Every discarded skyline class has a direct all-q
    strictly greater class among the final survivors. Complete pairwise
    Sturm root inventories delimit disjoint root brackets. On each open
    root-free cell, a rational sample fixes all pairwise signs for the
    *whole cell*; the only possible ordering change inside an isolated
    bracket is the one pair whose unique root is enclosed.

    Identical polynomial classes retain all their tied sequences. Shared
    or overlapping roots from distinct differences require additional
    algebraic comparison not implemented here: the checker abstains after
    bounded refinement, rather than guessing a total root order.
    """
    n = landscape.length
    if (type(n) is not int or not 2 <= n <= 20
            or type(h_count) is not int or not 0 <= h_count <= n):
        raise ValueError("invalid chain length / H composition")
    if type(width) is not Fraction or width <= 0:
        raise ValueError("root width must be a positive exact Fraction")
    if any(type(v) is not int or v < 1 for v in
           (max_sequences, max_classes, max_events, max_refinements)):
        raise ValueError("resource limits must be positive integers")
    masks = tuple(sum(1 << i for i in indices)
                  for indices in combinations(range(n), h_count))
    if len(masks) > max_sequences:
        raise ValueError("fixed-composition enumeration budget exceeded")

    cache: dict[int, tuple[int, ...]] = {}
    conformations = None

    def polynomial(mask: int) -> tuple[int, ...]:
        nonlocal conformations
        if mask not in cache:
            seq = sequence_from_mask(mask, n)
            p = _trim(tuple(landscape.polynomial(seq)))
            if not p or any(type(v) is not int or v < 0 for v in p):
                raise ValueError("invalid partition polynomial")
            if conformations is None:
                conformations = p[0]
            if p[0] != conformations or conformations <= 0:
                raise ValueError("inconsistent q=1 conformational count")
            cache[mask] = p
        return cache[mask]

    skyline: list[int] = []
    for mask in masks:
        p = polynomial(mask)
        if any(dominates(polynomial(old), p)
               and polynomial(old) != p for old in skyline):
            continue
        skyline = [old for old in skyline
                   if not (dominates(p, polynomial(old))
                           and p != polynomial(old))]
        skyline.append(mask)

    # Independently replay final coefficientwise coverage for *every* input.
    skyline_p = tuple(polynomial(s) for s in skyline)
    if not all(any(dominates(s, polynomial(m)) for s in skyline_p)
               for m in masks):
        raise ArithmeticError("final skyline fails to dominate original set")

    groups: dict[tuple[int, ...], list[str]] = {}
    for mask in skyline:
        groups.setdefault(polynomial(mask), []).append(sequence_from_mask(mask, n))
    class_polys = tuple(groups)
    classes = tuple(tuple(g) for g in groups.values())
    if len(classes) > max_classes:
        raise ValueError("too many skyline polynomial classes")

    comparisons = {}
    eliminated = set()
    for i in range(len(classes)):
        for j in range(i + 1, len(classes)):
            cert = certify_sequence_comparison(
                landscape, classes[i][0], classes[j][0], width=width
            )
            comparisons[i, j] = cert
            if cert.positive_root_count == 0:
                if cert.cell_signs == (1,):
                    eliminated.add(j)
                elif cert.cell_signs == (-1,):
                    eliminated.add(i)
                else:
                    raise ArithmeticError("duplicate/nonzero mismatch in polynomial classes")
    survivors = tuple(i for i in range(len(classes)) if i not in eliminated)
    if not survivors:
        raise ArithmeticError("no global polynomial candidates")

    def strictly_greater(i: int, j: int) -> bool:
        cert = comparisons[min(i, j), max(i, j)]
        return (cert.positive_root_count == 0
                and cert.cell_signs == ((1,) if i < j else (-1,)))

    # Recheck direct all-parameter coverage by FINAL classes, not a chain
    # of previously eliminated intermediate candidates.
    if any(not any(strictly_greater(i, j) for i in survivors)
           for j in eliminated):
        raise ArithmeticError("final classes do not cover all discarded classes")

    finals = tuple(classes[i] for i in survivors)
    polys = tuple(class_polys[i] for i in survivors)

    def regenerate(current_width: Fraction):
        roots = []
        for i in range(len(survivors)):
            for j in range(i + 1, len(survivors)):
                old_i, old_j = survivors[i], survivors[j]
                cert = (comparisons[old_i, old_j] if current_width == width
                        else certify_sequence_comparison(
                            landscape, finals[i][0], finals[j][0],
                            width=current_width
                        ))
                for root in cert.roots:
                    roots.append((
                        root.lo, root.hi, i, j,
                        root.sign_left, root.sign_right
                    ))
        roots.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
        if len(roots) > max_events:
            raise ValueError("too many pairwise positive-root events")
        disjoint = ((not roots or roots[0][0] > 0)
                    and all(roots[k - 1][1] < roots[k][0]
                            for k in range(1, len(roots))))
        return roots, disjoint

    for attempt in range(max_refinements):
        pair_roots, separated = regenerate(width / 2**attempt)
        if separated:
            break
    else:
        raise ValueError(
            "overlapping or coincident roots: cannot certify multiway ordering"
        )

    # Every open cell contains NO root of ANY surviving pair difference.
    # Consequently exact rational evaluation at one point fixes all signs
    # and the unique winning polynomial class over the *whole* open cell.
    samples = []
    previous = Fraction(0)
    for lo, hi, *_ in pair_roots:
        if previous >= lo:
            raise ArithmeticError("root-free cell unavailable")
        samples.append((previous + lo) / 2)
        previous = hi
    samples.append(previous + 1)
    owners = []
    for t in samples:
        scores = tuple(_eval(p, t) for p in polys)
        best = max(scores)
        matching = tuple(i for i, score in enumerate(scores) if score == best)
        if len(matching) != 1:
            raise ArithmeticError("unisolated polynomial equality in root-free cell")
        owners.append(matching[0])

    events = []
    for k, (lo, hi, i, j, sign_left, sign_right) in enumerate(pair_roots):
        before, after = owners[k], owners[k + 1]
        if before != after:
            if {before, after} != {i, j} or sign_left == sign_right:
                raise ArithmeticError("global switch inconsistent with isolated root")
            kind = "global_switch"
        elif before in (i, j):
            if sign_left != sign_right:
                raise ArithmeticError("sign-changing root cannot preserve pair winner")
            kind = "global_tangent"
        else:
            kind = "irrelevant"
        events.append(CriticalEvent(
            lo, hi, i, j, kind, finals[before], finals[after]
        ))
    return MultiwayEnvelope(
        n, h_count, len(masks), len(skyline), len(classes),
        finals, tuple(finals[w] for w in owners), tuple(events)
    )
