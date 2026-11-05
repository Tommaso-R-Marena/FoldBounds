"""Conservative complete all-q phase checker for arbitrary exact HP geometry.

Unlike the faster skyline/Sturm filter, this independent reference oracle
retains ALL coefficientwise-undominated polynomial classes and checks
EVERY positive root of EVERY pairwise difference using rational Sturm
isolation. It works with 2D or cubic 3D exact UniversalLandscape tables.

It refuses coincident/overlapping roots, invalid input, or budget excess.
Conditional Python certificate, not a Lean-kernel-refined executable.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Literal

from .exact_phase import _eval, certify_sequence_comparison
from .fixed_composition_dominance import dominates, sequence_from_mask


@dataclass(frozen=True)
class SkylineCriticalRoot:
    lo: Fraction
    hi: Fraction
    left_class: int
    right_class: int
    kind: Literal["global_switch", "global_tangent", "irrelevant"]
    before: tuple[str, ...]
    after: tuple[str, ...]


@dataclass(frozen=True)
class CompleteSkylineEnvelope:
    length: int
    hydrophobic_count: int
    population: int
    skyline_sequences: int
    skyline_classes: tuple[tuple[str, ...], ...]
    cell_winners: tuple[tuple[str, ...], ...]
    events: tuple[SkylineCriticalRoot, ...]

    @property
    def switches(self) -> tuple[SkylineCriticalRoot, ...]:
        return tuple(e for e in self.events if e.kind == "global_switch")


def certify_complete_skyline(
    landscape, hydrophobic_count: int, *,
    width: Fraction = Fraction(1, 1_000_000),
    max_sequences: int = 20_000,
    max_classes: int = 32,
    max_events: int = 512,
    max_refinements: int = 8,
) -> CompleteSkylineEnvelope:
    """Exact global maximizers of partition Z for ALL real q>1.

    This deliberately retains all coefficientwise skyline classes rather
    than pruning classes via additional all-parameter pair dominators.
    That supplies an independent reference root-inventory for comparing
    the faster optimizers on the same 2D or 3D landscape.
    """
    n = landscape.length
    if (type(n) is not int or not 2 <= n <= 20 or type(hydrophobic_count) is not int
            or not 0 <= hydrophobic_count <= n):
        raise ValueError("invalid chain length or hydrophobic composition")
    if type(width) is not Fraction or width <= 0:
        raise ValueError("root width must be an exact positive Fraction")
    if any(type(v) is not int or v < 1 for v in
           (max_sequences, max_classes, max_events, max_refinements)):
        raise ValueError("positive integer work budgets required")

    masks = tuple(sum(1 << i for i in indices)
                  for indices in combinations(range(n), hydrophobic_count))
    if len(masks) > max_sequences:
        raise ValueError("sequence enumeration budget exceeded")
    cache: dict[int, tuple[int, ...]] = {}

    def polynomial(mask: int) -> tuple[int, ...]:
        if mask not in cache:
            p = tuple(landscape.polynomial(sequence_from_mask(mask, n)))
            if not p or any(type(c) is not int or c < 0 for c in p):
                raise ValueError("invalid nonnegative polynomial table")
            cache[mask] = p
        return cache[mask]

    skyline: list[int] = []
    for mask in masks:
        p = polynomial(mask)
        if any(dominates(polynomial(other), p)
               and polynomial(other) != p for other in skyline):
            continue
        skyline = [other for other in skyline if not (
            dominates(p, polynomial(other)) and p != polynomial(other)
        )]
        skyline.append(mask)

    # Final-state coverage: no transitive deletion history is trusted.
    skyline_vectors = tuple(polynomial(s) for s in skyline)
    if not all(any(dominates(p, polynomial(m)) for p in skyline_vectors)
               for m in masks):
        raise ArithmeticError("final skyline does not cover every original sequence")

    grouped: dict[tuple[int, ...], list[str]] = {}
    for mask in skyline:
        grouped.setdefault(polynomial(mask), []).append(sequence_from_mask(mask, n))
    vectors = tuple(grouped)
    classes = tuple(tuple(names) for names in grouped.values())
    if len(classes) > max_classes:
        raise ValueError("too many distinct skyline polynomial classes")

    def isolate(w: Fraction):
        roots = []
        for i in range(len(classes)):
            for j in range(i + 1, len(classes)):
                comparison = certify_sequence_comparison(
                    landscape, classes[i][0], classes[j][0], width=w
                )
                if comparison.positive_root_count != len(comparison.roots):
                    raise ArithmeticError("incomplete Sturm root inventory")
                for root in comparison.roots:
                    roots.append((root.lo, root.hi, i, j,
                                  root.sign_left, root.sign_right))
        roots.sort(key=lambda e: (e[0], e[1], e[2], e[3]))
        if len(roots) > max_events:
            raise ValueError("pair-root budget exceeded")
        disjoint = ((not roots or roots[0][0] > 0)
                    and all(roots[k - 1][1] < roots[k][0]
                            for k in range(1, len(roots))))
        return roots, disjoint

    for attempt in range(max_refinements):
        roots, safe = isolate(width / (1 << attempt))
        if safe:
            break
    else:
        raise ValueError("coincident or overlapping algebraic roots: abstaining")

    samples = []
    previous = Fraction(0)
    for lo, hi, *_ in roots:
        if not previous < lo:
            raise ArithmeticError("empty root-free cell")
        samples.append((previous + lo) / 2)
        previous = hi
    samples.append(previous + 1)

    owners = []
    for t in samples:
        results = tuple(_eval(p, t) for p in vectors)
        maximum = max(results)
        winners = [i for i, x in enumerate(results) if x == maximum]
        if len(winners) != 1:
            raise ArithmeticError("unisolated equality within a root-free cell")
        owners.append(winners[0])

    events = []
    for index, (lo, hi, i, j, left, right) in enumerate(roots):
        before = owners[index]
        after = owners[index + 1]
        if before != after:
            if {before, after} != {i, j} or left == right:
                raise ArithmeticError("global switch is inconsistent with isolated pair root")
            kind = "global_switch"
        elif before in (i, j):
            if left != right:
                raise ArithmeticError("sign-changing root cannot preserve the same winner")
            kind = "global_tangent"
        else:
            kind = "irrelevant"
        events.append(SkylineCriticalRoot(
            lo, hi, i, j, kind, classes[before], classes[after]
        ))
    return CompleteSkylineEnvelope(
        n, hydrophobic_count, len(masks), len(skyline), classes,
        tuple(classes[i] for i in owners), tuple(events)
    )
