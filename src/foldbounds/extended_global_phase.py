"""Complete all-q global optimum by exact skyline + Sturm partial orders.

For fixed H composition, first remove sequences dominated by nonnegative
partition polynomial coefficients in t=q-1. Then group ties and remove
whole polynomial classes whose order is strictly dominated for *all q>1*
by exact, zero-positive-root Sturm comparisons. If <=2 classes remain,
classify the global optimum on every root-free continuum cell. If more
remain, fail closed rather than infer a global optimum from sampled q.

No floating-point calculations; all outputs are finite-model conditional
Python certificates rather than Lean kernel-executed proofs.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations

from .exact_phase import PhaseCertificate, certify_sequence_comparison
from .fixed_composition_dominance import dominates, sequence_from_mask


@dataclass(frozen=True)
class ExtendedGlobalDiagram:
    length: int
    h_count: int
    population: int
    coefficientwise_survivors: tuple[str, ...]
    coefficientwise_classes: int
    globally_relevant_classes: tuple[tuple[str, ...], ...]
    # For each root-free open q-cell, all tied maximal sequences.
    cell_winners: tuple[tuple[str, ...], ...]
    crossover: PhaseCertificate | None


def certify_extended_global_diagram(
    landscape, h_count, *, width=Fraction(1, 100_000),
    max_sequences=20_000, max_classes=128,
):
    """Exact global phase diagram when <=2 Sturm-undominated classes remain.

    A direct recheck of *every* sequence ensures that the coefficient skyline
    covers the original composition class. A second direct recheck ensures
    that each discarded polynomial class is globally dominated by a final
    class (not just by an intermediate class removed during pruning).
    The result is exact for every real q>1 outside rational root brackets;
    each bracket contains exactly one equality root. No q-grid is trusted.
    """
    n = landscape.length
    if (type(n) is not int or not 2 <= n <= 20
            or type(h_count) is not int or not 0 <= h_count <= n):
        raise ValueError("invalid chain size or composition")
    if type(width) is not Fraction or width <= 0:
        raise ValueError("positive exact Fraction root width required")
    if (type(max_sequences) is not int or max_sequences < 1
            or type(max_classes) is not int or max_classes < 1):
        raise ValueError("positive integer work limits required")
    masks = tuple(sum(1 << i for i in combo)
                  for combo in combinations(range(n), h_count))
    if len(masks) > max_sequences:
        raise ValueError("composition enumeration budget exceeded")

    cache = {}

    def poly(mask):
        if mask not in cache:
            cache[mask] = tuple(
                landscape.polynomial(sequence_from_mask(mask, n))
            )
        return cache[mask]

    # Online antichain: preserve precisely the earlier undominated vectors,
    # including equal-polynomial sequences with different residue positions.
    skyline = []
    for mask in masks:
        p = poly(mask)
        if any(
            dominates(poly(other), p) and not dominates(p, poly(other))
            for other in skyline
        ):
            continue
        skyline = [
            other for other in skyline
            if not (
                dominates(p, poly(other))
                and not dominates(poly(other), p)
            )
        ]
        skyline.append(mask)

    # Independently recheck direct coverage by the FINAL skyline.
    survivors = tuple(sequence_from_mask(m, n) for m in skyline)
    survivor_vectors = tuple(poly(m) for m in skyline)
    for m in masks:
        if not any(dominates(candidate, poly(m))
                   for candidate in survivor_vectors):
            raise ArithmeticError(
                "coefficient skyline failed exhaustive coverage audit"
            )

    grouped = {}
    for m in skyline:
        grouped.setdefault(poly(m), []).append(sequence_from_mask(m, n))
    classes = tuple(tuple(names) for names in grouped.values())
    if len(classes) > max_classes:
        raise ValueError("polynomial-class comparison budget exceeded")

    # A positive signed polynomial with no roots for t>0 dominates strictly
    # throughout the full continuum q>1. Root-bearing pairs are NOT pruned.
    dominated = set()
    proof_cache = {}
    for i in range(len(classes)):
        for j in range(i + 1, len(classes)):
            proof = certify_sequence_comparison(
                landscape, classes[i][0], classes[j][0], width=width
            )
            proof_cache[(i, j)] = proof
            if proof.positive_root_count == 0:
                if proof.cell_signs == (1,):
                    dominated.add(j)
                if proof.cell_signs == (-1,):
                    dominated.add(i)

    winners = tuple(i for i in range(len(classes)) if i not in dominated)

    # Recheck all discarded classes against the FINAL surviving set,
    # rather than relying solely on transitive pruning bookkeeping.
    for j in dominated:
        if not any(
            i != j and (
                (
                    proof_cache[(i, j)].positive_root_count == 0
                    and proof_cache[(i, j)].cell_signs == (1,)
                ) if i < j else (
                    proof_cache[(j, i)].positive_root_count == 0
                    and proof_cache[(j, i)].cell_signs == (-1,)
                )
            )
            for i in winners
        ):
            raise ArithmeticError(
                "Sturm survivor set does not cover every removed class"
            )

    if not 1 <= len(winners) <= 2:
        raise ValueError(
            "cannot certify full global phase: >2 root-undominated classes"
        )

    final = tuple(classes[i] for i in winners)
    if len(final) == 1:
        cells = (final[0],)
        crossing = None
    else:
        crossing = certify_sequence_comparison(
            landscape, final[0][0], final[1][0], width=width
        )
        cells = tuple(
            final[0] if sign > 0 else final[1]
            if sign < 0 else final[0] + final[1]
            for sign in crossing.cell_signs
        )
    return ExtendedGlobalDiagram(
        n, h_count, len(masks), survivors, len(classes),
        final, cells, crossing
    )
