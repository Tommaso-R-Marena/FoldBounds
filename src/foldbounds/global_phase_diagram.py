"""Complete global all-q partition optimum after exact coefficientwise pruning.

Currently handles at most two coefficientwise-nondominated fixed-composition
sequences. For larger sets it abstains, rather than approximating winners.
"""
from dataclasses import dataclass
from fractions import Fraction

from .fixed_composition_dominance import certify_fixed_composition_dominance
from .exact_phase import certify_sequence_comparison, PhaseCertificate


@dataclass(frozen=True)
class GlobalPhaseDiagram:
    length: int
    h_count: int
    population: int
    candidates: tuple[str, ...]
    # One winner (or tied pair) for each open root-free parameter cell q>1.
    cell_winners: tuple[tuple[str, ...], ...]
    # Exact rational bounds on sorted positive roots in t=q-1.
    comparison: PhaseCertificate | None


def certify_global_phase_diagram(landscape, h_count, *,
                                 width=Fraction(1, 10_000)):
    """Certify the complete fixed-composition all-q partition envelope.

    Because the coefficientwise order is finite, every sequence is dominated
    by some maximal element. If there are <=2 maximal elements, comparing
    these using exact Sturm root isolation determines every global optimum
    in the root-free parameter regions. Equality lies in certified rational
    root brackets. For more than two maxima, fail closed.
    """
    if type(width) is not Fraction or width <= 0:
        raise ValueError("width must be a positive exact Fraction")
    dom = certify_fixed_composition_dominance(landscape, h_count)
    candidates = dom.coefficientwise_nondominated
    # Recheck that the *final* candidates cover every sequence directly.
    # This extra O(binomial(n,k)*candidate_count) pass does not rely on the
    # correctness of the partial-order pruning algorithm's transitivity.
    from itertools import combinations
    from .fixed_composition_dominance import sequence_from_mask
    winner_vectors = [landscape.polynomial(s) for s in candidates]
    for positions in combinations(range(landscape.length), h_count):
        mask = sum(1 << i for i in positions)
        row = landscape.polynomial(sequence_from_mask(mask, landscape.length))
        if not any(all(
            (w[j] if j < len(w) else 0) >= (row[j] if j < len(row) else 0)
            for j in range(max(len(w), len(row)))
        ) for w in winner_vectors):
            raise ArithmeticError("survivor dominance does not cover every sequence")
    if len(candidates) == 1:
        return GlobalPhaseDiagram(
            landscape.length, h_count, dom.sequences_checked,
            candidates, (candidates,), None
        )
    if len(candidates) != 2:
        raise ValueError(
            "more than two coefficientwise candidates: global phase not certified"
        )
    a, b = candidates
    c = certify_sequence_comparison(landscape, a, b, width=width)
    winners = tuple(
        (a,) if sign > 0 else (b,) if sign < 0 else (a, b)
        for sign in c.cell_signs
    )
    return GlobalPhaseDiagram(
        landscape.length, h_count, dom.sequences_checked,
        candidates, winners, c
    )
