# Adaptive conformational-certificate exploration

The FoldBounds checker makes model-relative exact-rational statements about
finite 2D lattice HP conformational probabilities. A search policy only
selects which prefix to expand. *The policy is not trusted*: the complete
resulting frontier is independently validated and its bounds recomputed.
The current Python checker is not Lean kernel verified.

## Three deterministic policies

- breadth: expand the shortest, lexicographically earliest bounded prefix.
- mass: expand the subtree with the largest aggregate event/complement
  uncertainty in its computed exact-rational mass intervals.
- probability: evaluate every legal one-step expansion, selecting the one
  which produces the narrowest certified event-probability interval. Break
  ties on partition-function interval width and canonical prefix order.

The probability-greedy policy can spend substantially more CPU time per
refinement than breadth-first expansion. None of these methods is
proven optimal or asymptotically superior to published search methods.

For each refinement, the program independently checks the complete
frontier again, checks interval nesting, and records partition/probability
bounds. A valid refinement cannot invalidate either enclosure. There is
no assumption that the choices made by the untrusted policy are correct.

## Zero-completion certificates

An anchored self-avoiding prefix can reach a vertex with no unoccupied
adjacent lattice sites before all residues have been placed. It then has
*zero valid completions*. Refinement must not simply drop this prefix:
doing so would fail the independently checked coverage condition.

A terminal frontier node of kind dead asserts that no legal immediate
successor exists. Both reference verifiers independently recompute this
condition, rejecting fake claims. Such a prefix contributes precisely
zero partition mass and can remain covered by the frontier.

## Reproduce

    python -m pip install -e '.[test]'
    python -m pytest
    python -m foldbounds HHHHHHHH --q 2 --adaptive-steps 12 --adaptive-policy probability --matching

The output contains interval history and the exact number of expensive
frontier expansions, not a promise of polynomial-time performance.

## Novelty test

Research contribution is not established by being adaptive or using max
flow alone. Compare against published algorithms on matched models and
explicitly report work-normalized certified interval widths. The greedy
policy's own scoring cost must be counted, not omitted.

Future work: prove a generic theorem that an independently validated
refinement of a disjoint covering prefix frontier preserves soundness,
including dead-end prefixes, then relate the Python implementation to
a Lean-refined checker.
