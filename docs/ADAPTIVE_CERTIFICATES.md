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
frontier again. **A newly refined frontier can give a looser upper bound
than its parent frontier** because it may discard useful parent-level
constraints. We therefore retain every earlier checked frontier and
intersect their *individually sound* event, complement, partition and
probability mass bounds. This intersection produces a nested, sound
history even when a single frontier-only estimate temporarily widens.
The audit trail is replayable using verify_refinement_run; stored bounds
are not trusted. A valid refinement never invalidates any enclosure.
No assumption about the quality of the search policy is required.

## Zero-completion certificates

An anchored self-avoiding prefix can reach a vertex with no unoccupied
adjacent lattice sites before all residues have been placed. It then has
*zero valid completions*. Refinement must not simply drop this prefix:
doing so would fail the independently checked coverage condition.

The initial frontier generator also records trapped branches encountered
*before* the requested frontier depth; omitting them would invalidate
coverage. A terminal frontier node of kind dead asserts that no legal immediate
successor exists. Both reference verifiers independently recompute this
condition, rejecting fake claims. Such a prefix contributes precisely
zero partition mass and can remain covered by the frontier.

## Reproduce

    python -m pip install -e '.[test]'
    python -m pytest
    python experiments/compare_adaptive.py > adaptive_comparison.csv
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

## Benchmark interpretation

The CSV experiment reports exact rational enclosure widths and elapsed
wall time for equal numbers of frontier expansions, across three fixed HP
sequences and three search policies. Exact enumeration supplies only the
post-hoc correctness oracle. Wall time explicitly includes policy scoring
and full audit replay. Equal expansion budgets are **not** equal CPU budgets;
report both before drawing performance conclusions. This does not claim
comparison to published lattice-polymer methods.
