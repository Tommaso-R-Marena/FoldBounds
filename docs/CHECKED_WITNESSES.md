# Checked conformational witnesses: constructive partition lower bounds

**Current assurance:** exact-rational software-checked results for the finite
anchored 2D HP square-lattice toy model, NOT a Lean-executable refinement theorem
or any physical validation of real proteins. Explicit conformational witnesses
and partial enumeration are established ideas, not claimed novel.

## Why this matters

Previous upper bounds used perimeter, parity matching, and a relaxed
prefix-obstacle transfer DP. The lower bound was often the Boltzmann
weight of only one known conformation in each event class, leaving
wide probability intervals.

A new *untrusted proposer* now explores self-avoiding continuations within
a deterministic DFS node budget and proposes explicit completed paths.
A small verifier independently checks every path.

For a prefix p, let W(p) be a set of proposed complete walks.
The verifier checks that each walk:

- has the correct fixed first edge and exactly n lattice vertices,
- is a 2D square-lattice self-avoiding walk,
- extends p exactly,
- appears **once**, not multiple times.

It recomputes the HP contact score k(c) directly. For exact rational q>=1,
the separately checked event and complement lower bounds are

    L_A(p) = sum_{c in W(p), k(c)>=t} q^k(c),
    L_notA(p) = sum_{c in W(p), k(c)<t} q^k(c).

Because W(p) is a set of genuine completions with no duplicates,
each sum is <= the corresponding true conformational mass. No
probabilistic assumption or acceptance of claimed proposer scores is needed.

The existing straight-ray witness bounds may overlap with these new
witnesses. We therefore use a **maximum** of the lower bounds, never
add these potentially overlapping sums. Summing event and complement
is sound because the classes are disjoint.

Within a complete prefix frontier, different bounded prefix subtrees
represent disjoint sets of completions, so their independently verified
lower mass bounds may be summed. The existing probability ratio uses
both event/complement upper and lower masses.

## Budget and trust boundaries

The new argument witness_nodes specifies how many DFS nodes to expand
*per bounded subtree*. It bounds proposer effort, not the verifier's cost
for checking emitted walks, nor overall work if the frontier/number
of recursive lookahead nodes is large. This implementation is not a
proof of efficient folding or a polynomial-time exact solver.

The generator may stop at any point. Missing conformations cannot make
these lower bounds too large. Malformed, duplicate or wrong-prefix
walks cause certificate rejection, not silent acceptance.

The verifier always recomputes values; proposed scores, probabilities
and weights are not trusted. The implementation has not yet been
proved equivalent to Lean kernel definitions.

## Reproduce

    python -m pip install -e '.[test]'
    python -m pytest
    python -m foldbounds HHHHHHHH --q 2 --enhanced-lookahead 0 --matching --transfer --witness-nodes 200
    python experiments/compare_witnesses.py

The ablation reports exact oracle values from separate full enumeration,
plus the interval width and elapsed wall time with and without
proposer effort. The node budget is not an equal-CPU benchmark.

## Formalization targets for Aristotle

1. Formalize canonical self-avoiding walks and prefix extension.
2. Prove all accepted witnesses are valid members of the relevant
   finite model's conformation space.
3. Prove distinctness and disjoint event/complement partition.
4. Prove sums over witness subsets are lower bounds for the full
   Boltzmann partition and event masses.
5. Prove distinct prefix subtrees in the covered frontier are disjoint,
   justifying addition of their lower bounds.
6. Prove combining potentially overlapping witness schemes by maximum
   (not sum) is sound.
7. Bridge the actual executable parser, walk checker and rational
   arithmetic to Lean, explicitly documenting remaining trusted code.

**Novelty:** software robustness and faster convergence are research
objectives. These elementary lower-bound arguments are not original
mathematical theorems. Published methods remain essential benchmarks.
