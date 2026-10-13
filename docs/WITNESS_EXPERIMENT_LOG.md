# Experiment log: certified witness lower-mass ablation

**Date:** 2026-10-13
**CI:** https://github.com/Tommaso-R-Marena/FoldBounds/actions/runs/38461558201
**Reproduce:** python experiments/compare_witnesses.py
**CI verification:** 620 tests passed on Python 3.11/3.12/3.13.

## Design

All results below concern the *finite 2D square-lattice HP toy model* with
anchored first bond, not observed proteins. The same independently checked
upper-bound machinery (matching, perimeter, transfer) is used at each witness
budget. Each budget is the number of DFS nodes visited **per bounded prefix
subtree**, not a global wall-time or number of complete conformations.

A separate exhaustive enumerator computes exact values for these very
small models. The checker itself only uses concrete accepted witnesses
and sound upper bounds, and receives no oracle values.

## Representative comparisons

| Sequence (q, event) | Depth | Lookahead | Witness node budget | Certified Z interval | Exact oracle Z | Certified event probability interval |
| --- | ---: | ---: | ---: | --- | --- | --- |
| HHHHHHHH (2, >=1) | 4 | 2 | 0 | [249,1251] | 1251 | [208/443,1016/1057] |
| HHHHHHHH (2, >=1) | 4 | 2 | 16 | [1251,1251] | 1251 | [1016/1251,1016/1251] |
| HPHHPHHP (3/2, >=1) | 4 | 2 | 0 | [253/2,1349/2] | 1349/2 | [159/833,675/769] |
| HPHHPHHP (3/2, >=1) | 4 | 2 | 16 | [1349/2,1349/2] | 1349/2 | [675/1349,675/1349] |
| HHPHHHPHHH (2, >=2) | 4 | 2 | 0 | [131,8477] | 8477 | [40/4733,3784/3875] |
| HHPHHHPHHH (2, >=2) | 4 | 2 | 256 | [8477,8477] | 8477 | [3784/8477,3784/8477] |

Example runtimes (single 3.11 GitHub Actions runner; noisy sample, not a
rigorous timing comparison): HHHHHHHH at depth 4 / lookahead 2 took
19.3 ms at budget 0 versus 26.8 ms at budget 16. The 10-residue
example took 31.8 ms at budget 0 versus 129.7 ms at budget 256.

## Critical interpretation

These examples are **small enough for exhaustive enumeration**, and
bounded DFS plus recursive lookahead may in practice cover the entire
finite conformation space. Bound coincidence is mathematically
informative but NOT a demonstration of exact partition computation
without equivalent enumeration work. It is not evidence of a solution
to realistic protein-folding dynamics or of a new complexity result.

No published-method speed or scaling baseline has been matched yet.
Increased witness budgets add computational effort and need to be
compared on equal wall time/CPU use, not only equal node budgets.

The genuinely difficult next experiment is scaling the interval width
and computation cost on larger sequences where exact oracle enumeration
is impractical, and documenting cases where the bounds remain broad.
Potential alternatives include structural decomposition and
generating-function lower bounds; any claim of novelty must pass the
existing PRIOR_ART.md research gate.

## Formal-assurance limits

Python tests and exact-rational arithmetic are not Lean proofs of the
executable. Aristotle must eventually prove the abstract subset-sum
lower-bound theorem and an actual checker refinement theorem, with
walk validation, duplicate rejection and arithmetic in the theorem's
scope, before reporting an end-to-end Lean guarantee.


## Scaling diagnostic without a full oracle

The GitHub Actions run also evaluated experiments/scaling_witnesses.py,
using the same q=2 and event >=1-contact definition, a depth-2 frontier,
lookahead 0, matching/perimeter/transfer upper bounds, and per-prefix
witness budgets of 0, 64 and 256. Larger sizes were deliberately
evaluated **without** full enumeration or claimed independent exact
verification; the small-model oracle tests remain the only empirical
soundness checks.

| All-H sequence length | Witness nodes | Verified Z interval | Partition upper/lower ratio |
| ---: | ---: | --- | ---: |
| 10 | 256 | [367,70992] | ~193 |
| 12 | 256 | [359,2555712] | ~7119 |
| 14 | 256 | [351,46002816] | ~131,062 |
| 16 | 256 | [351,1656101376] | ~4,718,238 |
| 18 | 256 | [343,29809824768] | ~86,909,110 |

The probability interval widths remain near 1 for the longer
all-H sequences: these are almost vacuous enclosures. The CPU times
in the CSV and job logs are single-run, cache-sensitive measurements
and should not be interpreted as scaling laws.

**Falsification-oriented conclusion:** brute-force budgeted witnesses
are a useful small-system correctness mechanism but are NOT a
competitive large-system thermodynamic solver. A research-significant
next step must improve mass bounds without enumerating comparable
numbers of conformations, establish end-to-end Lean assurance, and
compare rigorously with published state-of-the-art methods under
matched model and work budgets.
