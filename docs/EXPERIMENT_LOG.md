# Experiment log: first adaptive-refinement comparison

**Date:** 2026-10-13
**CI run:** https://github.com/Tommaso-R-Marena/FoldBounds/actions/runs/38457325742
**Command:** python experiments/compare_adaptive.py
**All 356 tests passed in the same workflow.**

## Protocol

A fixed 2D HP lattice model with canonical first edge, exact-rational q,
and event defined by a minimum number of H-H contacts. Every reported
interval was independently replay-checked against a fresh exact enumerator.
The search policy did not use oracle values.

The three policies use equal numbers of bounded-prefix expansion actions
(32 below), but NOT equal execution time. The greedy policy tries all
candidate expansions at every step, making it substantially more expensive.

## Results, 32 prefix expansions, matching contact bounds

Event-probability interval width (smaller is better):

| Sequence, q, event | Breadth | Mass | Probability-greedy |
| --- | ---: | ---: | ---: |
| HHHHHHHH, q=2, >=1 contact | 0.969206 | 0.939310 | 0.913276 |
| HPHHPHHP, q=3/2, >=1 contact | 0.871374 | 0.857260 | 0.854101 |
| HHPHHHPHHH, q=2, >=2 contacts | 0.999589 | 0.999492 | 0.999080 |

Selected elapsed wall times from the same GitHub Actions 3.11 runner:

| Sequence | Breadth | Mass | Probability-greedy |
| --- | ---: | ---: | ---: |
| HHHHHHHH | 171.641 ms | 251.567 ms | 3812.972 ms |
| HPHHPHHP | 126.337 ms | 174.937 ms | 2540.879 ms |
| HHPHHHPHHH | 191.168 ms | 267.993 ms | 3673.523 ms |

The script records *all* 4, 8, 16 and 32-expansion comparisons as
CSV output, including partition widths, exact oracles, and durations.

## Interpretation and falsification

The greedy policy improved certified probability width under the fixed
32-expansion budget for these three examples, but it was much slower.
The mass policy typically improved over breadth at small extra cost.
All probability intervals remained disappointingly broad. These
experiments do NOT show a superior practical algorithm for HP
thermodynamics, nor superiority over published methods.

The next research task should prioritize major improvements to the
weighted subtree mass upper bound (e.g. generating-function,
transfer-matrix or constrained contact-count relaxations), and compare
cost against existing published algorithms. Continuing to tweak search
heuristics alone is unlikely to establish an important result.

These are tiny, deliberately chosen benchmark cases on one runner.
No statistical significance or generalization claim is supported.
