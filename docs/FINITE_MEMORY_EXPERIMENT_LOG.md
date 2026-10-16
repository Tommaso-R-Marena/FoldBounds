# Finite-memory transfer ablation (2026-10-15)

**Python implementation:** src/foldbounds/transfer_bound.py
**Reproduction:** experiments/compare_memory_transfer.py
**Measured CI run:** https://github.com/Tommaso-R-Marena/FoldBounds/actions/runs/38842065194

## Question

How much does excluding bounded-history self-intersections shrink a
weighted transfer upper mass while preserving all genuine square-lattice
HP self-avoiding walks?

The previous transfer method kept only enough information to forbid
immediate reversals, allowing later return cycles in its relaxation.
The new finite-memory option retains the last K sites and rejects
new suffix steps revisiting any of them. Every genuine self-avoiding
completion remains counted, and the relaxed sets are nested.

## Method

All runs use depth-2 frontier, exact rational Boltzmann factor,
matching/perimeter score caps, transfer mass upper caps, and exact
monotone-family lower masses. No additional DFS witnesses were used.
Only recent-site memory K changes: 2 (legacy), 4, 6.

The benchmark also independently computes exact Z and event
probability for the eight-residue case, checking enclosure.

### All-H models, q=2, event threshold >= 1

| Length | K | Z lower | Z upper | Upper/lower ratio | Measured milliseconds |
|---:|---:|---:|---:|---:|---:|
| 8 | 2 | 221 | 4,680 | 21.18 | 1.37 |
| 8 | 4 | 221 | 4,344 | 19.66 | 1.55 |
| 8 | 6 | 221 | 4,344 | 19.66 | 1.75 |
| 12 | 2 | 3,581 | 2,555,712 | 713.68 | 7.99 |
| 12 | 4 | 3,581 | 2,091,712 | 584.11 | 23.54 |
| 12 | 6 | 3,581 | 1,959,104 | 547.08 | 53.06 |
| 16 | 2 | 57,341 | 1,656,101,376 | 28,882 | 23.44 |
| 16 | 4 | 57,341 | 1,028,982,272 | 17,944 | 132.07 |
| 16 | 6 | 57,341 | 895,417,856 | 15,616 | 516.61 |
| 18 | 2 | 229,373 | 29,809,824,768 | 129,962 | 38.17 |
| 18 | 4 | 229,373 | 16,212,644,864 | 70,682 | 262.88 |
| 18 | 6 | 229,373 | 13,587,510,272 | 59,238 | 1,166.54 |

The n=18 upper mass is 54.4% smaller at K=6 than K=2.
The lower mass is unchanged by this option, as expected. The upper/lower
ratio remains around 59,238: still far too wide for accurate equilibrium
prediction.

For the eight-residue all-H case, the independent exact partition
function is Z=1,251. Every interval listed encloses it; the exact event
probability is 1016/1251.

### Mixed HP case

For HPHHPHHPHHHH, q=3/2 and contact threshold >=2, the Z upper changed:

| K | Z lower | Z upper |
|---:|---:|---:|
| 2 | 3,325 | 29,111,157/64 |
| 4 | 3,325 | 23,825,907/64 |
| 6 | 3,325 | 22,315,419/64 |

The event lower mass remained zero, so a useful lower probability
enclosure was **not** established by this change in that example.

## Verification and negatives

- All Python test jobs passed on 3.11, 3.12 and 3.13: 753 tests on the
  measured run, including a full independent small-model oracle.
- The bounded-memory relaxed spectrum at K=2 agrees with the legacy
  cross-contact DP over multiple conformation prefixes.
- The test suite checks monotonicity of partition, event, complement
  mass bounds and final probability enclosures as memory increases.
- The transfer DP avoids full path enumeration only by retaining a bounded
  suffix history. State cost increases sharply with memory; at n=18 the
  measured K=6 time was around 30 times K=2, not a speed improvement.
- No exact full-space reference is used for the 12-18 residue results;
  these are model-relative tested intervals, not independent full-space
  confirmations. Exact small cases and structural inclusion justify
  pursuing the relaxation, not making broad numerical physics claims.
- Relative wall-clock timings may depend strongly on interpreter, cache
  and runner load. No matched published solver comparison is claimed.
- The implementation is not yet a Lean-kernel-proved checker and the
  algorithmic concept of partial self-avoidance in transfer methods is
  not asserted to be novel.

## Next falsification priorities

1. Formalize R_K and its nestedness in Lean.
2. Prove the DP aggregation of weighted cross-contact classes.
3. Improve the upper *energy* scoring conditional on DP suffix history,
   not merely the number of relaxed walks.
4. Compare with published exact ISAW transfer methods on all-H sequences,
   with identical anchoring and no reflection quotient.
5. Establish meaningful event probability intervals (far narrower
   than near-vacuous [0,1]) on larger mixed sequences.

Aristotle's current independent proof work is left untouched.
