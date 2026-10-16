# Monotone-family lower-bound ablation: 2026-10-15

**Implementation:** src/foldbounds/monotone.py
**Controlled benchmark:** experiments/compare_monotone.py
**CI benchmark log:** https://github.com/Tommaso-R-Marena/FoldBounds/actions/runs/38840574473
**Status:** Python-tested exact rational **model-relative** bounds, not a
Lean proof, published comparative advantage, or realistic protein result.

## Experimental design

The Python checker received the same complete frontier at depth 2,
lookahead 0, matching/perimeter and prefix-obstacle transfer upper
bounds, and identical DFS witness node budgets (0 or 64).
The only changed feature was the new monotone-family lower-bound DP.
The benchmark reports exact interval endpoints and timings, with an
independent full-conformation oracle at n=8.

### Key results: all-H, q=2, event is >=1 H-H contact, zero DFS nodes

| Residues | Monotone DP | Partition lower | Partition upper | Upper/lower ratio |
|---:|:---:|---:|---:|---:|
| 8 | Off | 1 | 4,680 | 4,680 |
| 8 | On | 221 | 4,680 | 21.18 |
| 12 | Off | 1 | 2,555,712 | 2,555,712 |
| 12 | On | 3,581 | 2,555,712 | 713.7 |
| 16 | Off | 1 | 1,656,101,376 | 1,656,101,376 |
| 16 | On | 57,341 | 1,656,101,376 | 28,882 |
| 18 | Off | 1 | 29,809,824,768 | 29,809,824,768 |
| 18 | On | 229,373 | 29,809,824,768 | 129,962 |

At n=18, the event/complement lower masses under monotone DP were
respectively 65,536 and 163,837 (summing to 229,373). At n=8, the
independently enumerated exact partition function was 1,251, inside
the computed [221,4,680] interval.

At a fixed **64-node DFS budget**, the n=18 lower mass changed from
48 without monotone to 229,373 with monotone, with all other model
and checker settings unchanged. This shows that the new restricted
family DP counts many different valid conformations without explicitly
listing their paths.

For the 12-residue sequence HPHHPHHPHHHH at q=3/2 and event >=2,
the method improved the complement lower mass from 1 to 3,325,
but produced no positive event lower mass. It therefore did not
materially constrain the event probability in this example.

## What these results do and do not establish

1. The new method improves **constructive lower bounds** and does
   not tighten upper bounds on its own.
2. Independent exhaustive short-model checks establish empirical
   agreement with direct restricted-family enumeration. Those tests
   are not a formal mathematical or cross-language refinement proof.
3. The reported 18-residue upper/lower ratio (~130,000) remains far
   too wide for precise partition-function predictions. Probability
   intervals for some events are still close to [0,1].
4. The four monotone cones are a restricted class. Their weighted
   mass can be exponentially smaller than the full weighted HP
   ensemble even though each cone is compactly countable.
5. The benchmark's wall-time values are cache-sensitive and
   measured on one runner. Warm-cache paired rows are not an
   equal-independent-CPU comparison. No performance speedup over
   published methods is claimed.
6. Partial-enumeration and restricted-walk techniques have significant
   prior literature. The mathematical idea is not claimed original.
7. No algorithmic theorem about runtime scaling or soundness has yet
   been checked by Lean for this Python implementation.

## Research direction

The next useful gain may come from *combining* disjoint valid restricted
families or using obstacle-conditioned dynamic programming to improve
the weighted upper mass at the same time.

An Aristotle research assignment already targets the more fundamental
Lean-checked compact frontier certificate theorem. The monotone
family DP should remain a standalone Python research module until
its precise recurrence, inclusion-exclusion correction and weight
semantics are proved in Lean or checked through an independent,
kernel-verified certificate.
