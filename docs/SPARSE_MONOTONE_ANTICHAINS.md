# Sparse monotone antichain certificates — exact rational lower mass

**Status:** Python mathematical research prototype, tested against independent
exhaustive oracles. **Not Lean-verified** and not a claim about real proteins.

## Structural theorem (mathematical argument)

Let C be the fixed-anchor 2D square-lattice self-avoiding walks of length n,
with exact nonnegative rational weights w(c)=q^k(c), q>=1. For any actual
subset S of C, sum of w over S intersect event A is a valid *lower* bound on
the exact event mass; similarly for the complement. **No coverage of C is
necessary** for a lower certificate.

A certificate specifies prefixes p_1,...,p_t of genuine SAWs, all extending
the declared starting prefix, with no duplicate and no ancestor-descendant
pair. Therefore no completed walk has two listed prefixes. This prefix-free
antichain establishes that all selected completion families are **disjoint**.

At each selected prefix p, the previously tested
monotone_family_mass(sequence,p,q,threshold) exactly computes the union of
four non-opposing cardinal-direction continuation classes, with overlapping
cardinal rays subtracted. Each monotone suffix is self-avoiding because one
signed potential (+/-x +/-y) strictly increases on every move; any prefix
intersection is rejected. Suffix/suffix nonconsecutive contacts cannot occur,
since their signed-potential separation exceeds one. The integer dynamic
program sums exactly the prefix/suffix H-H contacts. This produces exact
restricted-family event and complement weighted masses without listing all
complete paths. Disjoint selected prefixes permit summing these masses.

## New recurrence: SKIP / STOP / EXPAND

At each SAW prefix, the untrusted optimizer has three options:

- SKIP: omit this entire subtree, cost 0, certified lower mass 0.
- STOP: select this prefix's exact monotone union, cost 1.
- EXPAND: allocate the remaining budget among distinct child subtrees.

Exact rational tree-knapsack convolution combines all child budget choices.
Structural induction shows that this recurrence maximizes the stated
event (or complement) lower mass among all **sparse prefix-free monotone
certificates** with at most B selected leaves and maximum core length m.
This is **not** an optimality claim about all protein conformations or all
certificate methods. The proposal algorithm can visit exponentially many
prefixes; a 25-leaf certificate does not imply a 25-node search.

The new family contains each extremal-quadrant family used by the earlier
complete shell-portfolio algorithm at the same prefix. Its allowed antichains
also include all complete frontier antichains. Consequently, the new exact
budget optimum cannot be smaller than the previous complete-shell optimum
at equal depth and leaf budget. The Python recurrence and checker have
independent brute-force optimality and conformation-enumeration tests;
the induction has not yet been expressed as a Lean theorem.

## Trust boundary

verify_sparse_monotone checks the model identity and types, canonical anchored
geometry, descendant relationship, duplicates, and prefix shadowing, then
recalculates exact rational masses using monotone_family_mass. It NEVER
trusts proposed weights or the optimizer's claimed scores. Skipped subtrees
are allowed because the result certifies a lower mass only.

Independent global lower bounds can overlap across separately submitted
certificates. FoldBounds intersects lower-mass claims using **max**, not sum,
and preserves independently computed upper masses. Exact rational coupled
probability processing remains optional. Adaptive replay saves certificates
and checks them afresh.

## Controlled budget-matched results

Event threshold is at least two H-H lattice contacts; max core depth 9.

| Sequence | q | Leaf budget | Prior optimal complete shell | Sparse monotone |
|---|---|---:|---:|---:|
| HPHHPHHPHHHH | 3/2 | 25 | 4455/4 = 1113.75 | 98289/32 = 3071.53125 |
| HPHHPHHPHHHH | 3/2 | 50 | 7587/4 = 1896.75 | 65601/16 = 4100.0625 |
| H^12 | 2 | 25 | 2848 | 13024 |
| H^16 | 2 | 25 | 46048 | 214624 |

At 25 leaves the mixed sequence improves the event mass lower bound
approximately 2.76x, **at matched certificate-leaf count, not matched CPU**.
The sparse optimizer examines 2,326 prefixes at maximum depth 9. For the
mixed 12-mer, selected complete paths count 1,632 with exact event mass
98289/32 and complement mass 1251/2; an independent exhaustive oracle over
all 30,073 anchored walks confirms these exact scores and multiplicities.

On a local reference CPU, the optimizer took roughly 0.3 seconds and the
25-leaf independent checker roughly 0.003 seconds. Runtime varies by host,
length, budget and memory; no unrestricted polynomial runtime theorem.

## Reproduce

    python -m pytest tests/test_sparse_monotone.py tests/test_sparse_monotone_integration.py
    python experiments/compare_sparse_monotone.py
    python -m foldbounds HPHHPHHPHHHH --q 3/2 --event-at-least 2 \
      --enhanced-lookahead 0 --matching --transfer --transfer-memory 4 \
      --monotone --turn-family --sparse-leaves 25 --sparse-max-depth 9 \
      --coupled-polytope

## Remaining work

The Python code does not refine Aristotle's Lean kernel by a machine-checked
theorem; the file parser and compiled binary also retain their independent
trust boundaries. Prior-art novelty, practical large-system interval
sharpness, real-protein predictive validity, and mathematical scaling against
state-of-the-art lattice-polymer transfer-matrix methods are unresolved.
