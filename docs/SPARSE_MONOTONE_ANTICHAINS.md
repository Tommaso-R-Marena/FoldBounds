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


## Stronger optional exact union: nonmonotone paths with at most two turns

The optional `--sparse-turn-union` mode expands each monotone restricted
family by adding genuinely **nonmonotone** self-avoiding suffix direction
words having at most two direction changes. Their complete coordinates and
contact scores are recomputed exactly. These extra words are disjoint from
the monotone family by construction: if all suffix moves lie in any one
non-opposing cardinal quadrant, the path is already included in the
monotone DP, so it is excluded from the extra enumeration.

For fixed maximum turn count 2, every suffix direction word has at most
three constant-direction blocks. Distinct positive block-length
compositions and adjacent-unequal block directions uniquely represent
all words in this restricted class. There are at most

    4 * [1 + 3*(r-1) + 9*(r-1)*(r-2)/2]

candidate words for suffix length r >= 2, hence the supplementary union
can be computed by a polynomial-size restricted enumeration. This is a
restricted-family algorithm, not global polynomial-time HP folding.

For event threshold 2, matched leaf budget 25, the stronger optional union
achieves these exact rational event lower masses:

| Sequence | q | Previous complete portfolio | Sparse monotone | Sparse monotone + nonmonotone <=2-turn |
|---|---|---:|---:|---:|
| HPHHPHHPHHHH | 3/2 | 1113.75 | 3071.53125 | **71631/16 = 4476.9375** |
| H^12 | 2 | 2848 | 13024 | **20840** |
| H^16 | 2 | 46048 | 214624 | **275300** (max core depth 8) |

For the mixed 12-mer, the 25-leaf union certificate at maximum core depth
8 selects 2,876 distinct complete SAWs, with event mass 71631/16 and
complement mass 1746. An independent oracle exhaustively enumerates every
one of the 30,073 anchored complete conformations and checks exact
restricted-membership, scores, multiplicities, and rational masses.

At this budget and sequence, maximum core depth 7 already yields the same
union event optimum. Deeper optimization does not increase it for depths 8
and 9 in this tested case, while requiring more candidate exploration.

Like the base sparse method, this adds lower masses only. The independent
checker validates the versioned mode flag, recomputes exact restricted
weights, and never accepts a claimed mass from the untrusted optimizer.
Mixed overlapping lower-bound certificates are combined using max, not sum.
This extension is still Python-tested, not Lean-formally checked.


## Large-n compression demonstration (not a full partition solution)

For an all-hydrophobic 64-residue chain with `q=2`, event threshold
one H-H contact, budget 25 leaves and max core length 6:

- Prefix nodes examined by the untrusted proposer: **109**.
- Selected restricted-family certificate leaves: **25**.
- **20,032,011,142,543,966,146** distinct full-length SAWs represented.
- Exact *restricted-family* event mass: **66,437,101,702,969,556,806**.
- Exact *restricted-family* complement mass: **2,882,303,761,517,117,437**.
- Example local reference run: approximately **2.6 seconds** to plan and
  independently recheck the restricted family using exact big-integer
  arithmetic.

Those are not the exact 64-mer total partition function, event probability,
or even a nonvacuous upper bound. They demonstrate that the restricted-family
certificate can encode and check many true conformations without enumerating
individual completed walks. Exhaustive 64-mer comparison is not feasible
here; independent complete-state oracles apply only to much shorter chains.


## High-coverage restricted certificates using 200 leaves

The leaf budget is a tunable certificate-size constraint, not a global
computational cost. In the mixed 12-mer with q=3/2, event threshold 2, maximum
core depth 8, and an event-focused 200-leaf union certificate:

- True event mass (from independent exhaustive full 12-mer enumeration): 23166.
- Restricted-family event lower mass: **421407/32 = 13168.96875**,
  approximately **56.85%** of the true event mass.
- The certificate covers **10,824 distinct valid completed walks**.
- The separate complement-optimizing 200-leaf union certificate has
  independently checked complement lower mass **19426**.
- Both certificates' complete restricted memberships and rational weights
  are checked against an independently implemented enumeration of all
  **30,073** anchored 12-mer conformations.

Unlike the previous same-budget 25-leaf comparisons, this is **not** a
matched-budget comparison with the fixed-core depth-10 shell containing
thousands of enumerated prefix seeds. Its purpose is to measure how much
of the exact event mass the compressed restricted family can capture.
The complete upper-mass algorithm is unchanged. An additional controlled
experiment combines these lower bounds with the best existing optional
endpoint-conditioned transfer bounds and reports the resulting probability
interval and its exact small-ensemble oracle.

Even 56.85% of the true event weighted mass is not a probability estimate:
to infer population, one must still bound the complementary mass above.
