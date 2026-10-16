# Finite-memory transfer upper bounds

**Status:** Tested Python research implementation, not yet a Lean-verified executable.

## Mathematical relaxation

The preexisting exact-integer weighted transfer DP counts walks avoiding
all fixed-prefix vertices and immediate reversals but still allows later
suffix intersections. These relaxed paths contain every valid self-avoiding
completion, permitting conservative upper bounds on partition mass.

The new suffix_memory=K, 2<=K<=8, forbids a proposed new site if it
appeared among the last K suffix/terminal-prefix visited sites.
Memory 2 reproduces the old nonbacktracking relaxation. K>=4 excludes
four-step square loops, and larger K forbids additional short returns.

Every true self-avoiding completion remains in every relaxed set, and
the relaxed sets are nested: R_(K+1) is a subset of R_K. Because each
retained walk contributes a nonnegative exact rational upper weight,
and each event/complement category upper mass is computed by summing
nonnegative contributions, all three upper mass bounds are nonincreasing
with memory. Lower bounds (including disjoint monotone families) are
unchanged. The nested property is checked on multiple examples in CI.

The DP state is a tuple of the last K sites and exact H-H contact score
with the fixed prefix; distinct histories leading to the same state
contribute additive **integer** multiplicities. Future/future contacts
remain bounded using the existing matching and perimeter relaxations.
No floating-point values enter the partition-mass calculations.

## Complexity and trust boundaries

Retaining K sites increases the number of distinct DP states exponentially
with K even though it avoids explicit complete-path enumeration.
The operational K<=8 limit bounds time/memory exposure. This is an
anytime hierarchy of safe relaxations, not a polynomial-time exact
algorithm for the full HP conformational problem.

The Python implementation is not Lean-kernel verified and is not a
prediction of actual protein folding. The mathematical containment
argument and exact integer recurrences must be formalized before
stronger proof-carrying assurances can be claimed.

## Reproduce

    python -m pytest tests/test_transfer_memory.py
    python experiments/compare_memory_transfer.py
    python -m foldbounds HHHHHHHHHHHHHHHHHH --q 2 \
      --enhanced-lookahead 0 --matching --transfer \
      --transfer-memory 6 --monotone

Benchmarks use matching+perimeter, transfer upper bounds, monotone-family
lower mass, identical frontier depth, rational q and event threshold.
The default memory 2 remains backward compatible. Measured runtimes
depend on cache state and are not directly comparable to published solvers
without a standardized hardware and computational budget.

## Formalization milestones

1. Construct the bounded-history relaxed-walk type as a finite model.
2. Prove inclusion of genuine completions into every R_K.
3. Prove nestedness for all K>=2.
4. Prove DP state aggregation preserves full relaxed-walk multiplicity
   and fixed-prefix H-H contact scores.
5. Prove correctness of event/complement mass bounds under exact rational
   weights, including the existing future/future score cap.
6. Give a Lean-native executable certificate checker whose accepted
   interval endpoints satisfy an end-to-end theorem.

No new mathematical novelty or real-protein accuracy is claimed.

## Observed suffix-contact refinement (experimental branch)

The transfer DP can additionally retain a *joint* integer histogram of
(prefix-to-suffix H-H contacts, guaranteed short-range suffix-to-suffix
H-H contacts). For a suffix index j, it observes contacts with earlier
suffix indices i when 2 <= j-i <= K. These are precisely the pairs that
are still visible in the bounded history while j is appended; each
unordered pair is counted once. The K=2 baseline has no observable
nonconsecutive suffix contacts, so it is backward compatible.

For every **genuine** self-avoiding completion c, let k_fixed be the
contacts wholly within the fixed prefix, k_cross those connecting prefix
to suffix, and k_local the suffix contacts visible within memory K.
Every true H-H contact belongs to only one of the three categories or
to an unobserved suffix/suffix pair. Thus

    k(c) >= k_fixed + k_cross + k_local.

The ordinary transfer upper Boltzmann weight remains unchanged:
the same (cross-conditioned) globally safe upper cap is used, so neither
the partition upper nor the event upper is inflated. But for event
A_t = {k(c) >= t}, **no genuine complement conformation** can satisfy
k_fixed + k_cross + k_local >= t. Dropping such joint-DP coefficients
from the complement upper sum is safe. The exact finite-memory DP
still counts every genuine path in the relaxed set with multiplicity.

Important: The number k_local is a lower contact-score witness **only**
for actual self-avoiding walks. Invalid relaxed walks may be handled
arbitrarily conservatively because they do not correspond to genuine
conformations. It would be UNSOUND to replace the energy upper score
with k_fixed+k_cross+k_local, which is merely a lower score. The
implementation does not do that.

This feature only *may* sharpen probability lower bounds; it does not
change partition upper caps and does not directly sharpen event upper
caps. It adds histogram dimensions and can increase runtime/memory.
The independent differential tests enumerate relaxed trajectories
without state aggregation, score complete coordinate lists from scratch,
and test enclosure of every genuine short-model conformation.
No Lean proof, large-model exact oracle, or novelty is asserted.

Reproduce with:

    python -m pytest tests/test_recent_contact_pruning.py
    python -m pytest tests/test_transfer_memory.py

The full-size performance and probability-width consequences need
controlled benchmarking before claiming a practical algorithmic gain.


## First contact-pruning ablation (2026-10-15)

**Reproduce:** `python experiments/compare_recent_contact_pruning.py`

**CI run:** https://github.com/Tommaso-R-Marena/FoldBounds/actions/runs/38844945406
(all three Python 3.11/3.12/3.13 matrix jobs passed).

The ablation changes **only** the complement membership classifier, leaving
the relaxed path population, exact rational upper weights, monotone lower
witnesses, memory K, and frontiers matched. These numbers are not a
published-method comparison or a benchmark of equal implementation runtime.

For all-H sequences, q=2, event threshold >=1, canonical root prefix:

| Length | K | Old complement upper | New complement upper | Reduction | Old probability lower | New probability lower |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 4 | 2,904 | 1,880 | 35.26% | 8/371 | 8/243 |
| 12 | 4 | 1,321,280 | 489,280 | 62.97% | 16/20,661 | 16/7,661 |
| 16 | 4 | 628,887,040 | 127,979,008 | 79.65% | 32/1,228,327 | 32/249,991 |
| 18 | 6 | 8,390,747,136 | 1,467,108,352 | 82.52% | 64/8,194,153 | 64/1,432,787 |

The number of genuine H-H contacts is not certified to be near these
bounds: the entries are conservative model-relative enclosures.
Even after refinement, the length-18 probability lower bound is about
0.0000447 and the upper bound can still be close to 1. This is therefore
a real but insufficient tightening of a previously near-vacuous interval.

For mixed sequence HPHHPHHPHHHH at q=3/2, event threshold >=2,
complement upper also decreased, but the monotone-family event lower
remained zero, leaving the certified probability lower equal to zero.
This negative result must accompany positive all-H examples.

Independent differential tests enumerate relaxed paths, re-score
complete words without using the production DP, and check the
conservative mass inequalities against genuine complete walks.
The formal Python-to-Lean bridge remains unproved. Computation cost
increases because the spectrum retains an additional score dimension;
we make no throughput or complexity improvement claim.
