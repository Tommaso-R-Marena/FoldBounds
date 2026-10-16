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