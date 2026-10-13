# Geometric bounds, constructive witnesses and a lookahead hierarchy

**Research status:** Mathematical argument plus tested Python reference,
**not** a Lean-kernel proof and **not** a claim of research novelty.
Everything below describes the finite anchored 2D lattice HP model in
SOUNDNESS.md. It does not validate a real-protein energy function.

## Score upper bound

Let a valid self-avoiding prefix contain m residues of an n-residue H/P
chain. Let k0 be the number of already realized nonconsecutive H-H lattice
contacts. Future contacts can only arise from candidate H-H residue pairs
with odd chain separation because the square grid is bipartite.

For a pair (i,j) with i already placed and j not yet placed, any final
contact must place j at an empty lattice neighbor of i. Starting from the
current tip, j can lie no more than j-m+1 steps away, so the distance
from the tip to i must not exceed j-m+2. Failure of these necessary
conditions eliminates a candidate pair, without eliminating any actual
conformation.

Each internal chain residue has two chain bonds and at most two
noncovalent lattice contacts; chain endpoints have at most three.
For each placed H residue, subtract observed H-H contacts and
conservatively limit further contacts by unoccupied neighboring sites.
For unplaced H residues retain their general 2-or-3 capacity. Each future
H-H contact consumes two capacities. Thus the number of new contacts is
at most the minimum of:

1. The number of candidate as-yet-unformed H-H pairs satisfying the
   necessary parity, vacancy and reachability conditions.
2. Half the total remaining H contact capacity, rounded down.

Call this minimum K. Every completion has contact score <= k0+K.

A square-lattice self-avoiding walk has at most three possible next
steps, so a prefix has no more than 3^(n-m) possible completions.
With rational Boltzmann factor q>=1 and weight q^k:

    0 <= Z(prefix) <= 3^(n-m) * q^(k0+K).

The same upper mass holds for each event class, except that a class whose
score interval excludes the requested contact threshold has mass zero.

## Constructive lower bounds

Extend the current tip along each of the four lattice directions in a
straight line until the full chain length is reached. Reject any walk
that intersects an occupied lattice position. Every remaining walk is
a mathematically valid completion with exactly calculated weight.
The largest witnessed event weight provides a sound event-mass lower
bound, and the largest witnessed complement weight provides a sound
complement-mass lower bound. These are distinct classes; separate lower
bounds can be summed. We do NOT sum weights from rays without checking
whether they are distinct.

## Nested lookahead bounds

For horizon h=0, use the geometric upper and constructive lower bounds.
For horizon h>0, enumerate all valid immediate continuations and compute
their h-1 bounds. The resulting conformational sets are disjoint.
For each event or complement mass, intersect the base bound with the
sum of all child bounds:

    L_h = max(L_base, sum child L_(h-1))
    U_h = min(U_base, sum child U_(h-1)).

An independent upper bound on the total partition function is retained
and similarly intersected with summed child upper bounds. Intersections
of sound intervals are sound. The h+1 enclosure is nested in the h
enclosure by induction, and h>=n-m computes all leaves exactly.

For global event mass E and complement mass N, suppose their bounds are
[EL,EU] and [NL,NU]. The straight full-length anchored path exists, so
E+N>0. Then

    EL/(EL+NU) <= P(event) <= EU/(EU+NL),

using the precise limiting cases P=1 when NU=0 and P=0 when EU=0.
No floating-point arithmetic is involved in the returned intervals.

## Reproduce the numerical experiments

    python -m pip install -e '.[test]'
    python -m pytest
    python experiments/compare_bounds.py
    python -m foldbounds HHHHHHHH --q 2 --depth 2 --enhanced-lookahead 4

The tests exercise all H/P sequences of length 2-8 for q values
1, 3/2 and 3 and thresholds 0, 1 and 2, using separate exact enumeration
as an oracle. Other tests check nested convergence, failed coverage,
duplicates, shadowed nodes and malformed self-avoiding walks.

For HHHHHHHH at q=2, the exact partition function is 1251. The
first reference checker gives a partition upper bound 373248. With
lookahead 4, the enhanced checker gives partition interval [249,23832].
This remains a broad interval and is only an internal comparison against
the initial baseline, not a state-of-the-art comparison.

## Pending formal assurance

The principal Lean lemmas needed are: square-grid parity; valid SAW
prefix generation and disjoint coverage; at most three successors;
the stated H coordination constraint; necessary contact reachability;
exact weight of a completed path; validity of ray witnesses; and
sound interval composition and refinement. Finally, an executable
checker refinement proof, including parsing and arithmetic, is necessary
before any statement of Lean-certified executable output.

No priority or publication novelty is claimed. Published lattice-protein
counting, bounding, and transfer-matrix algorithms must be compared before
a scientific novelty claim.
