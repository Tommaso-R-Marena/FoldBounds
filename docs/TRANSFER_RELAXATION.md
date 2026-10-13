# Weighted prefix-obstacle transfer relaxation

**Status:** Exact-rational upper bounds for a finite 2D HP lattice model,
with independently checked software tests. The implementation is NOT
Lean-kernel verified, physically validated, or claimed original.

## Set-theoretic relaxation

Let p be a valid anchored self-avoiding prefix of m residues of an H/P
sequence of n residues. Every legal completion avoids all previously
occupied sites, never immediately reverses a chain edge and follows the
square-lattice adjacency relation.

Relax suffix self-avoidance: suffix walks may intersect sites created
in their own suffix, but may NEVER visit fixed prefix sites, and may
NEVER immediately reverse. This creates a strict overapproximation:

    GenuineCompletions(p) subset RelaxedSuffixes(p).

Count relaxed suffixes with a transfer dynamic program whose states are
(x, y, incoming direction, cross-contact score). Cross contacts are
nonconsecutive H-H contacts from newly placed suffix residues to
already placed H residues. Integer state multiplicities give an exact
histogram of these cross-contact scores in the relaxed system.

Let k0 be fixed H-H contacts already within p, c a relaxed suffix's
cross-contact score, and Kflow the maximum feasible number of
future-future H-H contacts under the parity-bipartite b-matching
relaxation. Let Kiso denote the maximum total contact count from
prefix-based isoperimetry. Any genuine completion obeys

    futurefuture <= min(Kflow, Kiso-k0-c).

For invalid relaxed suffixes, clamp the second term below by zero; it
still gives a nonnegative overcount. For exact rational q>=1 define

    W(c) = q^(k0+c+min(Kflow,max(0,Kiso-k0-c))).

Every valid completion injects into exactly one relaxed suffix whose
weight W(c) is at least the genuine Boltzmann weight. If H_c denotes
the exact number of relaxed suffixes of score c, then

    Z(true) <= sum_c H_c * W(c).

For the event A_t defined by final contact count >=t, include H_c W(c)
in an event upper bound whenever k0+c+K(c)>=t. For the complement
include it whenever k0+c<t. The two mass upper bounds may count
some invalid relaxed suffixes twice; this is valid because they are
upper enclosures, not an exact partition.

The independent certificate checker intersects these new bounds with
existing geometric, max-flow, explicit-witness, and bounded-lookahead
enclosures. Enabling transfer computation cannot widen an interval.

## Limits and costs

The relaxation's state count is polynomial in sequence length for
fixed 2D lattice geometry, but numerical coefficients can become very
large. The DP is exact over relaxed paths, not exact over self-avoiding
conformations, and may still give loose thermodynamic bounds.

Similar transfer-matrix and exact weighted lattice-polymer methods are
established in published literature, e.g. Lee 2018:
https://doi.org/10.1016/j.cpc.2018.03.022.
No novelty or predictive-physics claim is made here.

## Reproduction

    python -m pip install -e '.[test]'
    python -m pytest
    python experiments/compare_transfer.py
    python -m foldbounds HHHHHHHH --q 2 --enhanced-lookahead 0 --matching --transfer

Tests independently enumerate small H/P walks and compare event,
complement and partition upper enclosures to exact oracle values.
Additional tests exercise the composed checker and replay audit.

## Aristotle formalization obligations

1. Formalize genuine and relaxed suffix sets; prove inclusion.
2. Prove the integer transfer recurrence computes relaxed histogram counts.
3. Prove geometric and graph upper bounds on future-future contacts.
4. Prove that every genuine completion's weight is bounded by its
   corresponding relaxed suffix weight for all q>=1.
5. Prove event/complement upper bounds and interval intersections.
6. Relate the executable checker and graph arithmetic to the Lean
   definitions without introducing trusted algorithmic assumptions.
