# Constructive monotone-family lower bounds for 2D HP thermodynamics

**Scope:** exact rational Python research code, independently tested against
an exhaustive short-model oracle, not currently proved equivalent to Lean.
Neither elementary monotone path counting nor contact generating functions
are claimed to be mathematically original.

## Mathematical construction

Fix a valid canonical prefix p of m sites of an HP sequence of n residues.
Let r=n-m be the remaining number of sites. Consider four **monotone cones**
whose continuations use only two directions each:

- East and North
- East and South
- West and North
- West and South

A valid prefix site may never be revisited. Within any chosen cone,
however, no *newly created* suffix site can ever be revisited: the
signed-coordinate potential sx*x+sy*y strictly increases by exactly one
at every suffix step.

For two suffix sites separated by at least two chain steps, that
potential differs by at least two. Since lattice adjacency changes any
signed-coordinate potential by at most one, **nonconsecutive suffix
sites cannot be in contact**. The contact count of a conformation
in a cone therefore equals

    fixed_prefix_contacts + prefix_to_suffix_contacts.

The latter is exactly known from the position and sequence index of
each newly placed suffix residue. A dynamic program with states
(position, prefix-to-suffix contact count) gives an exact integer
contact-score histogram for all valid completions *within that cone*,
without listing each path individually.

## Avoiding duplicates without enumerating completions

Two distinct cones intersect in precisely one all-cardinal-direction
ray when they share a direction. Opposite cones have empty intersection
for r>0. Each ray occurs in exactly two cones; triple intersections
are empty. Consequently, for a positive suffix length,

    weighted_union_mass
       = sum(weighted_mass_of_each_of_four_cones)
         - sum(weighted_mass_of_each_of_four_cardinal_rays).

The same identity applies separately to the event score>=threshold and
its complement. All path weights in the union are exact q^contacts,
and each represents a **genuine self-avoiding completion**.

Thus the two resulting weighted masses are valid *lower* bounds on
the event and non-event masses for the entire prefix completion set.

For r=0 the function returns the weight of that unique completed path.
Neither mass may be added to existing ray/DFS witness lower masses,
because these methods overlap. The integrated checker takes the maximum
of each independently sound event/complement lower bound instead.

The four monotone families are disjoint as a *union set* only after
subtracting overlaps. No subtraction is performed against any previous
witness family.

## Runtime and limitations

The DP state uses position and cross-contact count, with exact integer
multiplicities. The number of possible suffix positions is polynomial in
r, and cross contacts can only occur against a fixed prefix. The
restricted-family calculation therefore avoids 2^r individual
conformations and has a polynomial-size state space for a fixed prefix.

This does not imply polynomial-time partition computation for all
self-avoiding walks. The cone family is only a subset of the full space,
and its weighted mass can be exponentially smaller than the true
partition function. Further, checking all prefixes or creating the
frontier can itself require exponential work.

The Python function is neither Lean-kernel verified nor calibrated
against real protein structures. Importing it into a Lean proof would
require establishing a verified directional-cone DP recurrence and
a verified interface to the finite model.

## Reproduction

    python -m pytest tests/test_monotone.py
    python experiments/compare_monotone.py
    python -m foldbounds HHHHHHHHHHHH --q 2 --enhanced-lookahead 0 \
      --matching --transfer --monotone

The benchmark uses matched matching/perimeter/transfer bounds with
the monotone option on and off, both without DFS and with the same
finite DFS node budget. It includes an exact oracle at length 8;
longer cases report only verified software intervals, with no
claimed independent exact reference.

## Aristotle proof obligations (deferred intentionally)

1. Formalize each monotone-quadrant suffix type and prove inclusion in
   canonical self-avoiding completions after checking fixed-prefix
   obstacles.
2. Prove the potential argument excludes both suffix self-intersections
   and nonconsecutive suffix-suffix contacts.
3. Prove the DP recurrence computes integer contact-score multiplicities
   for all cone continuations exactly.
4. Prove quadrant intersections are precisely the cardinal rays for r>0,
   so inclusion-exclusion yields an exact union histogram.
5. Prove weighted subset sums are valid event and complement lower masses.
6. Prove that replacing existing potentially overlapping witness lower
   bounds with their maximum never invalidates a certificate.
7. Connect the actual Python output to trusted Lean via a verified
   parser or independently computed Lean certificate checker.

These proof obligations are not satisfied by the present test suite.
