# Isoperimetric upper bounds from the occupied lattice perimeter

**Status:** Elementary discrete geometry; implemented with exact integers.
The Python implementation is tested, but not Lean kernel verified.
This is a known mathematical principle, not a claim of original research.

## Geometry

Let S be the occupied square-lattice sites of a completed n-residue
self-avoiding walk, and let P(S) denote the number of lattice edges
leaving S. There are 4n site-edge incidences. If A(S) is the number
of adjacent *unordered* pairs of occupied cells, then

    4n = 2 A(S) + P(S),  A(S) = 2n - P(S)/2.

The n-1 consecutive chain bonds are among these adjacent pairs.
The number C(S) of all noncovalent grid contacts therefore satisfies

    C(S) = A(S) - (n-1) = n+1 - P(S)/2.

Since every H-H contact is a noncovalent grid contact, k_HH <= C(S).

If S has bounding-box width w and height h, S meets every row and
column between its coordinate extremes because the walk is connected.
Each occupied row contributes at least two horizontal boundary edges,
and each occupied column at least two vertical boundary edges. Hence

    P(S) >= 2(w+h).

It follows that

    k_HH <= n+1-(w+h).

For any n-cell occupied set, wh >= n. Therefore

    k_HH <= n+1-ceil(2 sqrt(n)).

This global inequality is evaluated with **integer square root**
and no floating point.

## A stronger, prefix-dependent upper bound

If the current self-avoiding prefix has bounding width w0 and height h0,
every completion must have w>=w0 and h>=h0. Together with wh>=n:

    B(n,w0,h0) = min {w+h: w,h positive integers,
                            w>=w0, h>=h0, wh>=n}.

Since a connected n-cell walk spans no more than n columns, it suffices
to search the finite set w in [w0,n] and set

    h = max(h0,ceil(n/w)).

Thus every completed conformation extending the prefix satisfies

    k_HH <= n+1-B(n,w0,h0).

The implementation intersects this inequality with the previous
parity/b-matching contact cap. This is never weaker than the old
matching cap; it can be much stronger for long, elongated prefixes.

### Important distinction

The isoperimetric cap bounds *contact count*, not the number of
conformations. It is valid for any H/P sequence in this precise
2D square-lattice geometry. It does not establish the true
Boltzmann mass, real protein physics, or a novel polynomial algorithm.

## Verification

- Independently enumerate anchored self-avoiding walks of n=2 through 10.
- Compute P(S) by independently counting boundary edges.
- Verify C(S)=n+1-P(S)/2, and confirm every complete conformation
  satisfies every prefix's proposed contact cap.
- Compare matched exact partition functions and conformational event
  probabilities with and without the new cap.
- Run full pytest on all three supported Python versions.

## Formal proof obligations for Aristotle

1. Prove the lattice incidence identity 4n=2A+P.
2. Prove a self-avoiding chain contributes n-1 distinct adjacent bonds.
3. Prove row/column boundary lower bounds using connectedness.
4. Prove prefix coordinate extrema cannot decrease in a completion.
5. Prove every complete n-cell set has wh>=n.
6. Prove the integer minimization algorithm gives a safe bound.
7. Connect the bound with max-flow relaxation and with exact-rational
   partition-function inequalities.
8. Ultimately relate the Lean definitions to the executable checker.

No step should be described as Lean verified until actual proof
files are checked in the pinned Lean toolchain.
