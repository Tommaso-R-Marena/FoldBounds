# Exact all-H ground-state degeneracy and finite-temperature guarantees

**Scope:** the anchored 2D square-lattice HP model, with all residues H, at a
perfect-square sequence length n=m² (m=2,...,8). This is a Python-checked,
exact-integer research result, **not a Lean-kernel-verified Python program**
or a prediction for physical proteins. Counting Hamiltonian paths is
established mathematics, including published [OEIS A096969](https://oeis.org/A096969).

## Exact result at n=64

The previous compact-tile construction in PR #65 counts 105,840
maximum-contact 64-site folds. Connected-frontier transfer instead counts
**every** such fold, including every permitted placement of the square:

    [q^49] Z_64(q) = 1,511,656,642.

There are 1,511,656,642 canonical-first-East maximum-contact walks, **14,282.47
times as many** as the compact tile subset.

### Why all maximum-contact walks fill a square

For an n-site anchored SAW let S be its occupied sites, P its edge perimeter,
and A its number of adjacent occupied pairs. Four incident edges per site imply

    4n = 2A + P.

The n-1 backbone bonds are occupied adjacency edges, hence every all-H
contact count is

    K = A-(n-1) = n+1-P/2.

If the occupied set spans w columns and h rows, connectedness gives
P >= 2(w+h), and its cardinality gives w*h >= n. At n=m²,
P >= 4m and thus K <= (m-1)². Equality forces w+h=2m and
w*h>=m². Arithmetic-geometric mean forces w=h=m. Since precisely m²
occupied sites lie inside this m×m box, the box is full.

Conversely, every Hamiltonian path through the filled m×m box contains
exactly (m-1)² non-backbone occupied contacts. Therefore the coefficient
of q^((m-1)²) in the **full unrestricted** partition polynomial is exactly
the number of canonical first-East Hamiltonian paths through m×m squares.

The rooted walk starts at (0,0) and continues to (1,0). The position of this
root inside the square determines the square's absolute placement uniquely.
Different starting positions cannot duplicate a normalized conformation.
An m×m grid supports m(m-1) legal first-East root positions. The transfer
algorithm counts every allowed final endpoint and sums over all such roots.

### Exact connectivity transfer

Process vertices row-major, maintaining labels for connected components
meeting the frontier half-edges, a left half-edge, and a Boolean recording
whether the second degree-one endpoint has been selected. At each vertex
choose right and down half-edges subject to degree 1 at the first-East root
and unique terminal endpoint, degree 2 elsewhere. Joining two legs from
the same component would close a cycle, and is prohibited. A component
losing its final frontier leg before the last vertex cannot be connected
to the remaining unvisited vertices, and is prohibited. Accept only a
single connected, cycle-free spanning component after the final vertex.
Canonicalizing component labels merges equivalent future states.

Every Hamiltonian path gives one accepted edge-choice trace, and every
accepted trace is one Hamiltonian path. Reflection in the horizontal
midline preserves the East root bond and halves the number of root
positions that need explicit transfer runs.

Published OEIS A096969 reports 6,046,626,568 *directed* Hamiltonian walks
on the 8×8 square with unrestricted initial direction. Rotation by 90°
partitions them equally among East, North, West and South initial steps:
6,046,626,568/4 = 1,511,656,642. The DP reproduces this number, and
independent small-grid DFS tests check the actual transitions.

## Finite-q population guarantee

Let C=1,511,656,642, Kmax=49. Use the previously verified finite-block
self-avoiding-walk count certificate to bound the number of *all*
canonical length-64 conformations by N. The current default uses at most
11-step exact rooted blocks:

    N = 5,231,553,077,432,424,492,299,000,247.

For event A={K>=t}, every non-event has at most t-1 contacts. Every one
of the C known ground states is in A. If any conformation has an arbitrary
independent multiplicative energy-weight correction a_s in [l,u] with
0<l<=u, then **without assuming independence of weights**:

    P_q(A) >= l*C*q^49 / (l*C*q^49 + u*N*q^(t-1)).

For all-H n=64, q=10, t=30, the exact-rational bound is:

* l=u=1: P(K>=30) >= 0.9665495769777...
* l=1/2, u=2: P(K>=30) >= 0.8784009290029...

These are **certified conservative lower probabilities**, not estimates of
the actual ensemble probabilities. For every real q'>=10 the same numeric
lower bound holds (possibly with different admissible state corrections):
the ratio of the unfavorable upper mass to favorable lower mass is
proportional to q'^(t-1-Kmax), a nonincreasing function for t<=Kmax.

## Exact disjoint union with the PR #56 partially directed family

The existing partially directed all-H n=64 contact polynomial contains
exactly TWO conformations with K=49. Every such path must already belong
to the complete ground-state family. Thus the **exact union** of the
partially directed subset and the entire ground-state subset has weighted mass

    Z_union(q) = Z_PDSAW(q) + (1,511,656,642-2)*q^49.

Subtracting the two duplicates prevents overlap. At q=10 this improves
the directed-only partition lower bound by **1,128,589.81x**; compared
with the earlier directed + 105,840 compact-tile lower bound in #65,
the improved lower bound is **14,103.99x** larger.

The union ratio is nondecreasing in q>0 because the added monomial has
the maximum possible exponent of the directed polynomial.

## Limits, assurance and prior art

* General HP sequences with P residues are **not** covered by the
  ground-state polynomial theorem. Width >8 is not currently supported.
* The checker replays the *same Python recurrence* as its producer.
  Exhaustive small-grid tests and independently published OEIS counts are
  distinct checks, but **not a machine-checked Lean proof**.
* The Hamiltonian transfer mechanism and OEIS values are prior art.
  We do not claim discovery priority for the combinatorial count.
* Numerical gains are bounds for a toy 2D lattice model, not improved
  experimental protein-folding accuracy or calibrated temperature physics.

Reproduce:

    python -m pytest tests/test_ground_state_transfer.py -q
    FOLDBOUNDS_LONG_TESTS=1 python -m pytest tests/test_ground_state_transfer.py -q
    python experiments/ground_state_majority_demo.py --full
