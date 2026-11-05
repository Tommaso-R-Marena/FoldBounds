# Infinite-subsequence concentration for a heterogeneous 3D HP lattice model

**Research status:** A reproducible exact-integer/rational Python-checked
sufficient theorem for a deliberately simplified synthetic polymer.
**Not** a Lean-kernel formalization of the 3D tile model, **not** an
experimental or atomistic protein folding result, and **not** a claim
of mathematical historical priority over lattice-polymer research.

The key improvement on the finite cubic family checker is an **analytic
constructive lower bound** for the number of valid Hamiltonian walks that
holds for *arbitrarily many* 2×2×2 cubes. This removes the need to run the
finite-state transfer separately at every large chain length.

## Fixed heterogeneous model and event

Let m be a positive integer multiple of 11, at least 11. Consider a
simple-cubic lattice box of physical side 2m, exactly filled by a
self-avoiding polymer with n=8m³ residues. The HP sequence is periodic:

    sequence[i] = P   when i % 11 == 0,
                  H   otherwise,  0 <= i < n.

Since m is divisible by 11, so is n. Consequently there are exactly
p=n/11 polar residues, the first residue is P and the last is H.
In the simple cubic HP model, K_HH denotes the number of
*nonconsecutive* nearest-neighbor H-H contacts. For exact contact
fugacity q>=1, a conformation has nominal weight q^K_HH. Every
microstate can be independently multiplied by any factor a_s in a
fixed positive interval [l,u]. No probabilistic independence or
correlation assumption is made, and these corrections may depend on
the chain size or on q.

The target event is **K_HH >= n/2**, a count of contact *edges*, not a
statement that half of the residues have a contact or that a unique
native structure exists.

## Local lemma (checked by exhaustive eight-site enumeration)

The 2×2×2 cube has precisely 144 ordered Hamiltonian paths.
For EVERY pair of an even-parity start vertex and an odd-parity
terminal vertex, there are **at least four** such walks.
For the fixed canonical root (0,0,0)->(1,0,0), there are at least two
paths finishing at the odd vertex (1,1,1).
For each even-parity starting vertex, summing over the four possible
odd-parity terminal vertices gives at least 16 paths.

These are only finite facts about the eight sites of one small cube;
they are recomputed by the checker from lattice adjacency, not trusted
from an external file. Independent tests reconstruct the paths.

## Compositional family theorem

Arrange B=m³ disjoint 2×2×2 macro-cubes in a
Hamiltonian macro-snake, with adjacent cubes joined by unit bridges.
For each bridge, select a local exit of **odd** parity on the
outgoing face. The corresponding next-cube entrance is **even**
parity. Thus:

* The first tile has >=2 valid paths with canonical first edge and
  a fixed odd exit on its +x face.
* Each of the B-2 internal tiles has >=4 valid Hamiltonian paths
  between its fixed even entry and fixed odd exit.
* The last tile, with a fixed even entry but free terminal choice,
  has >=16 Hamiltonian paths.

Every selection produces a **different** globally valid self-avoiding
Hamiltonian path through the full box, so the explicit family has size

    C >= 2 * 4^(B-2) * 16 = 2 * 4^B.        (B>=2)

This exponential lower bound holds for arbitrarily large m. It does
not count all Hamiltonian walks in the box; it counts a particular
constructively guaranteed subset.

## Contact geometry, competing mass, and exact rational bound

The full (2m)^3 cube contains

    G = 3*(2m-1)*(2m)^2 = 24m³-12m²

adjacent vertex pairs. The n-site path uses n-1=8m³-1
of them as covalent chain bonds, so every full all-H
Hamiltonian path has precisely

    K_allH = 16m³-12m²+1

nonconsecutive contacts.

Each internal P residue can remove at most four of these
H-H contacts, the first endpoint P can remove at most five,
and the terminal residue is H. Therefore EVERY member
of the constructed heterogeneous family has

    K_HH >= Kmin = 16m³-12m² - 4*n/11
                       = (144/11)m³-12m².

The event threshold is t=n/2=4m³. This construction is
nonvacuous when Kmin>=t.

Let c8=308981 be the independently recomputed number of free
eight-step continuations after a fixed first bond on the simple-cubic
lattice. Every other anchored n-site SAW has at most c8^B
possibilities in total: split its n-2=8(B-1)+6 steps into B-1
eight-step blocks and one six-step block, and observe c6<=c8.
This is an *upper bound on all compact and noncompact walks*.

For the event A={K_HH>=t}, every competing non-event has K_HH<=t-1.
Using the family lower mass and the universal complement upper mass:

    E(A) >= l * (2*4^B) * q^Kmin,
    E(not A) <= u * c8^B * q^(t-1).

Thus the complement/event mass ratio satisfies

    N/E <= F(m)
       = (u/(2*l)) * (c8/4)^(m³)
         * q^[-100m³/11 + 12m² - 1],

and consequently

    P_m(A) >= 1/(1+F(m)).

The exponent is an exact integer for every eligible m (multiple
of 11). The checker evaluates this inequality using only
integer arithmetic and Python's exact Fraction; no floating-point
energy simulation or very-long-chain enumeration occurs.

## A finite gate implies an *infinite* 3D tail theorem

The eligible side progression is m=m0+11k, for every k>=0.
Let s=11. Since

    Delta_2 = (m+s)^2-m²,
    Delta_3 = (m+s)^3-m³,

one directly checks that

    Delta_2 / Delta_3 <= 1/m <= 1/m0.

Therefore the single exact-integer inequality

    c8^(11*m0) <= 4^(11*m0) * q^(100*m0-132)

implies

    F(m+11)/F(m) <= 1

for ALL eligible m>=m0, not merely for a finite grid of sizes.
The strict positive asymptotic energy/entropy gap implied by the
same gate ensures F(m)->0 as m->infinity along that progression.
The initial finite exact-rational check of 1/(1+F(m0)) therefore
certifies a uniform lower probability bound for the WHOLE
infinite progression, and the lower bound converges to 1.

## Concrete reproducible results

All cases use independent statewise corrections in **[1/2,2]**
(fourfold maximum relative preference against the event).

| q | First certified m0 | Physical cube side | n residues | P count | Guaranteed HH Kmin | Event threshold t | Verified failure-odds upper |
|---|---:|---:|---:|---:|---:|---:|---:|
| 5 | 11 | 22 | 10,648 | 968 | 15,972 | 5,324 | < 10^-900 |
| 4 | 22 | 44 | 85,184 | 7,744 | 133,584 | 42,592 | < 10^-2700 |

At q=5, **every** m=11,22,33,... is certified to satisfy
P(K_HH>=n/2)>0.99. At q=4, the same holds for
**every** m=22,33,44,... . The bound approaches 1 along either
progression.

The first certified m0 is not asserted to be the true
thermodynamic transition or the smallest mathematically possible
starting size. For example, the present monotonicity gate
correctly abstains from claiming the whole tail starting at
m0=11 with q=4. Such an abstention does not mean that the
contact event fails in the actual finite model.

## Reproduce and independently audit

    python -m pip install -e '.[test]'
    python -m pytest -q tests/test_cubic_density_tail.py
    python experiments/cubic_density_tail_demo.py

The independent tests physically construct and inspect a complete
global Hamiltonian witness for multiple small 3D cube sizes,
verify every unit bridge and exact full-box coverage, compare
the 10,648-residue asymptotic certificate against the separate
finite-state family transfer checker, test tail recurrence
monotonicity at adjacent eligible sizes, and reject impossible
or malformed inputs. The finite-state checker performs only
small-tile enumeration and transfers, not 10,648-state or
85,184-state conformational enumeration.

## Research interpretation

This upgrades the earlier 2D all-H infinite-tail result into a
**genuinely heterogeneous H/P 3D lattice** infinite-family
certificate that tolerates adversarial weight uncertainty.

However, the sequence is still about 91% H, a synthetic periodic
arrangement chosen to make the counting lower bound informative.
The threshold event is merely contact density, not biological fold
selection. The fugacities q=4,5 are not experimentally calibrated
temperatures. Existing literature already studies 3D HP lattice
models and Hamiltonian walk families; this is a verified
composition of conservative inequalities, not a proved
historically original collapse law.

**Major remaining gates:** independent mathematical code review,
Lean formalization of the local tile library and macro-family
injection, adversarial performance against published heterogeneous
3D benchmark sequences with natural H/P ratios, quantitative
comparison to existing HP transfer/branch-and-bound methods,
and blinded experimental protein validation. Do not merge or
claim real-protein predictive success based on these results alone.
