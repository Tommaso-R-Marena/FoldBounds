# Three-dimensional cubic tile families: nonenumerative population certificates

**Status:** Exact-integer/rational Python research prototype, with independent
small exhaustive oracles. This is NOT a Lean kernel verification, a prediction
of natural protein structure, a new theory of protein collapse, or a discovery
of a historically novel lattice transfer matrix.

## What is new *in this FoldBounds implementation*

Earlier FoldBounds strong-coupling work used a **single fully verified
conformation** to lower bound the event population. Here a small reusable
cubic motif is independently validated, and an exact transfer computation
counts **exponentially many distinct globally self-avoiding Hamiltonian
conformations** using only finitely many local states per tile. This creates a
much stronger lower-mass certificate, and—for the first time in this
independent workstream—permits **sequence-specific H/P mixtures in 3D**.

This is a proof-of-concept method for a simplified physical model. Prior art
already contains 3D lattice HP optimization and enumeration (Schiemann,
Bachmann & Janke 2005, https://arxiv.org/abs/cond-mat/0501401), as well as
transfer-matrix enumerations of compact cubic lattice conformations (for
example, Schram et al., 2013,
https://physics-of-life.tu-dresden.de/fileadmin/core-group-files/schiessel/publications/jpa13b.pdf).
The contribution claimed here is a **reproducible certificate composition
and adversarially robust bound**, not invention of a new mathematical
research field or evidence about experimental folding.

## Construction

A cube has vertices x,y,z in {0,1}. The checker exhaustively generates
the *entire local library* of exactly **144 ordered Hamiltonian paths** that
visit each cube vertex once (and six have canonical first bond
(0,0,0)->(1,0,0)). It rejects the possibility of trusting an externally
supplied motif count.

Place cubes on a macro grid (nx,ny,nz), at offset
(2*ix,2*iy,2*iz), with **8 nx ny nz** chain sites in total. A deterministic
snake visits each macro cube once, moving by one macro coordinate at a
time. Within each small cube, choose any Hamiltonian path from the library
whose endpoints can be linked to the preceding/succeeding cube by a single
lattice bond. Every resulting global path:

1. Starts at (0,0,0)->(1,0,0).
2. Visits each of the **disjoint** 8-site macro-cube point sets once.
3. Uses exactly one unit bridge between consecutive macro cubes.
4. Is therefore a valid full 3D self-avoiding walk.

Different local-path selections give different global coordinate sequences,
so their weights can be added **without overlaps**. A dynamic programming
calculation over local endpoint states computes the EXACT number C of these
global walks, without enumerating the C walks. A parallel min-plus dynamic
program calculates a lower bound on intratile H-H contact counts for
**every** walk in the family. Intratile contacts are automatically
nonconsecutive H-H contacts in the full chain. Intertile contacts are
deliberately omitted, which only weakens the lower bound.

### An independent 3D geometric lower bound for H/P mixtures

All family walks occupy *every site* of the same rectangular solid with
side lengths sx=2nx, sy=2ny, sz=2nz. The lattice adjacency edge count is

    G=(sx-1)*sy*sz + sx*(sy-1)*sz + sx*sy*(sz-1).

The chain uses exactly n-1 of these edges, so every all-H member has

    K_allH = G-(n-1)

nonconsecutive contacts, independent of its path shape.

In the HP model, every internal **P** residue can make at most four
noncovalent contacts on the six-neighbor cubic lattice (two sites occupied
by covalent chain neighbors). Each endpoint P can make at most five.
Consequently every family member, with p total P residues, has

    K_HH >= max(0, K_allH - 4*p - P[first] - P[last]).

This is a genuinely sequence-dependent, guaranteed lower bound, even for
adversarial positions of every P. The checker takes the MAX of this bound
and its independently checked min-plus intratile-contact lower bound.

## A global upper bound that does not assume compactness

At any self-avoiding-walk prefix with a fixed first bond, old occupied sites
only REMOVE available continuations. By enumerating just the first 0..8
free-space steps on the cubic lattice and using a dynamic program to pick
the tightest product of these finite-step counts, the checker obtains a
rigorous universal cap W(n) on the number of complete 3D anchored walks.
Unlike C, W(n) includes noncompact and non-Hamiltonian conformations.

The event is A={K_HH>=t}. Each complement conformation has K_HH<=t-1.
For q>=1 its nominal Boltzmann weight is at most q^(t-1). If every
conformation's actual weight is its nominal weight times any factor
a_s in [l,u], 0<l<=u, with no independence assumptions:

    E = mass(A) >= l * C * q^Kmin       if Kmin>=t
    N = mass(not A) <= u * W(n) * q^(t-1)

and hence

    P(A) >= [l*C*q^Kmin] /
            [l*C*q^Kmin + u*W(n)*q^(t-1)].

This bound allows full **statewise adversarial** corrections (including
correlations and corrections varying from one q to another).

For Kmin>=t, the E_lower/N_upper ratio grows as q^(Kmin-t+1).
Therefore any checked strict majority at an exact rational q0 is a
*uniform conditional theorem for every real q>=q0*; these intervals
do not require numeric sampling over q.

## Reproducible exact arithmetic examples

All examples use the deterministic heterogeneous synthetic sequence
P at indices 0,11,22,... and H everywhere else. For each all-parameter
claim, q is a dimensionless contact fugacity and per-state corrections
lie independently anywhere in [1/2,2], a factor-four worst-case relative
weight shift. **These sequences are not biologically derived.**

| Occupied 3D box | Chain size n | Polar sites | Family size C | Guaranteed minimum HH contacts | Event threshold | q | Rigorous P(event) lower |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6×6×6 | 216 | 20 | 76,973,654,136,380,749,302,988,800 | 244 | 54 | 5 | > 0.999999999998 |
| 8×8×8 | 512 | 47 | > 1.8×10^62 | 644 | 256 | 6 | > 0.99999999999999 |
| 12×12×12 | 1,728 | 158 | > 10^212 | 2,391 | 864 | 5 | > 1 − 4×10^-96 |

The 216-residue result is particularly illustrative. At q=5,
the **same Kmin=244** and the **same global complement cap** but with
only ONE witness would give a probability lower bound roughly 1.08×10^-14:
completely uninformative. The exponentially large verified
family moves the rigorously established lower bound above
0.999999999998. This is a **certificate-strength improvement**;
it does not demonstrate improved accuracy for predicting an actual protein.

The *smallest certified integer q* for the stated three synthetic events
is 5, 6, and 5 respectively. This is neither the actual thermodynamic
crossover nor an experimental phase-transition estimate.

## Independent audit and failure modes

Run from the repository root with Python 3.11+:

    python -m pip install -e '.[test]'
    python -m pytest -q tests/test_cubic_tiling.py
    python experiments/cubic_tiling_demo.py

The tests use a **separately implemented enumerator** of the eight-vertex
cube Hamiltonian walks, independently enumerate and score all 36 canonical
paths in a pair of adjacent cubes, independently enumerate every anchored
8-residue 3D self-avoiding walk, compare exact statewise uncertainty extrema,
and test every 2^8 H/P sequence for the geometric polar contact lower bound.
They also reject malformed dimensions, invalid sequences, nonexact q,
invalid event thresholds, reversed/zero correction factors, and empty or
incoherent constructions. The exhaustive tests are for small or restricted
families; the huge full long-chain conformational spaces are **not** enumerated.

## Trust boundary and what is still missing

1. The *Python computation* and its combinatorial correctness argument
   have **not** been refined to or checked by Lean. No successful
   pre-existing Lean workflow changes this fact.
2. This remains the **3D simple cubic HP toy model** with only H-H
   nonbonded contact attraction. It has no atomistic sterics,
   hydrogen-bond specificity, solvent dynamics, realistic side-chain
   interactions, or an experimentally calibrated binding/contact energy.
3. The event is a **coarse number of contacts**, not a guarantee of one
   native fold, a unique geometry, biologically functional folding, or
   physically correct protein kinetics.
4. The examples are very hydrophobic (roughly 91% H) and use strong q
   values. A stringent next gate is to evaluate published diverse HP
   benchmark sequences at physically grounded q, ideally with blind
   protein measurements and false-discovery controls. The present
   examples do **not** meet that gate.
5. The universal non-event cap counts many geometrically impossible
   completions; better sequence-specific structural upper bounds
   could make the result useful at weaker q. Exact lower-mass transfer
   for **intertile** contacts is another natural improvement.
6. Scientific originality requires a deliberate comparison with prior
   exact compact-walk transfer-matrix and lattice-polymer literature.
   None is asserted here.

This is a potentially useful **verification method**: enormous *families*
of valid 3D states supply independently checkable population lower bounds
with polynomial-time transfer and very small, explicit untrusted inputs.
It is **not yet a massive breakthrough in real protein folding**.
