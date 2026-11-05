# 3D HP parity no-go: a provably flat folding energy landscape

**Trust boundary:** a simple-cubic *two-type HP lattice model*, not
an atomistic protein. A mathematical sequence-geometry obstruction
shows that certain synthetically chosen H/P sequences cannot make
any nearest-neighbor H-H contact in **any** lattice conformation.
This is standard bipartite graph mathematics rather than a newly
discovered physical phenomenon.

## Exact all-conformation obstruction

The simple cubic nearest-neighbor graph is bipartite, colored by

    color(x,y,z) = (x+y+z) mod 2.

Each legal chain step moves by one unit along a lattice coordinate
and therefore flips this color. Starting with the canonical root
(0,0,0), residue i must occupy lattice color i mod 2 in EVERY
self-avoiding conformation. A nearest-neighbor contact also crosses
the color bipartition.

Consequently **every H-H nearest-neighbor contact necessarily joins
an even-indexed and an odd-indexed H residue**. If a sequence has
all its H residues at even indices or all at odd indices, the
entire conformational space has zero H-H noncovalent contacts:

    for ALL valid 3D lattice conformations P:
      K_HH(P) = 0.

The condition includes all-P sequences and many nonalternating
H/P sequences; alternating HP is a particularly important example.

For the classical HP Hamiltonian with nonbonded contact weight q>0
for H-H and weight 1 for H-P and P-P:

    W_q(P) = q^K_HH(P) = 1  for EVERY P, and EVERY q>0.

The energy model is **completely structure-blind** for these
sequences. No change of q or temperature can recover a
preferred structure within the specified Hamiltonian.

## Exact exponential lower bound on the number of conformations

Fix the canonical first two sites (0,0,0)->(1,0,0). Each sequence
of the n-2 subsequent steps chosen from {+x,+y,+z} yields a
valid self-avoiding walk because x+y+z strictly increases at
every step. Distinct step words give distinct complete
indexed lattice conformations. Thus the global ensemble has
at least

    D(n) = 3^(n-2)

different states; no full enumeration is necessary.

In an uncorrected HP model all states have weight 1, hence for
any candidate target microstate

    P_q(target) = 1/N_full(n) <= 1/3^(n-2).

If actual weights are independently reweighted by arbitrary
positive microstate corrections a_P∈[l,u], the exact worst-case
target probability obeys

    P(target) <= u/[u+l*(3^(n-2)-1)].

The rigid-symmetry quotient of an indexed fold has at most 48
anchor-respecting coordinate representatives; a conservative
fold-class bound is

    P(target fold class)
      <= 48u/[48u+l*(3^(n-2)-48)], if 3^(n-2)>48.

This conclusion is valid for ANY q>0, including q→infinity.

## Concrete large synthetic counterexample

For n=216, alternating HP has exactly 108 hydrophobic and
108 polar residues, a 50/50 mixture. For ANY q>0 and arbitrary
statewise corrections in [1/2,2]:

    D(n) >= 3^214 ≈ 1.27×10^102,
    P(any exact indexed target) <= 4/(3^214+3)
                                 ≈ 3.14856×10^-102.

This is an exact-rational sufficient population **upper bound**,
not a Monte Carlo estimate or a measurement of natural protein
heterogeneity. The model cannot distinguish the conformations;
adding real direction-specific interactions or other physics
would invalidate the flat-landscape premise.

The example also illustrates why interpreting every H/P sequence
as a protein-designable sequence would be a scientific mistake.

## Independent tests and kernel-checked core

Run:

    python -m pip install -e '.[test]'
    python -m pytest -q tests/test_parity_no_go.py
    python experiments/parity_no_go_demo.py

The independent test enumerates ALL rooted 3D lattice self-avoiding
walks of lengths 2 through 8, directly checks that every spatial
nearest-neighbor contact joins opposite chain-index parities, and
tests every H/P sequence of lengths 2 through 7 that satisfies
the parity precondition. It separately enumerates all monotone
positive-axis step words for short lengths and confirms they
are valid distinct self-avoiding conformations. Additional tests
reject mixed H-parity sequences and malformed or nonexact inputs.

\`FoldBounds/ParityNoGo.lean\` checks the **abstract bipartite
graph proposition** (same-color hydrophobes cannot be adjacent).
The concrete cubic-lattice coloring, parity indexing of a walk,
and injective 3^(n-2) lower bound remain Python-checked reasoning,
NOT end-to-end Lean-proved geometry.

### Context and originality

Bipartite parity restrictions and HP degeneracy are classical and
must not be presented as new molecular physics. The value of this
component is a simple **fail-closed impossibility detector**:
before publishing a candidate fold's HP-model equilibrium
confidence, the verifier can show when the declared Hamiltonian
provides no structural information whatsoever. If H occurs on
both index parities, this particular checker **abstains** rather
than asserting folding is possible or impossible.

The model must not be conflated with directional hydrogen bonds,
solvent effects, atomic side chains, kinetics, measured equilibrium
populations, or state-of-the-art real protein predictors.
