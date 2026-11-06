# Concrete Lean 3D geometry for FoldBounds' HP parity no-go

## Scope

This isolated add-on targets a specific **previously unformalized premise**
of draft PR #60: the exact geometry of the simple cubic lattice
and a constructive lower count on the number of conformations.
It is stacked over PR #60 and does not modify its Python source or
any other active agent's research branch.

The Lean source is \`FoldBounds/CubicParityGeometry.lean\`.
The separate workflow \`.github/workflows/cubic-parity-geometry.yml\`
builds with the repository's pinned Lean 4.28/Mathlib environment,
checks exact named theorem declarations, scans for proof bypasses,
and rejects nonstandard axiom dependencies. The Python differential
oracle is \`tests/test_cubic_geometry_bridge.py\`.

## Actual integer coordinates—not an assumed graph coloring

A lattice site is a three-dimensional integer triple

    Site3 = {x : ℤ, y : ℤ, z : ℤ}.

The relation \`Adjacent p q\` means that q is obtained from p by one
of exactly six oriented coordinate displacements

    (+1,0,0), (-1,0,0), (0,+1,0),
    (0,-1,0), (0,0,+1), (0,0,-1).

These correspond to the six offsets independently specified by the
existing Python cubic-lattice contact checker; no diagonal bonds or
zero-length steps are admitted. Color a site by

    color(p) = (p.x + p.y + p.z) % 2.

The exact Lean theorem \`advance_color_flip\` proves that for ALL
integer triples p and ALL six direction choices d,

    color(advance(p,d)) = (color(p)+1)%2.

Its corollary \`adjacent_changes_color\` rules out an edge between
two sites of equal color. This is a concrete instantiation of
bipartite graph theory, not a new graph-theoretic result.

## Position-index parity for arbitrary walks

Let a walk be any function \`walk : ℕ → Site3\` that starts at
the origin and whose consecutive positions are joined by
\`Adjacent\` edges. It need **not even be assumed self-avoiding**
for the following parity property to hold:

    color(walk(i)) = i % 2   for every natural number i.

Lean proves this by induction on i using the specific cubic
step-color flip. Hence if a sequence's H residues occur only
on one parity class of sequence index, **no two of those H sites
can share a lattice edge in ANY such walk**.

This is the geometric core of the \`parity_no_go.py\` conclusion:
for an H/P HP Hamiltonian with only HH noncovalent contact
weights q, a parity-pure H sequence has K_HH=0 in EVERY
configuration, for all q>0. Thus the nominal energy is
completely uninformative about structure.

The physical assertion is restricted to the simple cubic lattice;
it says nothing about alpha helices, real protein stereochemistry,
hydrogen bonds, solvation or atomistic potentials.

## Explicit valid, distinct monotone self-avoiding walks

For any word of k symbols from the alphabet
\`{+x,+y,+z}\`, attach that word after the fixed first
bond \`(0,0,0)->(1,0,0)\`. This defines an (indexed) chain with
exactly k+2 sites, all adjacent by construction.

The Lean theorems establish:

1. Every positive-coordinate step increases x+y+z by exactly one.
2. All site positions in the resulting trace are distinct,
   hence the walk is self-avoiding.
3. The full list of coordinates uniquely determines its
   direction word. Different words give different walks.
4. There are exactly **3^k** different direction words of
   length k; the corresponding finite-set image has the
   same cardinality by the proved injectivity.

The key declaration is \`monotone_family_card\`, with a universal
natural index k rather than a long finite sample. Taking n=k+2
gives the exponential lower bound

    number of rooted, indexed, full n-site 3D SAWs >= 3^(n-2).

This does **not** enumerate the full n-site space and does not
claim that the monotone family contains every SAW.

## Reconnection to the earlier confidence inequality

PR #60 independently proved in Lean the following generic
rational-mass implication:

If a target has mass T<=u*W, and the other conformations
have mass O>=l*W*(D−1), then

    P(target) <= u / [u+l*(D−1)].

For the parity-flat HP model, all nominal conformational
weights W coincide. The present concrete geometry is the
basis for setting D=3^(n-2).

The following **integration obligations remain open** and must
not be described as Lean-verified:

- A checked conversion from serialized Python lattice coordinates
  and model sequence to the actual Lean Site3 and walk syntax.
- A fully formal construction of an n-site finite ensemble
  containing the positive family and the associated finite
  partition sum, rather than a separately supplied mass premise.
- A formal proof of the general cube rewire geometry, contact
  rescore factoring, rigid orbit handling and published cubic
  edge-isoperimetric theorem used elsewhere in PR #60.
- Real protein force-field calibration, experimental reference
  ensembles and genuinely independent scientific peer review.

The two Lean components together establish significant portions
of the abstract theorem and concrete geometric counting but do
**not** constitute an end-to-end trusted certificate receiver
for arbitrary user-provided protein models.

## How to reproduce

    python -m pip install -e '.[test]'
    python -m pytest -q tests/test_cubic_geometry_bridge.py

In GitHub Actions, run the \`cubic-parity-geometry\` job to
kernel-check all six exact named declarations and verify their
axiom reports.

The independent Python tests check all six lattice displacements
against independently computed L1 distance-one adjacency for
each integer site in [-2,2]^3. They independently enumerate
canonical 3D SAWs of lengths 2..8, verify parity at every index
and every nonbonded contact, exhaustively check short HP parity
sequences, and generate all 3^k monotone step words for k<=7.

**Scientific interpretation:** this is an assurance and formal
modeling improvement. Bipartite coloring and monotone SAW
construction are classical mathematics; it is not evidence of a
new experimental protein-folding phenomenon.
