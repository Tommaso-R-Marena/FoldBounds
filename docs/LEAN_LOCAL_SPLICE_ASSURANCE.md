# Lean local fold-splicing: exact occupancy and family injection

**Scope:** This is a list-theoretic kernel proof of important conditions
used in FoldBounds' Python local-cube rewiring certificate, not a
new protein-folding theorem or an end-to-end Lean checker for 3D HP energy.

A valid local 2×2×2 Hamiltonian rewire replaces eight consecutive
chain-indexed lattice sites with a different Hamiltonian route through
the *same set of eight sites*, with the same entry and exit coordinates.
The Python verifier directly enumerates the cube motif library, tests
global unit edges, and recomputes pairwise contact energy changes.

The Lean module \`FoldBounds/CubicSplice.lean\` formalizes the following
generic **site-list** obligations, without depending on a special cube:

- Given fixed prefix, original patch and suffix, and a replacement
  that is a list permutation of the original patch, the entire global
  walk's visited-site multiset is unchanged.
- If the original full walk had no repeated sites (was self-avoiding),
  the replacement also has no repeated sites. This requires the
  replacement to be a true permutation; merely having the same
  endpoint coordinates is NOT enough.
- Exactly the same global occupied-site set is retained.
- Distinct local replacements remain distinct after insertion
  between a fixed prefix and suffix. No accidentally duplicated
  full walk can inflate a lower partition-mass certificate.
- For any finite collection of distinct local replacement lists,
  their complete-path image has the **same cardinality**.
- Replacement cannot change global chain length.

The proofs use only the trusted \`List.Perm\`, \`List.Nodup\`,
append-cancellation, and finite-set image lemmas of Lean/Mathlib.
Six named declarations have explicit \`#print axioms\` reports and
are checked in \`.github/workflows/cubic-parity-geometry.yml\`.
The check fails closed on proof holes or unexpected axioms.

## What this does NOT prove

An arbitrary permutation of the eight occupied sites is **not**
automatically a valid rewire: it must preserve local unit bonds,
the first and last site coordinates, and global boundary covalent
bonds. These facts are checked by the Python motif enumerator and
independent brute-force tests, but their full concrete Lean
formalization is still open.

Similarly, this module proves injectivity for **one** replacement
at a time. The larger product formula for several spatially
nonadjacent rewiring cubes depends on separate geometrical
nonadjacency and contact-energy factorization obligations not yet
fully kernel checked.

The independent Python bridge constructs the complete 24-residue
three-cube Hamiltonian target, finds each legal 2×2×2 rewire,
and checks that every candidate replacement has the same eight
occupied coordinates, global uniqueness, correct chain length,
correct unit-edge adjacency and distinct resulting full paths.
This bridges the mathematical assumptions to concrete model
examples; a tested bridge is not an end-to-end formal refinement.

## Connection to geometric parity entropy

\`FoldBounds/CubicParityGeometry.lean\` adds the **specific**
three-dimensional integer-coordinate adjacency relation and proves
both infinite-walk and **finite-length** variants of the HP parity
obstruction. The finite version is crucial because some finite
self-avoiding chains have trapped endpoints and cannot be extended
to an infinite chain. The same module proves injectivity and
self-avoidance of positive-axis step words and exact cardinality
\`3^k\` for all k, not just sampled chain lengths.

Together the two modules reduce the central trust boundary behind
PR #60's exact-rational structural-confidence calculations.
Their honest assurance status remains **partially formal**:
the supplied Python model input, full energetic scoring, and
multi-patch product are not yet proved to refine the Lean definitions.

## Reproduce

    python -m pip install -e '.[test]'
    python -m pytest -q tests/test_cubic_geometry_bridge.py

Then inspect the \`cubic-parity-geometry\` GitHub Actions job, which
compiles the two Lean modules with the project's pinned Lean 4.28
and Mathlib versions, scans source for proof bypasses, and audits
all fourteen exact named declarations.
