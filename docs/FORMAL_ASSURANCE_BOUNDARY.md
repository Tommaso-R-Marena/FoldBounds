# Formal assurance boundary and integration map

## Status as of 2026-10-15

FoldBounds combines (1) Aristotle's Lean 4 / Mathlib proofs of abstract
finite-conformation thermodynamics, (2) a kernel-proved **executable
finite-rational certificate checker**, (3) a Lean-native concrete 2D HP
lattice model, and (4) a separately tested Python computational prototype.

Aristotle's run 8becb5ac-4678-49a8-8c98-84d3faa3f7b7 supplied the
original abstract library. Subsequent PRs #11-#16 added the executable
kernel, coupled-probability theorem, and concrete lattice geometry.
All three are formally compiled together. None of that by itself makes
the separate Python lattice program formally verified.

## Abstract and concrete mathematical semantics

Aristotle's abstract theory defines a finite conformation type C, real
energies E, inverse temperature beta, weights exp(-beta*E), and the
associated partition function and event probabilities.

The executable rational theory independently uses exact Q weights and
an explicitly finite model. It checks membership, positivity, and
per-state upper weight caps, then provides theorems about the returned
event/complement mass and coupled-probability intervals. Accepted
certificates cannot lie about membership or cap validity in this
explicitly enumerated finite universe.

The Lean lattice model constructs conformations from lists of four
absolute directions following a canonical first bond (0,0)->(1,0).
A finite subtype requires that a decoded conformation has no repeated
sites. Its proved theorems establish exact path length, that every
successive move is a square-grid-neighbor move, and that decoding
a list of directions is injective. It defines integer H-H contact
counts, exact rational q^contacts weights, and the contact-threshold
event. Its checkHP_coupled_sound theorem specializes successful
executable rational-checker decisions to the concrete Lean HP model.

The Python reference uses canonical anchored two-dimensional HP
self-avoiding walks, does not quotient reflections, and sums q^k for
exact rational q>=1. To interpret q physically, q=exp(beta*epsilon)
requires a separate physical/analytic bridge; arbitrary rationals
cannot be silently identified with real exponential weights.

## Proved vs. tested vs. unproved

**Lean-kernel-checked:**

- Abstract positivity, probability coupling and optimality, partition
  enclosure refinement, free energy, symmetry and variational results.
- The executable rational witness/upper-cap checker and its accepted
  mass and coupled probability enclosures.
- Concrete anchored HP conformation length, nonrepetition (by type),
  cardinal-edge adjacency, and injectivity of list-direction decoding.
- Model-relative theorem linking checkHP accepted outputs to coupled
  probability bounds for the Lean-defined finite conformation universe.

**Python-tested, not Lean-proved:**

- Prefix-frontier generation/coverage, certified witness DFS, geometric
  contact matching and reachability relaxations, discrete isoperimetry,
  weighted transfer upper bounds and adaptive bound refinements.
- Independent small-case exact-oracle benchmarks; the model still has
  extremely loose large-sequence intervals in documented cases.

**Still formally missing:**

1. A bijection/surjectivity result proving every canonical anchored
   self-avoiding path is represented by the Fin-indexed Lean direction
   model, with multiplicities matching the Python convention.
2. Verified parsing/serialization connecting Python-expressed paths,
   sequences, rational q, event thresholds and certificates to Lean.
3. A proof that Python's contact scores, path enumeration and
   production certificate algorithms refine the Lean definitions.
4. A compact non-enumerative certificate verifying mass upper/lower
   claims for families of paths without scanning the entire state
   universe; the present Lean checkHP uses a full Fintype traversal.
5. Published-method benchmarks on precisely aligned Hamiltonians,
   state-space conventions, compute and error budgets.
6. Any demonstrated physical predictive power for real proteins.

## Axiom and source policy

The current declaration-driven CI audit covers 50 selected declarations
(after PR #16) and reports only propext, Classical.choice, Quot.sound.
The scanner rejects sorry, admit, new explicit axioms and implemented_by
in the FoldBounds Lean modules. This is a focused audit, not a universal
guarantee about every dependency or third-party machine/runtime layer.

## Reproduction

    lake build
    lake env lean FoldBounds/AxiomAudit.lean
    lake env lean FoldBounds/ExecutableTests.lean
    lake env lean FoldBounds/LatticeModelTests.lean
    python3 scripts/check_lean_sources.py
    python -m pip install -e '.[test]'
    python -m pytest

Lean and Python checks are independent CI workflows. Their simultaneous
success is not a machine-checked cross-language refinement proof.

See docs/LEAN_LATTICE_MODEL.md and docs/EXECUTABLE_CERTIFICATE_KERNEL.md
for the formal definitions and exact theorem scope.