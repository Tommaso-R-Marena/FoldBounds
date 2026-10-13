# Aristotle formal library integration — staging

The uploaded Aristotle run `8becb5ac-4678-49a8-8c98-84d3faa3f7b7`
contains a reported clean Lean 4.28.0/Mathlib 4.28.0 library with
99 public declarations and 34 explicit axiom audits.

**Import status:** The full source archive was reviewed locally and is
prepared as a 21-file overlay ZIP. This staging PR adds CI and trust-boundary
documentation but **does not yet contain the uploaded Lean sources**.
Do not merge or claim independent verification until the source files and
Lake manifest are added and CI passes.

## Proof scope

The delivered Lean library models a finite conformation space, real energies,
Boltzmann weights, partition functions, probability intervals, free energy,
symmetry, optimality of coupled bounds relative to certificate data, and
rational interval arithmetic. It does **not** formalize the Python lattice
enumerator, geometry, max-flow, transfer DP, witness checker, parsing, or
an executable refinement theorem.

The extremal-landscape results show optimality **given only the specified
certificate information**. They do not establish an optimal protein-folding
algorithm or scientifically accurate protein dynamics.

## Requirements to complete the PR

1. Import `FoldBounds.lean`, all sources in `FoldBounds/`,
   `lean-toolchain`, `lakefile.toml` and `lake-manifest.json`.
2. Preserve the existing Python README and existing documentation. Import
   Aristotle's README as `docs/ARISTOTLE_TECHNICAL_REPORT.md` and the
   theorem list as `docs/THEOREMS.md`.
3. Run `lake build` under CI with Mathlib binary caches, and inspect
   the `FoldBounds/AxiomAudit.lean` outputs. The reviewing container
   has no Lean binaries, so the report cannot be marked reproduced yet.
4. Check that no claims about the executable Python checker are elevated
   into theorems without proving a refinement bridge.
5. Review scientific originality separately; standard thermodynamic
   inequalities and formal correct-by-construction infrastructure
   do not by themselves establish publication novelty.

The complete proposed overlay is available in this ChatGPT conversation
as `FoldBounds_Aristotle_Integration.zip`. If importing using local git:
unzip its contents **into the FoldBounds repository root**, then commit
the new Lean source files to this branch.
