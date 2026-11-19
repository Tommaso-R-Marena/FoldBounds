# Prior art and positioning (frontier kernel)

## Inside this repository

* `FrontierCheck` / `FrontierNext` (draft PRs #26/#32): kernel-proved
  non-enumerative checker with Mathlib-linked compiled runners
  (`foldcheck`, `foldcheck2`). Semantic lemmas reused by
  `FrontierKernelSound`.
* Python research stack (`src/foldbounds`, adaptive / v2 certificates):
  producers and empirical leaf families; not Lean-kernel verified.

## What is new here

* A **Mathlib-free** executable certificate kernel whose **string-level**
  acceptance function is proved sound, including the parser and multi-query
  evaluation.
* A **grid enumerator** (`enumFast`) that is definitionally the production
  path and is proved equal to the list enumerator (`enumGo`).
* An independent Python reference + oracle + adversarial mutation corpus
  with zero observed disagreements against the compiled kernel.
* Matched-budget benchmarks showing large wall-time and binary-size wins
  without accuracy loss on exact root-`X` workloads.

## External context (non-claims)

Self-avoiding walk enumeration and HP-model thermodynamics are classical;
this work does not claim novelty for the lattice model itself. The
contribution is a proof-carrying, fail-closed checker with an executable
kernel small enough to ship beside the proofs.
