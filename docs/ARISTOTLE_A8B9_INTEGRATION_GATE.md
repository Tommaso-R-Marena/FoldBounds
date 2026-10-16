# Aristotle a8b9ef96 integration gate

**Status:** PRE-IMPORT / SOURCE INCOMPLETE. This document records independent integration checks; it is not an endorsement or an alternative to running the Lean kernel.

## Branch and provenance

Integration target: `integration/aristotle-frontier-a8b9ef96`, created from `main` commit `e373a422eb785f5cb69163b0e114343b150abf08`. Aristotle run ID: `a8b9ef96-c289-471a-b7fc-f8a7cea9be37`.

The supplied `aristotle.tar.gz` includes the reports, old project snapshot, and recorded logs, but **omits** the new proof modules, `lean_export.py`, `scripts/foldcheck.lean`, differential scripts, and `FrontierTests.lean`. The archive contains 14 top-level Lean source files and the original 50-item `AxiomAudit.lean`, whereas the report claims 26 Lean modules and 117 audited theorems. Logs are not a substitute for sources. Do not merge until the actual sources have been imported and independently checked.

## Required proof source intake

Obtain the actual contents of these 12 new files:

- `FoldBounds/LatticeEquivalence.lean`
- `FoldBounds/LatticeContacts.lean`
- `FoldBounds/Frontier.lean`
- `FoldBounds/FrontierDegree.lean`
- `FoldBounds/FrontierFamily.lean`
- `FoldBounds/Transfer.lean`
- `FoldBounds/FrontierCheck.lean`
- `FoldBounds/FrontierSoundness.lean`
- `FoldBounds/FrontierPartition.lean`
- `FoldBounds/FrontierFlat.lean`
- `FoldBounds/FrontierCLI.lean`
- `FoldBounds/FrontierTests.lean`

Also import the updated `FoldBounds.lean`, `FoldBounds/LatticeModel.lean`, `FoldBounds/AxiomAudit.lean`, `lakefile.toml` (foldcheck target), `scripts/foldcheck.lean`, `src/foldbounds/lean_export.py`, `tests/test_lean_export.py`, `experiments/lean_differential.py`, `experiments/lean_frontier_benchmark.py`, source audit updates, and applicable report files. Inspect exact directory names before assuming these are exhaustive. Preserve all current main files, including `src/foldbounds/monotone.py`, `turn_family.py`, `transfer_bound.py`, and `docs/TERMINAL_BOX_TRANSFER.md`.

## Required build and assurance gates

- `lake build` using pinned Lean/Mathlib versions; check `lake build foldcheck`.
- `lake env lean FoldBounds/FrontierTests.lean`, `ExecutableTests.lean`, `LatticeModelTests.lean`.
- `lake env lean FoldBounds/AxiomAudit.lean` and `python scripts/check_lean_sources.py`; verify the **actual source** is auditing the intended 117 declarations and all reported axioms are standard.
- `python -m pytest` in an installed development environment: existing main previously passed a 1,260-test matrix, and all original features must remain intact.
- `python experiments/lean_differential.py`: independent Python geometry, weighted score, certificate acceptance/rejection, and endpoint comparison. Include corrupted input tests; do not claim a parser theorem.
- `python experiments/lean_frontier_benchmark.py`: compare equal inputs and resource budgets; track certificate bytes, checker wall time, mass widths, and pass/fail semantics.
- Retain `FrontierCheck.check_sound` and `checkFlat_sound` theorems as verified targets with their **exact** statements and hypotheses, not informal paraphrases. Ensure the Lean checker computes sound model-relative upper/lower envelopes without full conformation enumeration, while documenting any exponential certificate-size or verification-time regimes.

## Critical trust boundaries

- The parser and text/JSON exchange remain **unverified**. A valid theorem about `check` on parsed data does not establish correctness of the parser or claims about intended source inputs.
- The Python production algorithms (matching, isoperimetric, lookahead, finite-memory contact pruning, endpoint-conditioned geometry, two-turn witnesses, adaptive search) remain **unproved refinements** of the Lean model. Do not assert Python results are Lean-certified.
- The mathematical state space is a fixed-anchor 2D HP toy model with rational `q^k`; no real-protein calibration or physical predictive guarantee.
- The original transfer result can be non-enumerative over conformations while still incurring significant DP and certificate costs; formal polynomial-time complexity is not yet established.
- Treat publication novelty and direct performance comparisons as unresolved.

## Merge protocol

Prepare a cohesive PR from this branch to `main`. Avoid blanket overwriting Python sources with Aristotle's earlier snapshot, force pushing, recreating the branch from an old commit, or merging without all CI gates and a careful diff. Update outdated README/formal assurance statements only after verifying source and CI. 
