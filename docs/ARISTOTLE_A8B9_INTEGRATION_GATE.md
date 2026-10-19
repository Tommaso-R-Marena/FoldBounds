# Aristotle a8b9ef96 integration gate

**Status:** Complete source overlay imported into a private-to-PR staging branch, pending independent integrated CI. Source availability alone is not formal verification.

## Branch and provenance

Integration target: `integration/aristotle-frontier-a8b9ef96`, created from `main` commit `e373a422eb785f5cb69163b0e114343b150abf08`. Aristotle run ID: `a8b9ef96-c289-471a-b7fc-f8a7cea9be37`.

The **first** Aristotle archive lacked the new proof modules. The **second** supplied integration bundle contains all 35 listed source/report files, with every SHA-256 matching the upstream manifest. All 12 new proof modules, source imports, the expanded 117-theorem audit, foldcheck runner, exporter and experiments have now been imported onto the dedicated staging branch. Existing Python research modules were preserved. Actual kernel and cross-language CI on the integrated GitHub tree is still required before PR #26 can be merged. Full upstream digests are in `docs/ARISTOTLE_A8B9_SOURCE_MANIFEST.md`.

## Source intake: imported on staging, awaiting integrated verification

The following 12 new files were imported verbatim from the verified Aristotle overlay:

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

The updated `FoldBounds.lean`, `FoldBounds/AxiomAudit.lean`, `lakefile.toml` and existing README were merged with `main` rather than blindly replaced. `FoldBounds/LatticeModel.lean` and `scripts/check_lean_sources.py` matched the existing branch and were preserved unchanged. The compiled runner, exporter, experiments, reports, and tests were added. Current `main` Python methods, including `monotone.py`, `turn_family.py`, `transfer_bound.py` and `docs/TERMINAL_BOX_TRANSFER.md`, were left intact. Verify these facts again from the final PR diff before merging.

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

The source staging branch will first be merged via a dedicated PR into the existing `integration/aristotle-frontier-a8b9ef96` branch, then PR #26 will remain draft against `main` until final CI and review. Do not force push, overwrite Python research modules, or merge to `main` without every substantive gate. The README and assurance map now identify both the proved frontier checker and the unproved Python integration boundaries. 
