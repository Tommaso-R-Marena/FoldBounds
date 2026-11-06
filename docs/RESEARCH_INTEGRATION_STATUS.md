# FoldBounds selective research integration — DRAFT

## What is integrated

Source snapshot: 6 November 2026. The integration branch starts from **Aristotle v2 PR #32**, stacked on PR #26, and incorporates **19 additional research-tip trees**. The manifest lists 20 selected tip PRs and **61 covered PRs out of 62 then-open PRs**. Obsolete/staging PR #9 was deliberately excluded because the verified Aristotle import comes from PRs #26/#27/#32 instead.

This is a **file-level Git tree synthesis**, not a Git merge of every history. Source branches and draft PRs remain intact; `main` is untouched. Changes that add new files are included preferentially. Shared main files are reconciled separately; they were not overwritten by unrelated branch snapshots.

- **Foundational model:** Aristotle v2 Lean/Mathlib thermodynamic and frontier certificates, with compiled `foldcheck` and `foldcheck2`; the original 2D HP Python code and experimental bound algorithms remain available.
- **HP / statistical mechanics:** exact contact map atlas, rational temperature transport, spectral moment bounds, all-sequence design diagrams, 3D cubic lattice geometry, ground-state counting, mutation ensembles, parity, robust conformational bounds.
- **Structural biophysics:** geometry-aware GNM, spectral response and exact graph inverse certificates, coordinate uncertainty, NMR restraint checks, GA/GB fold fingerprint studies.
- **Kinetics:** independent kinetic barrier and pivot plateau checker experiments.
- **Experimental decision support:** S669 group-conditional evaluation, conformal/assay regret and adaptive budget verifiers, signed sparse geometry.
- **Core Python CLI:** PR #28 entry-point/adaptive/lookahead changes are selected to fully expose sparse monotone, optimized shell portfolios, coupled mass polytope and disjoint shell families.
- **Formal audit:** Aristotle v2 AxiomAudit imports are preserved while the additional protein-motion imports and theorem checks from PR #84 are unioned without duplicate declarations. A new `FoldBounds/IntegrationImports.lean` targets joint compilation of additional research Lean modules.

## How to run the checks

```bash
python -m pip install -e '.[test]'
python scripts/check_integration_sources.py
python -m pytest -q
python -m foldbounds HHHH --q 2 --depth 4 --event-at-least 1 --exact-spectrum
python -m foldbounds.research_catalog --json
lake build
lake env lean FoldBounds/IntegrationImports.lean
lake env lean FoldBounds/AxiomAudit.lean
```

Specialized tests and experiment scripts remain adjacent to each research module; the dedicated PR workflows remain in `.github/workflows`. The integration workflow is an additional gate, not a replacement for independent replication of the biological benchmarks.

## Trust and validation gates

1. **Joint compilation and regression testing are required before any merge.** A successful source-level synthesis alone is not proof that imports, package extras, or combined tests pass. No successful integration-wide CI run is claimed merely by this commit.
2. **Lean-proved vs Python-checked:** theorems are model-relative and the Python producer/checker/parser and physical-model calibration are not uniformly formally refined to the Lean kernel.
3. **Protein experiment limitations:** registered weighted-ANM RNase A evaluation against 2AAS **FAILED**; later source-only GNM robustness bounds and geometry certificates do not reverse that failure. Held-out structural fingerprints and S669 budget experiments have their own statistical and physical assumptions.
4. **No one magic predictor:** the suite includes distinct mathematically specified systems, not a new end-to-end experimentally validated predictor of protein folding. The core HP CLI intentionally does not silently apply NMR, kinetic, GNM or S669 methods to HP input.
5. **Shared CI/Lean conflicts:** PR #32's Lean workflow was retained, PR #28's main Python CLI and Python workflow were selected, and PR #84's axiom audit additions were merged. Separate specialized workflows remain in place. Changes in other PR-specific baseline workflows were **not blindly overlaid**.
6. **External validation:** verify hashes, actual data provenance, exact checker traces, group leakage controls and negative-holdout protocols separately. CIs are reproducibility aids, not automatically independent experimental validation.

See `research/integration_sources.json` for exact head SHAs and dependency chains. Draft PR must remain unmerged until required checks pass and an independent reviewer verifies mathematical and experimental claims.
