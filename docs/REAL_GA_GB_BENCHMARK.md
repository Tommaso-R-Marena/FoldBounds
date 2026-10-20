# FoldBounds: real-protein GA95/GB95 structural analysis and experimental population benchmark

**Status: experimental validation infrastructure and necessary thermodynamic constraints — NOT an atomistic fold prediction, NOT a discovered solution to protein folding.**

## Why this matters

The experimentally characterized GA95/GB95 fold switch offers a hard real-protein benchmark. PDB **2KDL** is the three-helix, human-serum-albumin-binding GA95 protein and **2KDM** the four-beta-plus-helix, IgG-binding GB95 protein. Both are 56 amino acids and differ at exactly residues 20, 30, and 45: L20A, I30F, L45Y going from GA95 to GB95. The 2009 PNAS paper reports an experimentally observed fold-state reversal after a single substitution along a designed mutational trajectory with >90% populations of the respective folds under reported conditions. This does **not** mean that a static NMR ensemble with 20 submitted coordinate models encodes a 20-sample equilibrium distribution.

In August 2026, Felbinger et al., *Benchmarking Deep Learning Predictions of Mutation-Induced Fold Switching* (bioRxiv DOI `10.64898/2026.08.02.742283`), introduced a wider, measured NMR dataset including 60 GA/GB variants with 3-alpha, 4-beta-plus-alpha, mixed, and unfolded state populations. A reproducible prospective leaderboard should predict these unseen mutant populations using a frozen model and compare to data the model did **not** see. This PR builds only the honest input/provenance and evaluation foundations.

References:
- 2KDL: https://www.rcsb.org/structure/2KDL
- 2KDM: https://www.rcsb.org/structure/2KDM
- 2009 paper: https://doi.org/10.1073/pnas.0906408106
- 2025 aromatic interaction study: https://doi.org/10.1039/D4SC04951A
- 2026 NMR mutant benchmark: https://doi.org/10.64898/2026.08.02.742283

## Verifiable program behavior

1. A dedicated CI job retrieves **the official experimental PDB structures** from RCSB. The program checks PDB IDs, NMR metadata, precisely 20 complete 56-residue chain-A C-alpha traces per entry, every sequence identity, and the expected point substitutions. It records SHA-256 hashes of all downloaded coordinate files; no guessed coordinates are permitted.
2. For each NMR structural model, it computes a rotation-independent geometric C-alpha contact map, with separation at least 4 residues and an 8 Å distance threshold. It defines **75%-model-consensus contacts**, measures the contact-map overlap and helix/extended local C-alpha distance proxies, and reports aromatic ring-centroid proximity at positions 29, 30, 45, and 52. These are *descriptive structural metrics*; they are not a force field or energy scoring function.
3. The scorer accepts three-state numerical predictions (three-alpha, four-beta-plus-alpha, unfolded), with required model Git SHA and training-freeze metadata, and matching labels. It rejects invalid probabilities, duplication, mismatched variant IDs, and missing model provenance. Metrics: multiclass Brier score, cross-entropy and dominant-state accuracy. It does NOT claim predictions were blind merely because a metadata field was filled.
4. A separate exact-rational thermodynamic module derives sharp necessary constraints and uncertainty effects, conditional on experimentally specified population lower bounds.

## Necessary energy/odds shift for a real fold reversal

Suppose mutant `X` has an A-fold population **at least 90%**, and mutant `Y` has a B-fold population **at least 90%**, at the *same physical conditions*. Since probabilities sum to one, their relative A/B odds satisfy

    (p_A(X)/p_B(X)) / (p_A(Y)/p_B(Y)) >= 81.

This statement remains correct if **unfolded or other states** are populated, since the competing fold population is no larger than the residual probability. Under equilibrium statmech, the difference in the two-fold free-energy gap between X and Y therefore has magnitude at least

    R*T*ln(81) = 2*R*T*ln(9) ≈ 2.56 kcal/mol at T=293.15 K.

This is a necessary constraint derived from measured populations, not a computed GA/GB energy, mutation prediction, or claim of measured stability differences for specific mutants without inspecting their labels. For a path of 3 residue changes, at least one step must contribute a gap shift of absolute magnitude ≥(2.56/3) kcal/mol if the sum over that specified path accounts for the endpoint difference (triangle inequality). This is not a prediction of which site will produce the switch.

## Certification is fragile under arbitrary conformer energy errors

If a nominal state population is `p`, and **every** state's statistical weight may be multiplied independently by a factor `f∈[lo,hi]`, then its sharp worst-case real population lies within

    lo*p/(lo*p + hi*(1-p)) <= p_true <= hi*p/(hi*p + lo*(1-p)).

This is the same general statement kernel-checked abstractly in prior FoldBounds uncertainty work. For an energy error ±delta per state at T, factors are exp(±delta/RT). A nominal 95% fold population cannot be guaranteed above 90% unless the per-state energy error budget satisfies

    delta <= (RT/2) * ln( (0.95/0.05) / (0.90/0.10) ) ≈ 0.218 kcal/mol at 293.15 K.

Such *independent, adversarial per-state* errors are a worst-case model, not a validated error distribution for modern molecular force fields. Common-mode errors and correlated perturbations can cancel, and exploiting their structure is a key research direction. Report high-confidence guarantees only if error bounds are defensible and validated externally.

## Run locally

Install 

    python -m pip install -e '.[test]'

then download actual PDB files from RCSB (or let the dedicated GitHub workflow do this), and run:

    python experiments/analyze_ga_gb_nmr.py --ga /path/2KDL.pdb --gb /path/2KDM.pdb --out /path/nmr_analysis.json

NMR files must originate from public PDB data, and their hashes appear in the report. The structural metrics are insensitive to superposition and pose. No model training, unknown experimental population fitting, or stochastic reweighting occurs.

To benchmark *precommitted* population predictions:

    python experiments/score_precommitted_fold_predictions.py --predictions predictions.csv --labels held_out_nmr.csv

This is a scoring interface, not a pre-trained model. A credible prospective test must independently enforce a training-data cutoff, code/data manifest, held-out mutant groups or independent future experiments, and versioned prediction submissions locked **before** any test labels become accessible. GA95/GB95 and even many of the 60 NMR variant labels are already published; fitting to them and scoring on the same variants is not blind validation.

## Scientific gaps still open

- Training or deploying a physical ensemble sampling system and obtaining converged relative free energies for *both folds* (and unfolded states).
- Calibration and falsification of statistical, numerical and model error bars on held-out conformational ensembles.
- Running current 2026 structure predictors / generative ensemble samplers as comparable baselines, not just contact heuristics.
- Finding and ethically validating genuinely unseen mutation cases, with pre-registered wet-lab NMR and/or CD measurements.
- Formal end-to-end bridging of independently computed atomistic energy bounds and conformational populations with the Lean checker.

**Scientific success criterion:** blinded and reproducible improvements in properly scored *probabilistic fold population estimates* with calibrated uncertainty and out-of-distribution mutation generalization. Anything less is a useful step, but not a revolution in protein folding.

## Independently measured contact-map separation and proposed FRET reporters

The authenticated source has now been analyzed using independently computed ring
plane normals as well as ring centroid distances. The 20 GA95 and 20 GB95
NMR models show a clean contact-network separation: within-GA95 contact-map
Jaccard is ≥0.673, within-GB95 ≥0.849, whereas *all 400 between-fold pairs*
are ≤0.106. These are **retrospective** comparisons among already known NMR
structures; no fold prediction from sequence is involved.

GB95's 20-model F30/F52 ring centroids are close and their planes nearly
perpendicular in all 20, while the Y45/F52 pair is close but does not satisfy
the strict parallel-plane or perpendicular-plane criteria used in our geometry
analysis. This directly illustrates why a proximity threshold alone cannot be
called aromatic π-stacking.

We also obtain 349 unmutated residue pairs whose C-alpha distances fall within
12–40 Å in both folds while having nonoverlapping ranges between folds. The
most strongly separated is Q11–V21 (14.49–15.16 Å versus 31.53–33.33 Å).
See `docs/EXPERIMENTAL_DISTANCE_REPORTER.md` for the **conditional** experimental
reporter proposal, hypotheses, chemistry and calibration caveats, and exact
three-state population-inversion calculations.