# GA95/GB95: a mutation-independent structural reporter and robust three-state population inversion

**Status: computational selection of experimental distance-reporter sites from authenticated 2KDL/2KDM NMR models. THIS IS NOT A VERIFIED FRET EXPERIMENT OR A SEQUENCE-TO-FOLD PREDICTOR.**

## Motivation

GA95 and GB95 are 56-residue proteins differing at just three sites (L20A, I30F, L45Y) but adopting distinct folded topologies. An experimentally reliable optical reporter could allow much higher-throughput tests of mutation-dependent conformational equilibria than NMR alone, *if* labeling and photophysical effects are calibrated. A 2026 bioRxiv dataset contains 60 NMR-characterized GA/GB variants (Felbinger et al., DOI 10.64898/2026.08.02.742283). A multi-state optical assay could be valuable for testing structural models on additional variants and environment conditions.

## Authenticated structural inputs and nonoverlap result

The official RCSB PDB files are exactly SHA-256-pinned in the adjacent reference artifact. For each of 20 GA95 and 20 GB95 NMR models, the algorithm computes all nonlocal C-alpha contact edges. Using an 8 Å contact threshold and sequence separation >=4, it obtains the following pairwise model-level Jaccard statistics:

| Comparison across submitted NMR models | Minimum | Median | Maximum |
| --- | ---: | ---: | ---: |
| Within GA95 (190 pairs) | 0.6731 | 0.8242 | 0.9500 |
| Within GB95 (190 pairs) | 0.8492 | 0.9160 | 0.9828 |
| Cross GA95–GB95 (400 pairs) | 0.0695 | 0.0884 | 0.1064 |

These intervals are disjoint, so the experimentally deposited folds possess a robust geometric contact-map separation **within the structural uncertainty represented by those specific NMR models**. This does not establish all-solution thermal ensemble separability or prospective sequence-to-fold classification.

## Optimizing pairs for geometric distance readouts

Scan all unchanged amino-acid positions (exclude 20, 30, and 45), sequence separation ≥4, and demand *every* deposited model places both folds' C-alpha distances in the assay-planning interval 12–40 Å. Retain pairs with disjoint distance ranges, then rank by the WORST-case separation between the two complete model ranges, not mean distance difference.

**349 candidate pairs pass this strict geometric criterion.** Top examples:

| Residue pair | GA95 C-alpha distance range | GB95 C-alpha distance range | Worst-case separation |
| --- | --- | --- | --- |
| **Q11–V21** | 14.49–15.16 Å | 31.53–33.33 Å | **16.37 Å** |
| K10–V21 | 15.83–16.23 Å | 31.90–33.14 Å | 15.67 Å |
| A12–V21 | 13.50–14.12 Å | 29.20–30.38 Å | 15.09 Å |
| **Q11–D22** | 15.40–16.28 Å | 30.31–32.39 Å | **14.03 Å** |

The Q11–V21 C-alpha separation is strong enough that any **hypothetical** combined dye-localization error less than 8.18 Å *per state pair-distance bound* would preserve nonoverlap in this structure-only model. There is no empirical justification that a dye at Q11 or V21 stays that close to the C-alpha atom, nor that substituting dye-attachment residues is fold-neutral. Q11–D22 is a useful secondary candidate because both are polar/charged, not because accessibility or nonperturbation has been established.

For a purely illustrative Förster radius of 25 Å and **unrealistic assumption that dye-dye distance exactly equals C-alpha separation**, the Q11–V21 model would imply emission efficiency proxy bands around 95.3–96.3% (GA95) versus 15.1–19.9% (GB95). These are NOT predicted experimental FRET efficiencies: physical dyes have linkers, orientation factors, kinetics and fluorescence changes.

## Proper three-state populations, not just two-fold fractions

The experimental 2026 benchmark measures alpha-helical, beta-rich and unfolded states. Treating fluorescence as a two-state mixture can be catastrophically wrong whenever unfolded species exist. A new **exact-rational, three-state constrained-population inverse** solves this transparently.

Suppose the measured intensity is in [y_lo,y_hi], with independently *calibrated* emission intervals for three states A, B and U and an externally measured bound p_U≤U_max. Fractions obey p_A+p_B+p_U=1. Linear inequalities ensure the lower possible mixture emission is ≤y_hi and upper possible mixture emission is ≥y_lo. Substitute p_B=1−p_A−p_U, giving a bounded two-dimensional rational polytope. Enumerating intersections of its linear constraints returns the **sharp compatible range** for both p_A and p_U.

As a **completely synthetic demonstration**, if state emission calibrations were A=[0.9,1], B=[0.1,0.2], U=[0,1], U≤0.1, and measured y∈[0.6,0.65], the exact inverted state-A population lies in **[2/5,7/10] = [40%,70%]**. These are not GA95/GB95 observations.

The computation refuses inconsistent measurements and remains conservative when unfolded emissions are broad. If U_max is not justified by measurement, using U_max=1 may produce nearly vacuous bounds; that is correct behavior, not a failure to be patched with an arbitrary prior.

## Verification and experimental acceptance criteria

- 20+20 official-source NMR coordinate models, full raw data hashes, 400 cross-model and 380 same-fold structural comparisons.
- Distance separation ranked across *all* valid contact pairs, not cherry-picked from one NMR model.
- Exact rational three-state linear-programming inversion with brute-force forward Monte Carlo rational containment regressions; no floating-point optimization is trusted to certify population bounds.
- Required wet-lab follow-up: model solvent-accessible sidechain positions; avoid buried or essential residues; assess label attachment mutants by NMR/CD; choose fluorophores calibrated at relevant 15–35 Å distances; calibrate FRET efficiency distributions using fold-locked standards; measure p_U independently; verify against held-out mutants and physical conditions.
- Publication-level criterion: prospective out-of-sample mutant population estimates with uncertainty coverage, plus controls for labeling-induced changes in folding free energies. In the absence of this prospective validation, this is only a structure-guided reporter proposal.

**Novelty:** choosing FRET reporters from structural distance differences and using linear population inference are established ideas. We have not performed a systematic prior-art search for this exact GA95/GB95 reporter pair, and neither novelty nor revolutionary impact is claimed.