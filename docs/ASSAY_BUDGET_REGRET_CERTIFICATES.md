# FoldBounds: experimentally calibrated assay-budget regret certificates

**Status:** real S669 mutation-stability benchmark, independently testable Python
reference module. **Not** an atomistic folding breakthrough, new energy law,
prospective clinical certification, or Lean-kernel-verified implementation.

## Why this is a distinct parallel branch

The existing FoldBounds workstreams cover rigorous lattice partition functions,
protein-motion graphs, exact 3D toy models, and S669 calibration of *absolute*
mutation free energies. This branch addresses a different and more directly
useful research question: **If I can assay only k out of n candidate mutations,
how far might the best measured result be from the best available candidate?**
It changes no existing files and does not depend on any other agent's draft PR.

## Real experimental data and science boundary

We use the manually curated experimentally measured single-mutant stability
change benchmark **S669** (Pancotti et al., *Briefings in Bioinformatics* 23,
bbab555, 2022, DOI `10.1093/bib/bbab555`) and the published GeoDDG-3D direct
prediction values (2024, DOI `10.1038/s43588-024-00716-2`) in:

- Repository: `Gonglab-THU/GeoStab`.
- Exact Git commit: `be8e6d079ea792396f3a4343a8de54273d3d3828`.
- File: `data/ddG/S669.csv`.
- Git blob SHA-1: `f7a19c2ac0821371a947469231eb03d728ab884f`.

The reference loader rejects any modified source bytes and checks every row is
a single amino-acid substitution. It is restricted to the original reported
experimental delta-delta-G convention: **negative means destabilizing**;
larger means more stable. The source contains **669 measurements**, **93 exact
wildtype sequences** and records of temperature and pH. To avoid comparing
physically unlike free-energy changes, form a **condition-matched panel** only
among substitutions with identical wildtype, reported pH and temperature
strings/numbers. Do not pool pH 6 and pH 7, or 25 and 56 in a single optimum.

There are **66 eligible proteins** with at least one panel of >=2 measured
variants, totaling **69 matched-condition panels** and **638 candidate
measurements**. Missing full assay chemistry and unreported conditions could
still render even matching pH/temperature panels incomparable. This is a
limited retrospective study; original predictors' training overlap was not
independently reaudited here, and no physical or biophysical equations were
fit to new measurements.

## Core mathematical object: decision regret

For a protein with a fixed experimental panel `P`, let `y_i` be the measured
experimental stability change in kcal/mol and `s_i` be the fixed published
prediction. For assay budget `b` (a rational fraction of candidate count),
choose the k = max(1, ceil(b*|P|)) top-ranked variants by s_i (ties by
variant ID, with no access to measured truth). Define

    R(P, b) = max_{i in P} y_i - max_{i in TOP_k(s, P)} y_i >= 0.

For each protein g with multiple matched-condition panels, take **maximum**
regret across panels `R_g=max_{P in g} R(P,b)`. This is a different target
from whether the model predicts accurate individual absolute free energies.
Systematic per-protein offsets in s cancel from ranking.

### Proof-carrying confidence certificate

Given m calibration protein groups, a ranking rule frozen independent of their
observed mutation outcomes, and chosen `alpha in (0,1)`, recompute each group
regret exactly (the proposed selection is independently replayable). Let
`q=ceil((m+1)*(1-alpha))` and let Q be the qth sorted calibration regret. If
q>m, **abstain** with an infinite cutoff; never invent a finite guarantee.

Assuming **entire new protein groups** (including their candidate sets, assay
conditions, available mutation counts, predictions, outcomes, and selection
rule) are exchangeable with calibration groups, the standard conformal rank
argument yields

    P(R_new <= Q) >= 1-alpha.

This is a **finite-sample marginal guarantee for a randomly exchangeable new
whole protein**, not a deterministic physical inequality for an identified
protein or a conditional probability guarantee for chosen scores. It measures
regret only with respect to the variants in the candidate library, not among
all theoretically possible sequences. A large prospective domain shift, assay
chemistry shift, or a new library design can invalidate group exchangeability.

Because the published predictor is entirely frozen, **4/5 of the protein
groups are used for calibration for each held-out fold**, leaving 1/5 for
strictly held-out assessment; no model-training labels are needed. A given
protein's outcome is never used for *its own* fold's cutoff, although its
outcome can calibrate a different fold. This is cross-fitted retrospective
evaluation, not independent external validation.

### Exact, fair uniform-random baseline

For a panel with sorted experimental values `y_1 <= ... <= y_n`, if k
variants are selected **uniformly without replacement**, the rank of the
maximum selected value is j with exact probability

    C(j-1, k-1) / C(n, k), j = k ... n.

Thus the **expected** regret of random selection is exactly

    max(y) - sum_{j=k}^n y_j C(j-1,k-1)/C(n,k).

The Python code evaluates this as **exact Fraction arithmetic**, converting
original decimal measurements to exact rationals. There are no lucky random
comparison seeds. Aggregate the model regret and random expectation as
**within-panel means**, then average equally across eligible proteins.
Compute an exploratory paired protein-cluster bootstrap CI; it is separate
from the conformal rank guarantee and is not a claim about new data.

## Independently retrieved real S669 findings

The data below were independently reanalyzed using the pinned exact source
CSV and 60 predefined, reused fivefold partitions (seeds 20261105..20261164).
`b` is the *per-panel fractional assay budget*, so ceiling rounds the actual
fraction of assayed candidates upward.

| Fraction b | Mean assays / protein | Model mean regret | Random expected regret | Paired mean gain | Mean 90%-nominal protein-level regret Q |
|---:|---:|---:|---:|---:|---:|
| 5% | 1.24 | 0.797 | 1.282 | +0.485 | 2.572 |
| 10% | 1.59 | 0.717 | 1.189 | +0.472 | 2.249 |
| **20%** | **2.41** | **0.536** | **1.039** | **+0.503** | **1.630** |
| 50% | 5.06 | 0.192 | 0.529 | +0.336 | 0.913 |

All energies are **kcal/mol**. The 20% budget actually selects 159 of 638
available eligible variants (24.92%) because of ceiling, over 69 panels.
For 20%, an exploratory 10,000-resample paired protein-cluster bootstrap
of the model's advantage over **uniform random K-subsets** gave a 95%
percentile interval of approximately **0.219–0.827 kcal/mol**. The 60 splits
are **not** 60 independent datasets and should not be treated as such.

With 20% candidate budget, nominal 90% **protein group** conformal coverage
across the 60 repeated fivefold studies averaged **90.23%** (held-out
seed-20261105: **60/66** groups covered). The mean cutoff was **1.630
kcal/mol** (seed-20261105 1.660). These figures assume new protein-group
exchangeability, as stated above. Covering 90% does not establish a
conditional bound for the specific chosen protein.

### Sequence-overlap stress test

A stricter, predeclared label-free stress test links wildtype sequences if
they share an identical **12-amino-acid contiguous motif**, then forms
transitive connected components. These are NOT true evolutionary homology
families; common affinity tags can cause false links, and remote homologs
can escape the rule. From the 66 eligible proteins there are **61 eligible
motif components**, formed using all 93 wildtype sequences before eligibility
filtering. Simultaneously validate max regret across all eligible member
proteins in a component, not a single mutation.

At 20% budget, across 60 split seeds: **90.16%** of complete held-out motif
components had regret <= calibrated cutoff, with mean **1.696 kcal/mol**;
for seed 20261105 **55/61** components were covered (mean cutoff 1.759).
The result is weaker than a real pre-registered homology-aware experiment.

### Important distinction from S669 absolute-energy assurance

Prior FoldBounds absolute-error studies found average widths of ~8 kcal/mol
for simultaneous mutation free-energy predictions. A 1.63 kcal/mol **regret
upper bound is NOT an equally valid 1.63 kcal/mol absolute-energy interval**:
it bounds a different operational estimand (relative quality of selected
candidate experiments). Do not advertise this as a new free-energy predictor
or a fivefold reduction in molecular prediction error.

## Test instructions

Existing Python 3.11+ stdlib and pytest:

    PYTHONPATH=src python -m pytest tests/test_assay_regret_certificates.py -q

To rerun pinned full real-data experiment, use the exact external source:

    PYTHONPATH=src python experiments/benchmark_assay_regret.py --download --seed 20261105 --repeat 60 --output real_assay_regret.json

Or supply `--csv /path/to/S669.csv` (SHA-1 still enforced).
The GitHub version must reproduce independently checked source-blob results;
CI on the draft branch and an independent second dataset are necessary
before any strong predictive or prospective claim.

The adversarial tests specifically cover: missing/wrong data provenance,
invalid mutations and duplicated IDs, condition mismatch, omitted tiny
panels, label-blind rank selection, all budget edge cases, exact
combinatorial K-subset enumeration for short panels, quantile index and
infinite-abstention correctness, rank-exchangeability argument, motif
transitive closure and tag collisions, no train/cal/test leakage,
permutation-invariant results, and source-outcome tampering attacks.

## What is NOT established

- No physical folding rates, transition-state energetics, realistic MD, or
  atomically resolved structures are predicted here.
- S669 is a much-studied retrospective benchmark, not a prospectively sealed
  independent second dataset. Original model training overlap is not audited.
- Matching pH and temperature is not enough to prove fully identical assay
  conditions (salt, denaturant, instrumentation, buffer, protocol).
- The conformal theorem is a **classic method applied to a new decision
  quantity**. We do not claim originality of conformal prediction, rank
  theory, selection-by-prediction or biochemical free-energy estimation.
- Exact-rational combinatorics and Decimal measurements avoid floating
  arithmetic errors, but do **not** turn noisy laboratory measurements into
  physically exact values or a kernel-verified program. Lean formalization
  and independently reviewed verification remain future goals.

### Research next step

A genuinely consequential breakthrough would need validation on separately
collected, matched-assay mutation panels, stronger sequence-homology
clustering, measurement-repeat uncertainty, cost-aware experimental design,
and protein-specific domain shifts. Tsuboyama et al. 2023 Nature mega-scale
stability scanning offers an independent high-throughput dataset, but the
underlying cDNA-display proteolysis assays and two-state assumptions are
not directly interchangeable with literature-curated S669 free energies.
See DOI `10.1038/s41586-023-06328-6`.
