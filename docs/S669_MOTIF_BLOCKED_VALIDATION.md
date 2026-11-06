# Stress-testing protein stability confidence with sequence-motif blocks

**This is a stricter, experimentally grounded validation stress test, not
a homology search or a claimed new folding-prediction algorithm.**

This additive workstream uses the exact S669 experimental file already
pinned in the parent experimental assurance module: GeoStab Git commit
be8e6d079ea792396f3a4343a8de54273d3d3828, path data/ddG/S669.csv,
blob f7a19c2ac0821371a947469231eb03d728ab884f.
No new source observations, prediction scores or assay data are invented.

## Why literal exact wildtype sequence grouping is insufficient

The 669 experimentally measured single substitutions represent **93**
distinct exact wildtype sequences; however different wildtype sequences
may be closely related, be different constructs of the same protein, or
include a shared affinity tag. Grouping only exact identical sequences can
therefore leave obvious sequence overlaps between training, calibration
and test proteins.

This stress test forms an undirected graph on all 93 unique wildtypes:
two distinct sequences are linked if they share ANY exact 12-residue
contiguous motif. Then it uses the **complete connected components** of
that graph as indivisible units, including transitive links.

This test produces **81 motif components** from the 93 wildtypes.
Some motif links are probably affinity tags, not evolutionary homology,
which is why we do NOT label these as biologically inferred families or
claim absence of distant homology leakage. The rule is simple, fixed
before examining outcomes, conservative, reproducible and label-free.

## Statistical claim and limits

Inside each calibration motif component, compute the MAXIMUM absolute
prediction error across EVERY measured mutation in ALL member sequences.
Use standard split-conformal quantile k=ceil((m+1)*(1-alpha)).
Assuming **the entire new motif components are exchangeable**, this
produces marginal simultaneous coverage for the full mutation panel of a
random future component. It is NOT a guarantee for a particular protein,
for independently chosen mutations within a component, or when sequences,
assays and experiment conditions have changed distribution.

The 5-fold split is identical in logic to the parent study but applied
to motif components: seeded xorshift32 Fisher–Yates shuffle, fold i
as test, (i+1)%5 as calibration, remaining three folds training/selection
only. For source traceability, the already published 3D model is kept
fixed, with one train-only model-selection sensitivity experiment.
A mutation-level calibration is used only as a negative control.

## Actual experimental results: S669, q is NOT involved

Using seed 20261105, nominal 90% **whole motif component** coverage:

| Approach | Entire motif components covered | Mutations covered | Mean full interval width |
|---|---:|---:|---:|
| Naive mutation-level calibration | 54/81 (66.67%) | 599/669 | 4.334 kcal/mol |
| **Correct motif-component calibration**, fixed GeoDDG-3D | **74/81 (91.36%)** | 652/669 | **8.783 kcal/mol** |
| **Correct motif-component calibration**, train-selected reciprocal | **76/81 (93.83%)** | 662/669 | **10.701 kcal/mol** |

The selected reciprocal strategy is wider and calls the sign in only
one case, despite better observed coverage on this split. It is not a
better practical mutant predictor. The fixed strategy calls sign in
23/669 cases, all 23 retrospectively correct. Again, **these are not
individual correctness guarantees**.

Across 60 predetermined seeds 20261105..20261164, each reusing the same
source observations, not 60 independent test datasets:

- Fixed model: mean **76.317/81 components** covered, min 70, max 80,
  average interval width **9.728 kcal/mol**.
- Selected reciprocal: mean **76.383/81 components**, min 71, max 80,
  average width **9.978 kcal/mol**.

Compared with the literal-sequence split, this stress test exposes an
additional, measurable increase in the uncertainty price of being more
conservative about sequence overlap.

## Methods and falsification checks

The independent Python implementation constructs motif-connected components
using a disjoint-set union algorithm, validates transitive closure,
allocates entire blocks to folds, refuses over-small calibrations and
preserves exact Decimal residual arithmetic and integer conformal order
statistics. It is packaged in separate NEW source/tests/docs/experiment
files and does not change the original experimental assurances.

Tests include transitive chains A~B~C, exact sequence aliases, multiple
predeclared motif lengths, split disjointness, row-vs-block quantiles,
independence of input ordering, invalid k/mode rejection, and shared
affinity-tag examples.

Reproduction with the pinned source CSV:

    python -m pytest tests/test_experimental_motif_blocks.py -q
    PYTHONPATH=src python experiments/benchmark_s669_motif_blocked.py --download --k 12 --seed 20261105 --repeat 60 --output s669_motif_audit.json

## Next scientific tests

A real homology-sensitive evaluation should replace the exact motif
surrogate with sequence alignment or structure/domain clustering, classify
and remove shared purification tags, exclude high-identity structural
homologs, and then evaluate on genuinely prospective new data measured
under matched thermodynamic conditions. A 12-residue shared motif cannot
establish evolutionary relatedness or absence of distant homologs.
This work tests assurance limitations, **not a new biophysical mechanism,
a breakthrough mutation predictor, or calibrated deterministic physics**.
