# Experimental protein-stability assurance: physical reciprocity + whole-protein coverage

## Scope and provenance

**This is a serious real-data evaluation, not a claim to have solved protein folding.**
It adds only an isolated research module, test suite, benchmark, and evidence.
Neither the existing formal thermodynamics, kinetic work, motion benchmarks,
nor any other parallel research branch is changed.

Data: Pancotti et al. (2022), Briefings in Bioinformatics 23:bbab555,
DOI https://doi.org/10.1093/bib/bbab555, measured S669 stability changes.
Prediction table: Gonglab-THU/GeoStab, Git commit
be8e6d079ea792396f3a4343a8de54273d3d3828, file data/ddG/S669.csv,
Git blob SHA-1 f7a19c2ac0821371a947469231eb03d728ab884f.
The loader refuses any other source bytes. No third-party raw measurements
are rehosted in this PR.

669 source variants are single substitutions, spanning 93 unique wild-type
sequences and 94 PDB chain IDs. Temperatures range 8–74 as recorded in the
source; pH spans 5–8. Source ddG sign conventions remain unchanged.
The source contains 22 published pairs of forward/reverse model predictions.
The group coverage protocol preselects five methods with no missing fields:
GeoDDG-3D, GeoDDG-Seq, ACDC-NN, DDMut, Dynamut.

## Theorem A: exact physical reciprocity projection

With thermodynamic free energy of folding defined at the SAME conditions,
the true forward and reverse stability differences are y and -y.
For submitted forward and reverse prediction (a,b), the Euclidean projection
onto physically reciprocal prediction pairs is (t,-t), t=(a-b)/2.

For ALL real a,b,y the exact identity is

    (a-y)^2+(b+y)^2 = 2*(t-y)^2+(a+b)^2/2.

Thus enforcing reciprocity cannot worsen COMBINED squared forward/reverse
prediction error. It CAN worsen forward-only squared error or absolute error.
Pancotti et al. already studied antisymmetry; the law itself is prior art.

Observed direct vs projected-forward MAE, published S669 table (kcal/mol),
using the same nonmissing rows for each method:

| Model | Rows | Direct | Reciprocal projection |
|---|---:|---:|---:|
| GeoDDG-3D | 669 | 0.9600 | 0.9600 |
| GeoDDG-Seq | 669 | 0.9683 | 0.9683 |
| DDMut | 669 | 1.0612 | 1.0525 |
| PremPS | 659 | 1.0806 | 1.0409 |
| FoldX | 659 | 1.5582 | 1.3560 |
| MAESTRO | 669 | 1.0619 | 1.2080 |

This is not a new benchmark leaderboard or a state-of-the-art accuracy claim.
The forward MAE can increase, as MAESTRO illustrates.

## Theorem B: finite-sample whole-protein split conformal

Group every observed variant by the entire wild-type amino-acid sequence,
NOT by a randomly assigned mutation or only a PDB structure identifier.
For a frozen prediction function f chosen using training proteins alone,
compute for each calibration protein g:

    R_g = max over its measured mutations j of |y_gj - f(x_gj)|.

With m distinct calibration protein groups and alpha in (0,1), let
k=ceil((m+1)*(1-alpha)) and Q be the kth ordered group residual. If k>m,
Q is INFINITY (abstain). Otherwise each new protein's mutations receive
the interval [f(x)-Q, f(x)+Q].

Assuming whole protein groups are exchangeable with calibration proteins
and the predictor was frozen independently of calibration and test labels,
the exchangeable-rank proof gives:

    Pr(all measured mutants of a random NEW protein fall inside
       their intervals) >= 1-alpha.

This guarantee is **marginal across new protein groups**, not conditional
for any identified protein, not a deterministic physical interval, and not
automatically valid under pH/temperature, assay or evolutionary shift.
It does not require independence among mutants within a protein.

A usual mutation-level conformal quantile is included as a NEGATIVE CONTROL:
it does NOT provide the simultaneous whole-protein guarantee.

## Locked reproducible evaluation

Five folds of whole protein sequences are assigned after seeded Fisher-Yates
xorshift32 shuffle (seed 20261105). For fold i, TEST uses i, CALIBRATION uses
i+1 modulo 5, while TRAIN consists of the other three folds. The same
wildtype never appears in two partitions, including if it has two PDB IDs.

Predictors:
1. Frozen published GeoDDG-3D direct (no fitted choices).
2. Train-only selected reciprocal projection among five complete models;
   selection optimizes macro-averaged protein MAE on TRAIN only.
3. Train-only selection among direct and projected models (sensitivity).
Both nominal 90% and 80% group coverage levels are evaluated.
All inputs are Decimal and quantiles use exact integer ranks.

For seed 20261105, **nominal 90%**:

| Calibration and predictor | Whole proteins fully covered | Variants covered | Mean interval width |
|---|---:|---:|---:|
| Mutation-only negative control, GeoDDG-3D | 62/93 (66.67%) | 593/669 (88.64%) | 4.762 kcal/mol |
| **Protein group**, frozen GeoDDG-3D | **86/93 (92.47%)** | 655/669 (97.91%) | **8.291 kcal/mol** |
| **Protein group**, train-selected reciprocal | **86/93 (92.47%)** | 654/669 (97.76%) | **8.471 kcal/mol** |

With the frozen model, only 18/669 intervals exclude zero (2.7% sign
decision yield); all 18 have the experimentally observed sign. This is
retrospective coverage, NOT proof that future 18 calls are infallible.

At nominal 80%: frozen group coverage 75/93 (80.65%) and mean width
6.425 kcal/mol; selected reciprocal group coverage 70/93 (75.27%) and
mean width 6.360 kcal/mol. A single held-out realization falling below
the nominal fraction does not refute an exchangeable-rank probability bound.

Across 60 FIXED predefined shuffle seeds 20261105..20261164
(reuse of the same data, **not** 60 independent prospective experiments):
- Frozen group conformal: mean 85.017/93 groups covered (min 77, max 90);
  mean width **8.481 kcal/mol**.
- Train-selected reciprocal: mean 85.033/93 (min 77, max 90);
  mean width **8.518 kcal/mol**.

There is **no robust predictive benefit** from train-set model selection.
The wider intervals are a very real problem: at T=298.15 K, a total
8.29 kcal/mol free-energy interval spans an approximately 1.2-million-fold
range of Boltzmann equilibrium constant ratios. The exponent sign depends on
whether folding vs unfolding free energies are used, but the magnitude of
the uncertainty ratio does not.

## Reproduction

Install existing FoldBounds pytest dependency and run:

    python -m pytest tests/test_experimental_reciprocal_conformal.py -q

Then run with the exact original external data:

    PYTHONPATH=src python experiments/benchmark_experimental_reciprocal_conformal.py --download --seed 20261105 --repeat 60 --output s669_audit.json

Offline: replace --download with --csv /path/to/original/S669.csv.
Incorrect source Git blob bytes cause a fail-closed error; no silent
data alteration and no data reshuffling to optimize the test set.

## Major limitations and next tests

These results are not new experimental measurements, prospective
prediction of unseen proteins, or an atomistic physics model. Even though
the source table is experimental, the model predictors are published
black boxes with training sets not reaudited here; this is an assurance
layer using existing predictions. Grouping by literal sequence eliminates
exact duplicates but is NOT a complete homolog-cluster exclusion test.
A true prospective, homology-aware benchmark with additional datasets,
assay consistency, uncertainty-aware experimental replicates, and real
structural energetics is required before proposing clinical or real protein
engineering use.

Mathematical split-conformal coverage relies on group exchangeability;
real biology can violate it severely. When calibration groups are too few
for a target guarantee, the algorithm correctly returns an unbounded
interval. This is Python code with unit and adversarial tests, not
Lean-kernel verified.

**Conclusion:** The result finds a large, measurable gap between ordinary
mutation-level confidence and scientifically safer whole-protein coverage,
while explicitly exposing the cost in practically unusable interval
width. It is a serious auditable research finding—not yet a major discovery
in protein folding thermodynamics or prediction accuracy.
