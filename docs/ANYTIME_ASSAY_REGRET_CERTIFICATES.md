# FoldBounds: jointly certified adaptive protein assay budgets

**Status:** independent exact-rational Python research prototype and retrospective
S669 evaluation. **Not** a Lean-kernel proof, proof of biological folding
kinetics, model of atomistic protein structure, or evidence of a new physical
law. The group-exchangeability assumption is statistical and is not
machine-checkable. This work adds files to a dedicated research branch and
leaves all other FoldBounds branches and existing methods unchanged.

## Question, novelty boundary and realistic experimental task

A protein-engineering laboratory has a fixed **finite library** of candidate
single-site mutations under a fixed reported pH and temperature, a published
model's predictions, and a limited assay budget. It can measure candidates
sequentially and decides how many to test *after observing the measurements*.
Can it have a valid, prospectively auditable **selection-regret guarantee**
that survives that data-dependent stopping decision?

Individual 90% confidence statements for each assay budget do not in general
retain 90% coverage if the researcher stops by inspecting the experiments.
The statistical device used here is the established **split conformal maximum
score over a prespecified family of decisions**. Similar pathwise,
selection-safe conformal ideas exist in the literature (e.g. Habib & Shihab,
2026, arXiv:2610.20773; existing online selective conformal work by Bao et al.,
JMLR 2025). We claim no originality of the general conformal rank theorem.
The distinct FoldBounds prototype instantiates the technique for physical
mutant-selection regret, reproduces complete protein/motif group verification,
accounts for proposed assays and conservative laboratory errors, and supplies
a second fully independent Python certificate checker.

## Source and physical scope

S669 is a published curated experimental dataset of 669 single-point
mutations; the pinned **GeoStab** prediction table contains original
experimental ddG labels and GeoDDG-3D forward model estimates. Reference:
Pancotti et al., *Briefings in Bioinformatics* 23, bbab555 (2022),
DOI `10.1093/bib/bbab555`. Prediction source: Gonglab-THU/GeoStab Git commit
`be8e6d079ea792396f3a4343a8de54273d3d3828`, file `data/ddG/S669.csv`,
**Git blob SHA-1 `f7a19c2ac0821371a947469231eb03d728ab884f`**.
The loader refuses other bytes. Raw data is intentionally **not** bundled.

A valid experimental comparison panel contains only variants with the same
entire wildtype sequence, reported pH and reported temperature. This is a
necessary but not sufficient condition for a comparable folding free energy;
assay chemistry, reference state, pressure, aggregation, aggregation-mediated
stability, fold switching and experimental uncertainties are not independently
resolved by this table. Higher ddG means more stable under the source's
published sign convention. Within each matched-condition panel, the target
is **best experimentally measured candidate**, not native structure, true
biological function, or the best of all possible sequences.

Of 93 distinct wildtype sequences, 66 have usable panels of at least two
variants: 69 matched-condition panels containing 638 experimental mutants.
A label-blind, exact contiguous 12-residue shared-motif connected-component
blocking rule reduces these 66 eligible proteins to 61 assay-group units
(after building motifs over **all** wildtypes, including ineligible bridge
proteins). This is a leakage stress test, not a homology-family detector:
affinity tags can cause false grouping and distant homologs can be missed.

## Mathematical definitions

For each complete panel P, a *fixed* published model produces scores s_i.
For any predeclared budget fraction b, take k_b = max(1,ceil(b|P|)) highest
scoring candidates with ID ties broken lexicographically. Let S_b be the
selected candidates. Observed experimental stability change is y_i. Define

    R_P,b = max_{i in P} y_i - max_{i in S_b} y_i >= 0,

    M_P,b = max(0, max_{i not in S_b} s_i - max_{i in S_b} y_i).

If all candidates are tested, set M=0 and R=0 exactly. M is a *measurement-
conditioned competitive gap*; it does not use the untested y labels.
Although it combines model scores and experimental assays in the same kcal/mol
units, it is sensitive to unaccounted additive assay baseline shifts. That is
an important limitation, not a new microscopic physical equation.

A complete group g can include multiple pH/temperature panels or, for the
sequence-overlap sensitivity analysis, multiple proteins. Fix in advance:

    B = {1/20, 1/10, 1/5, 1/2};
    lambda = 1;  w_b = 1 for every b in B;  alpha = 1/10.

The entire *one-dimensional, signed* nonconformity score of group g is

    S_g = max_{P in g, b in B} ( R_P,b - lambda M_P,b )/w_b.

One scalar score covers every budget, protein and compatible assay panel
in the group. From m complete calibration groups, take the exact rational
rank k=ceil((m+1)(1-alpha)) and cutoff Q=S_(k). If k>m, Q=INFINITY
(abstention). Signed scores and cutoffs must not be silently truncated to
zero: doing so can reduce efficiency and invalidate verifier replay.

### Simultaneous-selection theorem (standard conformal rank proof)

Assume the *entire future group* and the complete calibration groups
(including their candidate sets, assay conditions, published model scores,
selection rules and underlying measured outcomes) are **exchangeable**;
the score formula and all budget and model parameters were chosen without
looking at calibration/test outcome labels. Then

    Pr[ FOR EVERY P in new group AND EVERY b in B:
         R_P,b <= max(0, lambda M_P,b + w_b Q) ] >= 1-alpha.

Proof: all S_g are exchangeable, so the rank of S_new among m+1 scores is
uniformly distributed under a random tie break, hence at most alpha chance
of exceeding the kth calibration score. On the event S_new<=Q, each
individual (R-lambda M)/w <=Q. Since R>=0, the max-with-zero bound is safe.
The one event covers **all four budgets**: any stopping rule selecting a
budget in B based on previously observed selected assays can choose one of
the already-certified inequalities without further selection correction.
This is not a conditional guarantee for the particular stopping event or
an identified protein; nor does it hold under arbitrary protein-family,
assay, temperature or laboratory distribution shift.

### Noise-robust sensitivity theorem (stated, not estimated, error bars)

Suppose each original *calibration assay* and each future **selected** assay
has unknown adversarial absolute error at most epsilon in kcal/mol. It is
not an empirically measured value in S669 (there are no same-mutation,
same-condition experimental replicates in the pinned table). For each budget,
regret computed from potentially noisy complete source measurements changes
by <=2epsilon, and its competitive gap changes by <=epsilon. Therefore

    eta = max_{b in B} (2+lambda)*epsilon / w_b,
    Q_true <= Q_observed + eta.

For a future panel with measured selected candidates and their observed gap
M_obs, the additional error in M is <=epsilon. Consequently

    R_true(P,b) <= max(0, lambda*(M_obs+epsilon)
                              + w_b*(Q_observed+eta)),

on the same >=1-alpha **marginal group event** for latent true values, if
the underlying *true* complete protein groups are exchangeable and the
stated deterministic assay error caps really hold. With lambda=w=1,
this adds at most **4 epsilon** to the zero-noise bound (before clipping).
For epsilon=0.25 the contingency budget is +1 kcal/mol. It is a symbolic
sensitivity, not proof that any laboratory's epsilon is 0.25.

This logic is deliberately conservative: calibration score errors are
accounted for in the quantile *and* selected measurement errors in the
prospective margin. Arbitrary replicate-specific offsets exceeding the cap
invalidate it. The same test also checks failure on invalid/nonfinite assay
values, unrelated assay conditions, and unseen/counterfeit group identities.

## Actual experimental results (published data, independent numeric replay)

Fixed lambda=1, uniform weights and four budget fractions above. 60
**predefined** whole-group, fivefold partitions (seeds 20261105 through
20261164) with 4/5 of groups used for calibration and the other 1/5
held out on each fold. Model predictions are frozen before the study and
no held-out outcomes enter that fold's quantile. The source data are
*reused* across the 60 partitions; these are not independent replications
or prospective studies. Real-data summary numbers below were independently
recomputed with a JavaScript reference using IEEE double values; the
checked-in Python reproduction uses **exact Fraction/Decimal** arithmetic.
Both are subject to source data labels and model-training overlap caveats.

**Simultaneous 90% four-budget group coverage:** mean 90.66% for exact
wildtypes (66 groups), mean **90.27%** for the stricter shared-motif
components (61 groups). In either case, this is observed held-out coverage,
not evidence of an unconditional guarantee outside group exchangeability.

| Budget fraction | Joint multi-budget upper regret, mean (61 motif blocks) | Four-way Bonferroni upper regret, mean | Individually calibrated 90% only (NOT safe to select retrospectively) |
|---:|---:|---:|---:|
| 5% | **2.186** | 4.465 | 2.119 |
| 10% | **2.124** | 2.890 | 2.056 |
| 20% | **2.007** | 2.699 | 1.755 |
| 50% | 1.807 | **1.774** | 1.023 |

All values in kcal/mol, averaged equally over held-out protein/motif groups
within repeated cross-fits, with each group's panels averaged within group.
The Bonferroni comparator uses alpha=0.025 at each of four budgets, then
applies the union bound. Unlike four individually calibrated alpha=0.1
thresholds, *both* the joint and Bonferroni alternatives support adaptive
stopping across the four options. The joint certificate is markedly narrower
at 5, 10 and 20% budgets, but slightly wider at 50%.
No claim of uniform superiority or a newly discovered statistical principle.

### Retrospective adaptive stopping simulation

With exact wildtypes (66 eligible protein groups), epsilon=0, a **2 kcal/mol
regret target**, and all four budgets predeclared, stop at the earliest
budget whose jointly certified worst-panel regret is <=2; otherwise test
through the 50% budget but record **not certified below 2**.

Across the same 60 reused partitions: an average **57.65/66** groups
reached the target, approximately **112.13 assays** were used per 66-protein
cohort versus **334 assays** when always testing the top 50% within every
matched panel, a **66.4% lower simulated assay count**. Mean observed
coverage of the issued *regret bounds at the chosen stops* was ~62.52/66
(94.7%); this is not a coverage guarantee conditional on stopping. At a
stricter target of **1 kcal/mol, no groups achieved the target** in these
retrospective partitions; no false success claim is made.

This is NOT an actually conducted sequential wet-lab experiment. It assumes
all published predictions, candidates and retrospective measured values were
already available, then hides untested values only from the decision rule.
Laboratory conditions, observation delays, assay selection bias and actual
measurement error have not been assessed independently. The target cutoff
is interpreted as a *marginal* new-group assurance, not a reliable statement
for any particular protein.

## New negative/no-go result: naive common-offset correction cannot help

A physically appealing response to systematic assay offsets is to adjust
all untested model scores by the **mean measured residual** in selected
mutants. For top-k predicted selection,

    U = max unselected model score,
    sbar = mean selected model score,
    ybar = mean selected measurement,
    ymax = max selected measurement.

Then necessarily U<=sbar and ybar<=ymax, and therefore

    max(0, U + (ybar-sbar) - ymax) = 0.

So the supposedly smarter competitive-gap feature vanishes **identically**
after that straightforward common-mode correction. This mathematical no-go
is tested on many synthetic candidate libraries. It is an important
warning against assuming that the extra efficiency from the uncentered
feature is immune to systematic measurement offsets. This follows from
basic inequalities and is not asserted as original mathematical novelty.

## Independent assurance and threat model

- The original prior work's trusted input path checks source Git blob hash,
  single-residue substitutions and assay conditions. There is no hidden
  retraining of the published 3D prediction model.
- The untrusted **certificate producer** emits a predeclared budget grid,
  lambda/weights, sorted exact group scores, group identifiers and exact
  signed split-conformal cutoff.
- An independently written **verifier** rebuilds motif connected components
  from all wildtype strings using pairwise motif overlap (different from the
  producer's union-find), reconstructs condition panels, ranked candidates,
  measured regret, competitive gap and all group scores, and recomputes the
  quantile. It rejects tampering, duplicates, missing groups, malformed
  budgets, inconsistent source measurements and forged ranks.
- A safe deployment wrapper verifies the certificate before producing an
  interval, refuses unmeasured outcome leakage into the proposed candidate
  array, and rejects known exact-wildtype and direct shared-motif calibration
  overlaps. Distant protein family overlap cannot be ruled out by string
  equality, and human laboratory QA is still required.
- A budget not in the **predeclared** set is never labeled as covered by the
  joint certificate. Insufficient calibration gives an infinite/abstain
  result; unsupported Bonferroni comparisons do not silently become finite.
- Exact Fraction numerators/denominators and source decimal strings are
  checked. Normal rendered numbers are rounded; Python tests validate exact
  identities. **No Lean kernel or independently randomized wet-lab pilot**
  certifies the executable or its scientific assumptions.

## Reproduce and compare with prior FoldBounds work

Run 72 adversarial/functional/parameterized checks, including the inherited
24 assay-budget tests (from the parent stacked PR):

    PYTHONPATH=src python -m pytest \
      tests/test_assay_regret_certificates.py tests/test_anytime_assay_regret.py -q

The pinned actual S669 input is NOT distributed in this repository:

    PYTHONPATH=src python experiments/benchmark_anytime_assay_regret.py \
      --download --repeat 60 --seed 20261105 --epsilon 0 \
      --output anytime_s669_60_seed_audit.json

Or pass --csv path/to/the_pinned_original/S669.csv. The experiment uses
source hash verification and can audit all certificates independently.
For purely hypothetical worst-case 0.25 kcal/mol error, run --epsilon 1/4.

## Honest research status and next validation

Our contribution is an auditable **protein-engineering decision-assurance
workflow** and retrospective validation of its potential cost efficiency,
not new conformal mathematical priority, accurate physical protein folding,
a clinically or commercially ready assurance product, or an actually
performed biological experiment. Protein family grouping remains imperfect;
model training overlap, experimental repeatability, mutations from different
biochemical states, and truly independent 2026 protein-stability datasets
are not verified here. Recent independent benchmarks have warned that no
single UQ or acquisition method consistently wins on protein engineering
landscapes; see *Benchmarking uncertainty quantification for protein
engineering* (2025, PMC11741572).

Before scientific priority or publication claims, repeat this full pipeline
on independently obtained, well-labeled, **prospective** mutation libraries
with replicate assays, reliable laboratory error bars, protein-family
exclusion, multiple sequence/structure predictors, and prespecified budget
and stopping targets. A future Lean proof of the generic conformal rank
argument would verify mathematics, not assumptions about biology.
