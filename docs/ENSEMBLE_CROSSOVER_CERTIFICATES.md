# Exact robust conformational-ensemble crossover certificates

**Status:** independently tested Python research prototype, not Lean formalized, not a realistic protein force field, not a claim of publication-level novelty. The certificate covers a precisely stated two-dimensional finite lattice statistical-mechanical model.

## Scientific objective

A protein's *ensemble* can change with its environment even if its single most likely structure is not informative. Recent conformational-ensemble literature emphasizes that structure prediction alone does not supply equilibrium populations, conformational exchange kinetics or response to ligands/environment.

This study addresses a prerequisite mathematical problem: given a finite conformational model and a **continuum of uncertain energetic parameters**, can we prove that an event is a majority/minority population throughout an entire parameter region without testing finitely many sampled points?

## Toy model

An anchored 2D self-avoiding walk of N sites begins at `(0,0),(1,0)`, takes cardinal steps, and visits no lattice site twice. Each site has type H or P. For conformation x:

- K(x): count of non-backbone, lattice-adjacent hydrophobic/hydrophobic pairs.
- B(x): count of changes between consecutive backbone directions, *including* the first fixed bond.
- Weight: `w(x;q,b) = q**K(x) * b**B(x)`, with rational positive parameters `q`, `b`.
- Event: `A_t = {x : K(x) >= t}`.
- Probability: `P(q,b) = sum(x in A_t) w(x;q,b) / sum(x) w(x;q,b)`.

Here `q` parameterizes a hydrophobic-contact reward and `b` a propensity to form turns. Neither is calibrated to solvent, temperature or a physical protein. A 12-mer contains exactly 30,073 anchored SAWs under this convention.

Define exact bivariate generating polynomials

`E(q,b) = sum_{k >= t,r} n_(k,r) q**k b**r`

`C(q,b) = sum_{k < t,r} n_(k,r) q**k b**r`

where integer `n_(k,r)` is the number of actual conformations with K=k, B=r. Then `P=E/(E+C)`.

## Theorem 1 — monotonicity in hydrophobic contact strength

For **any** finite nonnegative collection of such states and positive q,b, an event defined by the contact-count threshold K>=t has probability nondecreasing in q. More precisely,

`q*(E_q*C - E*C_q) = sum_{x in A, y not in A} (K(x)-K(y)) w(x)w(y) >= 0`.

Every summand is nonnegative because `K(x)>=t>K(y)`. The derivative of `P` with respect to `q` has this sign. This theorem is general for threshold events and exponential-family contact weighting; it is not claimed new.

## Theorem 2 — sufficient exact-arithmetic bending monotonicity witness

For each fixed q, `∂P/∂b = (E_b*C - E*C_b)/(E+C)^2`. The numerator is an integer polynomial in q and b (up to rational common factors). Shift variables to nonnegative coordinates `u = q-q_min` and `v=b-b_min` and expand:

`E_b*C - E*C_b = sum_{i,j} a_(i,j) u**i v**j`.

If **all exact rational coefficients** a_(i,j) are nonnegative, then `P` is nondecreasing in b for every `q>=q_min`, `b>=b_min`. This is a sufficient, not necessary condition, and the checker fails closed if any coefficient is negative. No numerical sampling establishes the sign.

For `HPHHPHHPHHHH`, t=2, `q_min=3/2`, `b_min=1/4`, the computed expansion contains **133 strictly positive nonzero coefficients, zero negative coefficients**. The code recomputes this exact expansion from independently enumerated integer state counts. Thus bending monotonicity holds on this entire unbounded positive quadrant **conditional on the verified finite-model histogram**.

## Theorem 3 — sharp certified box probabilities

If both monotonicities are established and q and b vary independently in `[q_min,q_max] × [b_min,b_max]`, the global exact probability range is:

`P(q_min,b_min) <= P(q,b) <= P(q_max,b_max)`.

Both endpoints are **attained** at corners, so these are the sharpest possible lower/upper bounds for the specified model and rectangle.

For q ∈ [3/2,2], b ∈ [1/4,4], the exact bounds are:

- Lower: `84529143/508499807` = 0.166232399376...
- Upper: `1108499200/1368704059` = 0.809889612521...

For comparison, treating event and complement masses separately using interval extremization, without shared-parameter cancellation, gives:

- Naive lower `3130709/242528510002517`, approximately 1.29×10^-8.
- Naive upper `3487036971417600/3487037024413933`, approximately 0.9999999848.

The new interval is **exact** for the model and input parameter rectangle, not an estimate. The improvement is in the robustness certificate rather than in a physical folding prediction.

## Theorem 4 — unique 50% population crossover at fixed q

For each fixed rational q, define `F_q(b)=E(q,b)-C(q,b)`. Then `P(q,b)>1/2` exactly when `F_q(b)>0`. For the mixed 12-mer at t=2 and q in `{3/2,2,3}`, the nonzero coefficient sequence of `F_q` has exactly one sign change. By **Descartes' rule of signs**, there is at most one positive real root, counted with multiplicity; since the signs at the two rational endpoints of the computed interval differ, there is **exactly one** root. Bisection using exact rational Horner evaluation encloses this unique crossover:

- q=3/2: b* ≈ 1.173422028768
- q=2: b* ≈ 0.647013416866
- q=3: b* ≈ 0.263857423817

Exact rational brackets with widths less than 10^-10 are printed by `experiments/ensemble_crossover_benchmark.py` and preserved as machine-readable results.

Moreover, for **every** q ∈ [3/2,2], monotonicity plus exact endpoint evaluation gives:

- For b ∈ [1/4, 647/1000], **P(q,b)<1/2** throughout the entire rectangle.
- For b ∈ [1174/1000,4], **P(q,b)>1/2** throughout the entire rectangle.

These are two distinct, rigorously classified population-majority regions separated by a band where either status is possible. This is a *finite-system crossover*, **not** a singular thermodynamic phase transition.

## Tests and trust boundary

The project has independent 12-mer anchored-SAW enumeration in a separate pytest implementation, exact rational histogram/polynomial comparisons, root uniqueness/straddling checks, parameter-grid regression, rejection of invalid certificates, and a negative test where the coefficient-based sign witness is unavailable.

- Run tests: `python -m pytest tests/test_ensemble_crossover.py -q` (with Python import path pointing to `src`).
- Run benchmarks: `PYTHONPATH=src python experiments/ensemble_crossover_benchmark.py`.
- Check JSON output in `experiments/ensemble_crossover_results.json`.

**Not proved:** completeness/correctness of the Python path enumerator in Lean, correctness of a 3D or all-atom energy potential, protein thermodynamic population prediction, protein folding kinetics, force-field parameter identifiability, prior-art originality. The coefficient proof is an exact algebraic certificate *given the counted microstates*. No empirical assay or real protein data were used.

## Roadmap toward genuine protein-folding relevance

1. Upgrade from 2D HP to a small, realistic 3D peptide model with bonded/angle/dihedral and solvent terms; benchmark against an experimentally characterized fold-switching mini-protein or a small validated molecular-dynamics ensemble.
2. Replace exhaustive state histograms with compact **proof-carrying generating-function bounds** produced by the existing FoldBounds Lean frontier checker, and verify coefficient inequality certificates in Lean.
3. For energies with many correlated parameters, develop sound adaptive polynomial/interval subdivision rather than assuming monotonicity.
4. Validate probabilities and population switching against NMR, smFRET, or HDX-MS under varying conditions, blind to the experimental readouts during prediction.
5. Only consider it a field-changing breakthrough after demonstrating **prospective improvements in predictive accuracy and calibration** relative to contemporary ensemble-generation and simulation baselines.

## Related literature

- Bowman, *AlphaFold and Protein Folding: Not Dead Yet! The Frontier Is Conformational Ensembles*, Annual Review of Biomedical Data Science (2024), DOI 10.1146/annurev-biodatasci-102423-011435.
- *The landscape of machine learning approaches for modeling protein conformational ensembles*, Current Opinion in Structural Biology 98 (2026), DOI 10.1016/j.sbi.2026.103253.
- *Experiment-guided AlphaFold3 resolves measurement-consistent protein ensembles*, Nature Biotechnology (2026), DOI 10.1038/s41587-026-03166-5.

The basic mathematics of exponential families, partition functions, Descartes' rule, exact enumeration and interval coefficient bounds is well established. Novelty of the *combined workflow* has not been established; a dedicated prior-art search is required.


## Three-dimensional independent replication

The same exact-rational derivative witness also applies to the cubic-lattice
anchored **three-dimensional** HP model. A second path enumerator generates
all **13,565** conformations of eight-residue sequence `HPHHPHHH`, with a
threshold of one nonconsecutive H-H contact. An independent iterative
walker reproduces the full two-feature state histogram. After shifting
`q-3/2` and `b-1/4`, its derivative numerator has **55 strictly positive
coefficients and no negative coefficients**, certifying bending monotonicity
for all `q>=3/2, b>=1/4`. Contact-threshold monotonicity in q follows
from the same general covariance proof. The 50%-population crossover is
uniquely bracketed for q=3/2 and q=2 using Descartes' rule and rational
bisection. The 3D example is a qualitative replication of the method,
**not** atomistic protein thermodynamics and not independently validated
against biological measurements.

Reproduce with `python -m pytest tests/test_ensemble_crossover_3d.py -q`
and `python experiments/ensemble_crossover_benchmark.py`.