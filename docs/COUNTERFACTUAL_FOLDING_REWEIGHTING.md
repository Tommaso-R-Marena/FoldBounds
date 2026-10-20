# Counterfactual ensemble reweighting and the energy-accuracy gap

**Research status (2026-10-20):** Exact algebra and exhaustive reference checks in an 8-/10-residue **cubic H/P lattice model**. A small, independent Lean formalization covers the two-microstate covariance identity if and only if GitHub CI passes. No inference here is an experimentally predictive real-protein folding model; none of the statistical-mechanical identities is claimed to be new.

## The target application

In protein engineering, we want to evaluate how mutations change **relative** conformational populations while avoiding a new complete ensemble simulation for every possible mutation. AlphaFold-like structure predictors do not alone provide calibrated equilibrium populations. State reweighting, free-energy perturbation, and thermodynamic cycles are established approaches; the proposed contribution in this research line is making them certificate-friendly, including exact mutation-local sufficient statistics, fail-closed overlap assumptions, and independently tested uncertainty accounting.

## Exact counterfactual identities

Let a set of physically possible conformations be shared between parent and mutant. For baseline energy weights `w(s)>0`, and a mutation-specific likelihood ratio `R(s)=w_mut(s)/w_parent(s)>0`, exact identities are

- `Z_mut/Z_parent = E_parent[R]`.
- `odds_mut(A)/odds_parent(A) = E_parent[R | A] / E_parent[R | complement(A)]` for any mutation-independent event A and nonzero basin populations.

This holds for *any* positive reference weights, including an arbitrary common model error factor that affects both sequences. It does **not** mean the numerical ratio is unaffected by that error: the corrected reference distribution itself changes. Nor does it imply correctness if the mutant introduces new conformations outside the sampled reference support.

For two substitutions at positions i and j that **never contact directly**, suppose the microscopic mutation ratios factor `R_ij(s)=R_i(s) R_j(s)`. On a bipartite cubic lattice, nonconsecutive same-parity residues cannot be neighbors, so this factorization holds for pairs at positions of identical parity in the HP contact Hamiltonian.

Let `mu_g(·)` and `Cov_g(·,·)` denote expectation and covariance conditional on basin g (event or complement) in the reference ensemble. The thermodynamic interaction in each basin is

`I_g = E_g[R_i R_j] / (E_g[R_i] E_g[R_j]) = 1 + Cov_g(R_i,R_j)/(E_g[R_i] E_g[R_j])`.

The **odds epistasis factor** is exactly `I_A/I_complement`, not necessarily 1. Thus two mutations that never have a direct energetic interaction can exhibit strong nonlocal epistasis by redistributing the ensemble. In the fully general (direct-contact) case, the checker recomputes `R_ij` exactly; it does not silently assume the factorization.

## Controlled exact-lattice experiments

At `q=8`, `b=1`, mutation sites 1 and 3, `PHPHPHPH` and `PHPHPHPHPH`:

- 8-mer event probabilities: parent ~0.203465, site 1 mutation ~0.720332, site 3 mutation ~0.586778, both mutations ~0.863968. Direct contact impossible, yet odds epistasis factor ~0.44357.
- 10-mer event probabilities: parent ~0.132474, site 1 ~0.632882, site 3 ~0.518661, both ~0.815300. Direct contact impossible; odds epistasis factor ~0.36287.
- Exact full search over all 768 same-parity P-to-H double-mutation cycles on 8-mers: **504 nonzero** (140 positive and 364 negative), 264 zero. The common 13,565-state geometry contains only 171 lossless contact/bend fingerprints; at 10 residues, 308,981 states compress to 1,852 fingerprints. The optimizer uses exact rational weights; benchmark recomputes all 768 cycles.

The exact ratio identity is tested against **four independently computed sequence-specific partition functions**, not only algebraic self-consistency. Independent lattice enumerators separately validate the original 8-mer and 10-mer histograms in the parent projects. Direct force-field predictions of real proteins have not been made.

## A crucial physical calibration and falsification result

At 20 °C, `RT ≈ 0.58255 kcal/mol`. The previous 20%-weight perturbation `w' / w ∈ [5/6,6/5]` is equivalent to allowing **only ±RT ln(1.2) ≈ ±0.10621 kcal/mol** of arbitrary error in each microstate's energy (after a common additive energy shift). This is **not** evidence that the toy-model force field is physically accurate to that precision.

The previous continuum parameter-box mutation results can tolerate at most approximately **±0.22916 kcal/mol** (8-mer) and **±0.11188 kcal/mol** (10-mer) of independent statewise errors before their worst-case population-majority guarantees become inconclusive. These are *critical thresholds for the sufficiency of our enclosure*, not observed physical errors and not the precise robustness thresholds of the true physical folding systems.

Published GA95/GB95 experiments show that a single amino-acid substitution can reverse the relative populations of two radically different folds, leaving a >90% dominant state on each side at pH 7.2 and 20 °C (He et al., PNAS 2009, doi:10.1073/pnas.0906408106). Two-state equilibrium thermodynamics then implies a relative-fold free-energy shift **greater than 2RT ln9 ≈2.56 kcal/mol** for that particular experimental switch. This is an *inferred lower bound from experimental populations*, **not a model prediction** of GA95/GB95 by FoldBounds.

The 2025 Chemical Science experimental study identified residues 20, 30, and 45 and aromatic long-range interactions as relevant (doi:10.1039/D4SC04951A). Existing studies already compute relative fold free energies and mutation epistasis with mature physical techniques. We must benchmark against those methods before suggesting originality.

## What would be required for a breakthrough

1. Obtain provenance-checked experimental conformational ensembles (e.g. 2KDL/2KDM NMR), and full sequences and assay conditions; coordinates could not be fetched from this container, so none were reconstructed or guessed.
2. Compute meaningful **energy differences between two competing experimental fold basins** using a calibrated atomistic model, extensive sampling, and equilibrium tests, not static contact counts.
3. Test blinded mutations and temperatures, including negative controls and known fold-switching mutants; assess calibration, RMSE, coverage and confidence intervals against experimental NMR or thermodynamic populations.
4. Detect missing conformational support and sample correlation; in a realistic trajectory, the exact reweighting identity is not a guarantee that finite correlated samples approximate the required expectations.
5. Formalize the sum-to-conditional-expectation correspondence and executable certificate refinement inside Lean. The existing Lean covariance identity alone does **not** verify Python's geometry enumeration.

## Reproduction

```bash
python -m pytest -q tests/test_counterfactual_reweighting.py
python experiments/counterfactual_reweighting_benchmark.py > /tmp/reweighted.json
python - <<'PY'
import json
from pathlib import Path
assert json.loads(Path('/tmp/reweighted.json').read_text()) == json.loads(Path('experiments/counterfactual_reweighting_results.json').read_text())
PY
lake env lean FoldBounds/CounterfactualEpistasis.lean
lake env lean FoldBounds/CounterfactualEpistasisAudit.lean
```

The cited experiment code and published structures are external dependencies; this overlay is intended to sit on the unmerged FoldBounds mutation research branch PR #35. It is **not a drop-in replacement** of `main`.

## Real GA95/GB95 representation identifiability test

Authentic PDB FASTA sequences (2KDL GA95, 2KDM GB95; RCSB `https://www.rcsb.org/fasta/entry/2KDL/display` and corresponding 2KDM entry) are 56 residues and differ at **L20A, I30F, L45Y**. Under a coarse two-class encoding treating leucine, alanine, isoleucine, phenylalanine and tyrosine as a hydrophobic class, **both sequences encode to the exact same H/P string**. Consequently *any* algorithm whose sequence input contains only that encoding and whose remaining terms are sequence-independent produces identical energies and folding probabilities for both sequences, regardless of accuracy or proof complexity.

Tyrosine's polar hydroxyl makes its binary assignment ambiguous: marking Y as polar yields one visible H/P change at position 45, still losing L20A and I30F. Introducing a separate aromatic class reveals changes at 30 and 45 but still hides side-chain-size change L20A. Thus a chemically richer alphabet (and appropriate physical interactions) is required to even *represent* the GA/GB experimental mechanism. This is a simple representation no-go/negative control, not a prediction of the structures. See `src/foldbounds/ga_gb_indistinguishability.py`, fixtures and tests.

Real structure coordinates for 2KDL/2KDM could not be downloaded in the current tool runtime, so the research does not contain any invented NMR coordinates, contact maps, or calculated atomistic free energies.