# Coupled Ensemble Certificates: mathematical discovery report (Deliverable B)

**Status (2026-11-13):** complete for this session; independently reviewed in two
phases (ADVERSARIAL_REVIEW.md). Nothing here is claimed as novel
beyond what `PRIOR_ART.md` supports.

## 1. Question and answer in one paragraph

**Question.** Can rigorous, checkable bounds on a relative free energy
`ΔF = −β⁻¹ log(Z_B/Z_A)` be obtained without enumerating the conformational
space, when the individual partition functions cannot be bounded usefully?

**Answer (proved, machine-checked).** Yes, exactly when A and B share a
decomposition in which the *per-region ratio* is well determined.
- For a shared region decomposition, the sharp coupled enclosure has a
  single-inequality certificate. It is at least as tight as bounding Z_A and
  Z_B separately, and strictly tighter except in a precisely characterised
  degenerate case.
- For perturbations confined behind a separator ("Markov blanket"), local data
  confine the ratio to `[min_τ ρ(τ), max_τ ρ(τ)]`, which is sharp. No sound
  local procedure can do better. The same local data cannot bound Z_A at all:
  the realisable pairs (Z_A, Z_B/Z_A) are exactly (0, ∞) × [min ρ, max ρ].

**What this does not give.** An efficient algorithm in general. Whether local
certificates converge is exactly weak spatial mixing (T2.8), and below a phase
transition they provably cannot converge. For lattice polymers there is no
separator, and the method gives essentially nothing beyond the trivial
‖V‖∞ bound (measured; §5).

## 2. Results and their assurance status

| ID | Statement (informal) | Status | Novelty risk (PRIOR_ART.md) |
|---|---|---|---|
| T1.1/1.2 | λ ≤ Q on all admissible configurations ⟺ g₋(λ) ≥ 0 (exact test) | Lean: `forall_le_Q_iff`, `forall_Q_le_iff` | high (linear-fractional duality) |
| T1.3 | sharp enclosure attained at box vertices; g₋(λ₋) = 0 | Lean (attainment form); the sorted-threshold form is not formalized | high |
| T1.4/1.5 | coupled ⊇ separate dominance; equality ⟺ r_k(U_k−L_k) = 0 ∀k (needs S₋ > 0) | Lean: `coupled_dominates_lower`, `lamLo_eq_sepLo_iff` | medium as stated |
| T1.6/1.7 | uniform-ratio exactness; monotone refinement (any splitting map) | Lean | high |
| T1′ | LP-dual certificate soundness | Lean: `le_Q_of_lowerDualCert` | high |
| T2.0–2.1 | factorisation; local enclosure (support hypothesis) | Lean: `Zsys_eq_sum`, `local_enclosure` | high (DLR, Leisink–Kappen) |
| T2.2/2.2′ | sharpness (closed and open interval); instance optimality | Lean: `sharp_closed`, `sharp_open_pos`, `instance_optimal`, `separator_enclosure_tightest` | medium–high |
| T2.3 | absolute vs relative separation; achievable set = (0,∞) × [min ρ, max ρ] | Lean: `separation`, `achievable_eq`, `freeEnergy_unbounded` | medium |
| T2.4 | nested hierarchy. **Corrected:** needs the inner support condition | Lean: `hierarchy_nested`; counterexample `hierarchy_needs_support` | high |
| T2.5/2.6 | finite termination; TV contraction | Lean | high |
| §2.3 / T4.4 | interval- and soft-exterior reductions; partial ensemble knowledge | **refuted as originally stated** (review M3); corrected with support conditions; unformalized | high |
| T2.7 | soft separator for the exterior class E_s (single-boundary-variable couplings), sharp | **hand proof only**; numerically sound; hypothesis shown necessary by an explicit counterexample (MRF agent) | medium |
| T2.8 | local certifiability ⟺ weak spatial mixing | hand proof (an elementary identity); numerically illustrated | high |
| T4.1–4.3 | certified affinity and selectivity; non-identifiability. **Corrected:** T4.2 needs k₁ > 0 | Lean: `affinity_mem_Icc`, `selectivity_mem_Icc`, `nonidentifiable*`; counterexample `selectivity_needs_positive_k1` | high (identity) / medium (certification form) |
| AF | multi-site proton binding is an antiferromagnetic 2-spin system in protonation variables | symbolic check | high for the sign structure; the computational-complexity link is medium |

Lean facts: branch `research/cec-lean-20261112`; 93 audited declarations;
only propext, Classical.choice and Quot.sound; no sorry or native_decide;
kernel-evaluated checker tests (`decide +kernel`). Details are in
`LEAN_FORMALIZATION.md`.

## 3. Counterexamples found (all preserved)

1. **T2.4 without inner support (formal).** Φ_A = (1, 0), Φ_B = (1, 1),
   K ≡ 1 gives ρ = 1 on T⁺ but ρ′ = 2.
2. **T4.2 without k₁ > 0 (formal).** k₁ = (1, 0), k₂ = (1, 1), uniform
   receptor gives S = 2 > max over the support of k₂/k₁ = 1.
3. **T1.5 without S₋ > 0 (formal).** Regions (L,U,r) = (0,1,1) and (1,1,0).
4. **T2.7 outside E_s (computational, MRF agent).** An exterior factor
   touching two boundary variables gives an exact ratio of 2/3, outside the
   soft enclosure [4/9, 5/9] but inside the hard enclosure [0, 1].
5. **Methodological lesson.** The coordinator's random searches used strictly
   positive weights and missed items 1–2. Degenerate (zero-weight) cases must
   be part of every counterexample search.

## 4. Positive computational evidence: separator certificates on MRFs

Branch `research/cec-mrf-20261112`; independent minimal checker; exact oracles
(brute force, variable elimination, row DP up to L = 18).
- **Size independence (H4) holds by construction, so it is not an empirical
  finding.** This was corrected after adversarial review M2. For L ≥ 2R+1 the
  radius-R local data are identical, so the enclosure cannot depend on L
  (widths 0.5, 0.335, 0.205, 0.106 for R = 0..3).
  - The original instance had Z_B/Z_A = 3/4 for all L by symmetry.
  - The spanning-tree "separate" baseline is a strawman for zero-field planar
    Ising, which is exactly solvable.
  - A follow-up with a bulk field (exact ratio varying with L) gives 0
    violations in 36 enclosures.
  - What the experiment does establish is the informativeness of local
    certificates at a given radius, together with checked soundness.
- **Phase-transition obstruction (T2.8).** Width at R = 10 is 0.0007 at
  β/β_c = 0.5. It is 0.38 at β_c, and plateaus at about 0.48–0.50 below T_c.
- **Sharpness** is realised exactly (12/12).
- **Fairness caveat.** The separate baseline is naive. Correlation-decay
  FPTAS give good separate bounds in uniqueness regimes at poly(n) cost. The
  supported claim is *equal-cost size-independence*, not that separate
  methods are hopeless.

## 5. Negative computational evidence: HP lattice mutations

Branch `research/cec-algo-20261112`; 71,859 oracle rows; 0 soundness
failures. The independent checker accepted 6,945 certificates.
- At matched frontiers for n ≥ 12, the coupled ΔΔF width ≈ the trivial ‖V‖∞
  width (median ratio 0.96–1.00). At n = 30 and 50 it is exactly trivial.
- Gains appear only after exploring roughly 35–60% of the tree, at n ≤ 12,
  where exact enumeration is cheap anyway.
- **Conclusion.** For lattice polymers a residue's environment is set by the
  global fold, there is no separator, and CEC adds a checkable combination
  rule but no new capability. This matches the coordinator's prior
  (`hp_ratio_floor.py`).

## 6. Real-structure application: certified model titration

Branch `research/cec-pka-20261112` @ f0bb318. Full details are in
`PKA_APPLICATION.md`. The protocol was frozen (fbaf401, 02:42 UTC) before any
experimental comparison, which was computed at 08:18 UTC.

**Model.** A rational antiferromagnetic protonation MRF: PROPKA 3.5.1 model
pKa values, Coulomb with ε = 4r, a 12 Å cutoff, and 15-digit rational
weights. Certificates are boundary-cube decision trees with T1 tests. A float
proposer builds them; an exact Fraction checker verifies them.

**Results.**
- *Soundness.* 10,000 certificate pairs on BPTI and four 20-site subsets,
  against brute-force exact marginals: **0 violations, 0 rejections.** On
  4HHB, 3,480 certificates against a variable-elimination reference: 0
  violations.
- *Capability within the model.* GroEL 1OEL has 1,064 sites and a min-fill
  treewidth bound of 31, so exact inference is infeasible. At pH 7 the mean
  protonation-fraction width is 0.434 / 0.305 / 0.165 / 0.099 for
  R = 0 / 4 / 6 / 8 Å, and the median at R = 8 is 1e-4. The independent
  checker accepted 19,482 certificates.
- **19% of sites (205/1064) do not converge** (width > 0.1 at R = 8). They
  sit in 57 small clusters (largest 7 sites). This is consistent with, but
  not proven to be, the strong-coupling / non-uniqueness signature; budget
  and interior-cap limits were not ruled out.
- *Monte Carlo*, under the pre-registered criterion (mean inside the
  interval ± 3 SE): only **55–65%** of (site, pH) pairs are inside. An earlier
  "100%" figure used an unregistered +1e-3 slack, which review P2-M1 flagged
  as a post-hoc deviation; it is now disclosed. The likely benign cause is
  frozen sites with zero MC standard error against tiny certified intervals.
  Certificate soundness is unaffected.
- **Model vs experiment (pre-registered): negative.** On PKAD-2 (9 proteins,
  134 matched sites), the certified-midpoint RMSE is 1.59 vs 0.74 for the
  null (model-pKa-only) predictor on the 75 sites with a single bracket
  (Δ = +0.84, CI [+0.34, +1.34]). The correlation of predicted vs
  experimental *shifts* is 0.02. The Coulomb-only model has no skill. The
  certificates are about the model, and the model is poor.
- **An assumption corrected.** 16/100 exactly analysed sites titrate
  *non-monotonically* in pH, at a 1e-12 threshold. The BPTI increases are
  physically negligible (≤ 3e-5); the largest increase, 0.013, is in an RNase
  subset. Coupled multi-site titration is known to produce
  irregular titration curves. A certified crossing of 1/2 therefore asserts
  that a crossing exists, not that it is unique; multiple crossings are
  flagged.

**Interpretation.** This is the clearest demonstration of what CEC adds:
rigorous, independently checked bounds on equilibrium quantities of a
real-structure-derived model far beyond the reach of exact inference, with
non-convergence localized to identifiable clusters. It is a capability for
*model computation*, not a pKa predictor, and the model's predictive value is
negative.

## 7. Independent adversarial review

Phase 1 covered the Lean, MRF, spec, prior-art and audit material. Details
are in `ADVERSARIAL_REVIEW.md`. **0 CRITICAL, 3 MAJOR, 8 MINOR.**
- **No checker accepted an unsound certificate.** This held across 716 Lean
  `#eval` cross-checks, 50 crafted MRF attacks, and 6,241 degenerate random
  certificates.
- All recomputed MRF numbers match to at least 12 digits.
- The headline Lean statements quantify as their docs claim.
- **MAJOR findings:**
  - M1: a forged axiom named `«Classical.choice»` passed both axiom gates.
    Fixed in the new gate; recorded for the old one.
  - M2: the MRF size-independence headline was tautological and its baseline
    a strawman. Reframed.
  - M3: the unformalized §2.3 and Theorem 4.4 reductions were false without
    support conditions. Corrected.

**Phase 2 (HP, pKa): 0 CRITICAL, 1 MAJOR, 2 MINOR.**
- Independent reproductions:
  - 40 HP oracle rows: exact ΔΔF matches to 2e-6 and lies inside every
    interval.
  - All 42,069 n ≤ 12 rows are sound, and the coupled width never exceeds
    trivial.
  - 1,200 BPTI exact records match to 1e-9 and lie inside the certified
    intervals.
- Both checkers rejected every crafted attack.
- The pre-registration order is confirmed in git.
- MAJOR P2-M1: the MC "100%" claim relied on an unregistered slack.
  Corrected and disclosed.
- Minor: the non-monotonicity magnitudes were overstated for BPTI
  (corrected), and malformed numeric input raised the wrong exception type.
  It already failed closed; it is now a CertificateError, with a regression
  test.

## 8. Assessment against the program's success levels

**Level 1 (foundational consolidation): achieved.**
- An evidence-based audit of `main`, #91 and #92–#95.
- Integration-level Lean is unverified in CI (all jobs time out).
- A `native_decide` result is mislabelled as kernel-checked.
- **Two soundness bugs in the existing axiom gate**, both reproduced, with
  fixes proposed.

**Level 2 (strong contribution in a restricted setting): partially
achieved.**
- A sharp theory of local relative-free-energy certificates, machine-checked
  in Lean (93 declarations, standard axioms only, verified by CI on GitHub),
  with executable checkers proved sound and exact.
- A demonstrated capability *within explicit models*: rigorous, independently
  checked equilibrium bounds for a 1,064-site protonation model (GroEL)
  beyond exact inference, with non-convergence localised to small clusters.
- Limits:
  - The underlying mathematics is largely classical (PRIOR_ART.md).
  - The MRF size-independence is a construction-level corollary.
  - The polymer (HP) case is a clear negative.
  - The real-structure model has no predictive skill against experiment.

**Levels 3 and 4: not achieved.** There is no evidence of a generalisable new
capability on independently evaluated scientific problems. The physically
meaningful applications are bounded by model realism, and the theory says
precisely why local certification fails in the regimes that matter most:
phase coexistence, global folding, and strongly coupled clusters.

## 9. Manuscript readiness (Deliverable G)

**Not publication-ready, and no manuscript is claimed.** A defensible
short methods paper could be built on "machine-checked separator certificates
for multi-site titration", but only after:
1. a realistic electrostatic model (e.g. continuum-electrostatics pairwise
   energies, desolvation, conformational flexibility as in MCCE) and a
   head-to-head comparison with an established pKa engine;
2. resolution or characterisation of the 19% non-converging sites (true
   non-uniqueness vs. budget limits);
3. full-text verification of the prior-art items left [unverified]: OSPREY
   BBK*/MARK* ratio handling, Leisink–Kappen sharpness statements, and
   Onufriev et al. on irregular titration;
4. landing the fixes for the existing axiom gate (owner approval).

A paper based on the current evidence would overstate both novelty and
relevance.
