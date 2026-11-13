# Coupled Ensemble Certificates: mathematical discovery report (Deliverable B)

**Status (2026-11-13):** draft. Sections marked *[pending]* await the pKa
application and the adversarial review. Nothing here is claimed as novel
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

*[pending: pKa agent]*

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

*[Phase 2 (HP, pKa): pending]*

## 8. Assessment against the program's success levels

*[to be finalised after §6–7]*. Provisional: **Level 1** (foundational
consolidation, audit, gate flaws found) **plus a partial Level 2**. The Level 2
part is a machine-checked, sharp, local relative-free-energy certificate
theory with a verified executable kernel, and a demonstrated size-independent
capability on MRF models. The mathematics is largely classical, and the
polymer case is negative. There is no evidence of Level 3 or 4.
