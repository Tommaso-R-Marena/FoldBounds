module

public import FoldBounds.CoupledRatio

/-!
# FoldBounds.CoupledSeparator — separator (Markov-blanket) certificates (CEC §2)

Finite types: boundary states `T`, interiors `IA`, `IB`, exterior states `O`.  Interior weights
`a : IA → T → ℝ`, `b : IB → T → ℝ` and a **shared** exterior weight `e : T → O → ℝ` give
`Z_X = Σ_{(i,τ,o)} x(i,τ) e(τ,o)` (`Zsys`).  With `Φ_X(τ) = Σ_i x(i,τ)` (`Phi`) and
`Ψ(τ) = Σ_o e(τ,o)` (`Psi`) we prove, following `docs/coupled/CEC_THEORY.md` §2:

* Lemma 2.0 `Zsys_eq_sum`: `Z_X = Σ_τ Φ_X(τ) Ψ(τ)`;
* Theorem 2.1 `local_enclosure`: `min_{T⁺} ρ ≤ Z_B/Z_A ≤ max_{T⁺} ρ` (support hypothesis);
* Theorem 2.2 `sharp_closed`, `sharp_open_pos`: every value of `[min ρ, max ρ]` is realised by a
  nonnegative exterior concentrated on two boundary states, every value of the open interval
  by a strictly positive exterior (only `T⁺ ≠ ∅` is needed, not `T⁺ = T`);
* Corollary 2.2′ `instance_optimal`, `separator_enclosure_tightest`;
* Theorem 2.3 `separation`, `achievable_eq`, `Z_unbounded`, `freeEnergy_unbounded`;
* Theorem 2.4 `push_mem_Icc`, `rho_push_eq`, `hierarchy_nested` (needs the support
  hypothesis; `hierarchy_needs_support` is a counterexample without it);
* Theorem 2.5 `finite_termination`;
* Theorem 2.6 `abs_sum_sub_le_tv_mul`, `hierarchy_pair_contraction`, `osc_contraction`;
* an executable checker `checkSeparator` with soundness (`checkSeparator_sound`,
  `checkSeparator_sound_sys`) and exactness (`checkSeparator_attained`).
-/

@[expose] public section

open Finset

namespace FoldBounds
namespace CoupledSeparator

variable {T O : Type*} [Fintype T] [Fintype O]

/-! ### Definitions and Lemma 2.0 -/

/-- `Z = Σ_{(i,τ,o)} x(i,τ) e(τ,o)`. -/
noncomputable def Zsys {I : Type*} [Fintype I] (x : I → T → ℝ) (e : T → O → ℝ) : ℝ :=
  ∑ p : I × T × O, x p.1 p.2.1 * e p.2.1 p.2.2

/-- Interior marginal `Φ(τ) = Σ_i x(i,τ)`. -/
noncomputable def Phi {I : Type*} [Fintype I] (x : I → T → ℝ) (τ : T) : ℝ := ∑ i, x i τ

/-- Exterior marginal `Ψ(τ) = Σ_o e(τ,o)`. -/
noncomputable def Psi (e : T → O → ℝ) (τ : T) : ℝ := ∑ o, e τ o

/-- Local ratio `ρ(τ) = Φ_B(τ) / Φ_A(τ)` (meaningful on `T⁺`). -/
noncomputable def rho (φA φB : T → ℝ) (τ : T) : ℝ := φB τ / φA τ

/-- Support `T⁺ = {τ : Φ_A(τ) > 0}`. -/
noncomputable def Tplus (φA : T → ℝ) : Finset T := univ.filter (fun τ => 0 < φA τ)

/-- `min_{T⁺} ρ`. -/
noncomputable def rhoMin (φA φB : T → ℝ) (h : (Tplus φA).Nonempty) : ℝ :=
  (Tplus φA).inf' h (rho φA φB)

/-- `max_{T⁺} ρ`. -/
noncomputable def rhoMax (φA φB : T → ℝ) (h : (Tplus φA).Nonempty) : ℝ :=
  (Tplus φA).sup' h (rho φA φB)

theorem mem_Tplus {φA : T → ℝ} {τ : T} : τ ∈ Tplus φA ↔ 0 < φA τ := by simp [Tplus]

theorem rhoMin_le {φA φB : T → ℝ} (h : (Tplus φA).Nonempty) {τ : T} (hτ : 0 < φA τ) :
    rhoMin φA φB h ≤ rho φA φB τ := inf'_le _ (mem_Tplus.mpr hτ)

theorem le_rhoMax {φA φB : T → ℝ} (h : (Tplus φA).Nonempty) {τ : T} (hτ : 0 < φA τ) :
    rho φA φB τ ≤ rhoMax φA φB h := le_sup' (rho φA φB) (mem_Tplus.mpr hτ)

theorem rhoMin_le_rhoMax {φA φB : T → ℝ} (h : (Tplus φA).Nonempty) :
    rhoMin φA φB h ≤ rhoMax φA φB h := by
  obtain ⟨τ, hτ⟩ := h
  exact (rhoMin_le _ (mem_Tplus.mp hτ)).trans (le_rhoMax _ (mem_Tplus.mp hτ))

/-- **Lemma 2.0 (factorisation).** `Σ_{(i,τ,o)} x(i,τ) e(τ,o) = Σ_τ Φ(τ) Ψ(τ)`. -/
theorem Zsys_eq_sum {I : Type*} [Fintype I] (x : I → T → ℝ) (e : T → O → ℝ) :
    Zsys x e = ∑ τ, Phi x τ * Psi e τ := by
  simp only [Zsys, Phi, Psi, Fintype.sum_prod_type, sum_mul_sum]
  exact sum_comm

/-! ### Theorem 2.1: local enclosure -/

/-- Core weighted-average lemma: with nonnegative weights `φA ψ` and the support condition,
`Σ φB ψ / Σ φA ψ` lies between any bounds on the ratios `φB/φA` over `{φA > 0}`. -/
theorem ratio_mem_Icc_of_local {φA φB ψ : T → ℝ} (hφA : ∀ τ, 0 ≤ φA τ) (hψ : ∀ τ, 0 ≤ ψ τ)
    (hsupp : ∀ τ, φA τ = 0 → φB τ = 0) (hZ : 0 < ∑ τ, φA τ * ψ τ) {m M : ℝ}
    (hm : ∀ τ, 0 < φA τ → m ≤ φB τ / φA τ) (hM : ∀ τ, 0 < φA τ → φB τ / φA τ ≤ M) :
    (∑ τ, φB τ * ψ τ) / (∑ τ, φA τ * ψ τ) ∈ Set.Icc m M := by
  constructor
  · rw [le_div_iff₀ hZ, mul_sum]
    refine sum_le_sum fun τ _ => ?_
    rcases (hφA τ).lt_or_eq with h | h
    · have h' := hm τ h
      rw [le_div_iff₀ h] at h'
      calc m * (φA τ * ψ τ) = m * φA τ * ψ τ := by ring
        _ ≤ φB τ * ψ τ := mul_le_mul_of_nonneg_right h' (hψ τ)
    · rw [← h, hsupp τ h.symm]
      simp
  · rw [div_le_iff₀ hZ, mul_sum]
    refine sum_le_sum fun τ _ => ?_
    rcases (hφA τ).lt_or_eq with h | h
    · have h' := hM τ h
      rw [div_le_iff₀ h] at h'
      calc φB τ * ψ τ ≤ M * φA τ * ψ τ := mul_le_mul_of_nonneg_right h' (hψ τ)
        _ = M * (φA τ * ψ τ) := by ring
    · rw [← h, hsupp τ h.symm]
      simp

variable {IA IB : Type*} [Fintype IA] [Fintype IB] {a : IA → T → ℝ} {b : IB → T → ℝ}

omit [Fintype T] in
theorem Phi_nonneg (ha : ∀ i τ, 0 ≤ a i τ) (τ : T) : 0 ≤ Phi a τ :=
  sum_nonneg fun i _ => ha i τ

omit [Fintype T] in
theorem Psi_nonneg {e : T → O → ℝ} (he : ∀ τ o, 0 ≤ e τ o) (τ : T) : 0 ≤ Psi e τ :=
  sum_nonneg fun o _ => he τ o

/-- `Z_A > 0` forces a nonempty support `T⁺`. -/
theorem Tplus_nonempty_of_Z_pos {e : T → O → ℝ} (ha : ∀ i τ, 0 ≤ a i τ)
    (hZ : 0 < Zsys a e) : (Tplus (Phi a)).Nonempty := by
  rw [Zsys_eq_sum] at hZ
  obtain ⟨τ, -, hτ⟩ := CoupledRatio.exists_pos_of_sum_pos hZ
  refine ⟨τ, mem_Tplus.mpr ((Phi_nonneg ha τ).lt_of_ne fun h => ?_)⟩
  rw [← h, zero_mul] at hτ
  exact lt_irrefl 0 hτ

/-- **Theorem 2.1 (local enclosure), general bounds.** Nonnegative interior and exterior
weights, `Z_A > 0`, and the support condition `Φ_A(τ) = 0 → Φ_B(τ) = 0`.  No sign condition
on `b` is needed in this form. -/
theorem local_enclosure_of_bounds {e : T → O → ℝ} (ha : ∀ i τ, 0 ≤ a i τ)
    (he : ∀ τ o, 0 ≤ e τ o) (hsupp : ∀ τ, Phi a τ = 0 → Phi b τ = 0) (hZ : 0 < Zsys a e)
    {m M : ℝ} (hm : ∀ τ, 0 < Phi a τ → m ≤ rho (Phi a) (Phi b) τ)
    (hM : ∀ τ, 0 < Phi a τ → rho (Phi a) (Phi b) τ ≤ M) :
    Zsys b e / Zsys a e ∈ Set.Icc m M := by
  rw [Zsys_eq_sum] at hZ ⊢
  rw [Zsys_eq_sum]
  exact ratio_mem_Icc_of_local (Phi_nonneg ha) (Psi_nonneg he) hsupp hZ hm hM

/-- **Theorem 2.1 (local enclosure).** `min_{T⁺} ρ ≤ Z_B/Z_A ≤ max_{T⁺} ρ`. -/
theorem local_enclosure {e : T → O → ℝ} (ha : ∀ i τ, 0 ≤ a i τ) (he : ∀ τ o, 0 ≤ e τ o)
    (hsupp : ∀ τ, Phi a τ = 0 → Phi b τ = 0) (hZ : 0 < Zsys a e)
    (hT : (Tplus (Phi a)).Nonempty) :
    Zsys b e / Zsys a e ∈ Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT) :=
  local_enclosure_of_bounds ha he hsupp hZ (fun _ h => rhoMin_le hT h) (fun _ h => le_rhoMax hT h)

/-- Theorem 2.1 with the spec's literal hypotheses: `b ≥ 0` and no `τ` with
`Φ_A(τ) = 0 < Φ_B(τ)`. -/
theorem local_enclosure_spec {e : T → O → ℝ} (ha : ∀ i τ, 0 ≤ a i τ) (hb : ∀ i τ, 0 ≤ b i τ)
    (he : ∀ τ o, 0 ≤ e τ o) (hsupp : ¬ ∃ τ, Phi a τ = 0 ∧ 0 < Phi b τ) (hZ : 0 < Zsys a e) :
    ∃ hT : (Tplus (Phi a)).Nonempty,
      Zsys b e / Zsys a e ∈ Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT) := by
  have hT := Tplus_nonempty_of_Z_pos ha hZ
  refine ⟨hT, local_enclosure ha he (fun τ h => ?_) hZ hT⟩
  exact le_antisymm (not_lt.mp fun h' => hsupp ⟨τ, h, h'⟩) (Phi_nonneg hb τ)

/-! ### Realising exterior marginals -/

/-- Spread a boundary weight `ψ` uniformly over the exterior states. -/
noncomputable def spread [Nonempty O] (ψ : T → ℝ) : T → O → ℝ :=
  fun τ _ => ψ τ / Fintype.card O

omit [Fintype T] in
theorem Psi_spread [Nonempty O] (ψ : T → ℝ) (τ : T) : Psi (spread (O := O) ψ) τ = ψ τ := by
  have hc : (Fintype.card O : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  simp only [Psi, spread, sum_const, card_univ, nsmul_eq_mul]
  field_simp

theorem Zsys_spread [Nonempty O] {I : Type*} [Fintype I] (x : I → T → ℝ) (ψ : T → ℝ) :
    Zsys x (spread (O := O) ψ) = ∑ τ, Phi x τ * ψ τ := by
  rw [Zsys_eq_sum]
  simp only [Psi_spread]

omit [Fintype T] in
theorem spread_nonneg [Nonempty O] {ψ : T → ℝ} (hψ : ∀ τ, 0 ≤ ψ τ) (τ : T) (o : O) :
    0 ≤ spread ψ τ o := div_nonneg (hψ τ) (Nat.cast_nonneg _)

omit [Fintype T] in
theorem spread_pos [Nonempty O] {ψ : T → ℝ} (hψ : ∀ τ, 0 < ψ τ) (τ : T) (o : O) :
    0 < spread ψ τ o := div_pos (hψ τ) (Nat.cast_pos.mpr Fintype.card_pos)

/-! ### Theorem 2.2: sharpness -/

/-- Two-point boundary weight realising any value between two local ratios with `Σ φA ψ = 1`. -/
theorem exists_two_point {φA φB : T → ℝ} {τ₁ τ₂ : T} (h₁ : 0 < φA τ₁) (h₂ : 0 < φA τ₂)
    {lam : ℝ} (hlam : lam ∈ segment ℝ (φB τ₁ / φA τ₁) (φB τ₂ / φA τ₂)) :
    ∃ ψ : T → ℝ, (∀ τ, 0 ≤ ψ τ) ∧ (∀ τ, τ ≠ τ₁ → τ ≠ τ₂ → ψ τ = 0) ∧
      ∑ τ, φA τ * ψ τ = 1 ∧ ∑ τ, φB τ * ψ τ = lam := by
  classical
  obtain ⟨θ₁, θ₂, hθ₁, hθ₂, hθ, hcomb⟩ := hlam
  refine ⟨fun τ => (if τ = τ₁ then θ₁ / φA τ₁ else 0) + (if τ = τ₂ then θ₂ / φA τ₂ else 0),
    fun τ => add_nonneg (by split_ifs <;> positivity) (by split_ifs <;> positivity),
    fun τ hn₁ hn₂ => by simp [hn₁, hn₂], ?_, ?_⟩
  · simp only [mul_add, sum_add_distrib, mul_ite, mul_zero, sum_ite_eq', mem_univ, if_true]
    rw [mul_div_cancel₀ _ h₁.ne', mul_div_cancel₀ _ h₂.ne', hθ]
  · simp only [mul_add, sum_add_distrib, mul_ite, mul_zero, sum_ite_eq', mem_univ, if_true]
    rw [← hcomb, smul_eq_mul, smul_eq_mul]
    ring

/-- Every value of `[min ρ, max ρ]` is `Σ φB ψ` for a two-point `ψ ≥ 0` with `Σ φA ψ = 1`. -/
theorem exists_weight_of_mem_Icc {φA φB : T → ℝ} (hT : (Tplus φA).Nonempty) {lam : ℝ}
    (hlam : lam ∈ Set.Icc (rhoMin φA φB hT) (rhoMax φA φB hT)) :
    ∃ τ₁ τ₂ : T, ∃ ψ : T → ℝ, (∀ τ, 0 ≤ ψ τ) ∧ (∀ τ, τ ≠ τ₁ → τ ≠ τ₂ → ψ τ = 0) ∧
      ∑ τ, φA τ * ψ τ = 1 ∧ ∑ τ, φB τ * ψ τ = lam := by
  obtain ⟨τ₁, h₁, he₁⟩ := exists_mem_eq_inf' hT (rho φA φB)
  obtain ⟨τ₂, h₂, he₂⟩ := exists_mem_eq_sup' hT (rho φA φB)
  have hl : rho φA φB τ₁ ≤ lam := by rw [← he₁]; exact hlam.1
  have hu : lam ≤ rho φA φB τ₂ := by rw [← he₂]; exact hlam.2
  have hseg : lam ∈ segment ℝ (rho φA φB τ₁) (rho φA φB τ₂) := by
    rw [segment_eq_Icc (hl.trans hu)]
    exact ⟨hl, hu⟩
  exact ⟨τ₁, τ₂, exists_two_point (mem_Tplus.mp h₁) (mem_Tplus.mp h₂) hseg⟩

/-- **Theorem 2.2 (sharpness, closed interval).** For local data with `T⁺ ≠ ∅` (no sign or
support hypothesis is needed) and any nonempty exterior state space `O`, every
`λ ∈ [min ρ, max ρ]` is realised by a nonnegative exterior with `Z_A = 1` and `Z_B = λ`; the
exterior is supported on at most two boundary states. -/
theorem sharp_closed [Nonempty O] (hT : (Tplus (Phi a)).Nonempty) {lam : ℝ}
    (hlam : lam ∈ Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT)) :
    ∃ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) ∧
      (∃ τ₁ τ₂ : T, ∀ τ o, τ ≠ τ₁ → τ ≠ τ₂ → e τ o = 0) ∧
      Zsys a e = 1 ∧ Zsys b e / Zsys a e = lam := by
  obtain ⟨τ₁, τ₂, ψ, hψ, hsupp, hA, hB⟩ := exists_weight_of_mem_Icc hT hlam
  refine ⟨spread ψ, spread_nonneg hψ, ⟨τ₁, τ₂, fun τ o h₁ h₂ => ?_⟩, ?_, ?_⟩
  · simp [spread, hsupp τ h₁ h₂]
  · rw [Zsys_spread, hA]
  · rw [Zsys_spread, Zsys_spread, hA, hB, div_one]

/-- Positive boundary weight realising any value strictly between two local ratios. -/
theorem exists_pos_weight {φA φB : T → ℝ} (hφA : ∀ τ, 0 ≤ φA τ) {τm τp : T}
    (hm : 0 < φA τm) (hp : 0 < φA τp) {lam : ℝ} (hlo : φB τm / φA τm < lam)
    (hhi : lam < φB τp / φA τp) :
    ∃ ψ : T → ℝ, (∀ τ, 0 < ψ τ) ∧ 0 < ∑ τ, φA τ * ψ τ ∧
      (∑ τ, φB τ * ψ τ) / (∑ τ, φA τ * ψ τ) = lam := by
  classical
  rw [div_lt_iff₀ hm] at hlo
  rw [lt_div_iff₀ hp] at hhi
  set D := lam * ∑ τ, φA τ - ∑ τ, φB τ with hD
  set cp := φB τp - lam * φA τp with hcp
  set cm := lam * φA τm - φB τm with hcm
  have hcp0 : 0 < cp := by rw [hcp]; linarith
  have hcm0 : 0 < cm := by rw [hcm]; linarith
  set sp := max D 0 / cp with hsp
  set sm := max (-D) 0 / cm with hsm
  have hsp0 : 0 ≤ sp := div_nonneg (le_max_right _ _) hcp0.le
  have hsm0 : 0 ≤ sm := div_nonneg (le_max_right _ _) hcm0.le
  set ψ : T → ℝ := fun τ => 1 + (if τ = τp then sp else 0) + (if τ = τm then sm else 0) with hψ
  have hψpos : ∀ τ, 0 < ψ τ := fun τ => by
    simp only [hψ]
    have h1 : 0 ≤ (if τ = τp then sp else 0) := by split_ifs <;> linarith
    have h2 : 0 ≤ (if τ = τm then sm else 0) := by split_ifs <;> linarith
    linarith
  have hsum : ∀ φ : T → ℝ, ∑ τ, φ τ * ψ τ = ∑ τ, φ τ + φ τp * sp + φ τm * sm := fun φ => by
    simp only [hψ, mul_add, mul_one, sum_add_distrib, mul_ite, mul_zero, sum_ite_eq', mem_univ,
      if_true]
  have hZ : 0 < ∑ τ, φA τ * ψ τ :=
    sum_pos' (fun τ _ => mul_nonneg (hφA τ) (hψpos τ).le)
      ⟨τm, mem_univ _, mul_pos hm (hψpos τm)⟩
  refine ⟨ψ, hψpos, hZ, ?_⟩
  rw [div_eq_iff hZ.ne', hsum, hsum]
  have e1 : cp * sp = max D 0 := mul_div_cancel₀ _ hcp0.ne'
  have e2 : cm * sm = max (-D) 0 := mul_div_cancel₀ _ hcm0.ne'
  have e3 : max D 0 - max (-D) 0 = D := by
    rcases le_total D 0 with h | h
    · rw [max_eq_right h, max_eq_left (by linarith)]; ring
    · rw [max_eq_left h, max_eq_right (by linarith)]; ring
  rw [hcp] at e1
  rw [hcm] at e2
  rw [hD] at e3
  linear_combination e1 - e2 + e3

/-- **Theorem 2.2 (sharpness, open interval, strictly positive exteriors).** Needs only
`Φ_A ≥ 0` (through `a ≥ 0`) and `T⁺ ≠ ∅`. -/
theorem sharp_open_pos [Nonempty O] (ha : ∀ i τ, 0 ≤ a i τ) (hT : (Tplus (Phi a)).Nonempty)
    {lam : ℝ} (hlam : lam ∈ Set.Ioo (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT)) :
    ∃ e : T → O → ℝ, (∀ τ o, 0 < e τ o) ∧ 0 < Zsys a e ∧ Zsys b e / Zsys a e = lam := by
  obtain ⟨τm, hm, hem⟩ := exists_mem_eq_inf' hT (rho (Phi a) (Phi b))
  obtain ⟨τp, hp, hep⟩ := exists_mem_eq_sup' hT (rho (Phi a) (Phi b))
  have hl : rho (Phi a) (Phi b) τm < lam := by rw [← hem]; exact hlam.1
  have hu : lam < rho (Phi a) (Phi b) τp := by rw [← hep]; exact hlam.2
  obtain ⟨ψ, hψ, hZ, hq⟩ := exists_pos_weight (Phi_nonneg ha) (mem_Tplus.mp hm)
    (mem_Tplus.mp hp) hl hu
  exact ⟨spread ψ, spread_pos hψ, by rw [Zsys_spread]; exact hZ,
    by rw [Zsys_spread, Zsys_spread]; exact hq⟩

/-! ### Corollary 2.2′: instance optimality -/

/-- **Corollary 2.2′.** Any interval `[l, u]` that is sound for every nonnegative exterior on
some nonempty exterior space `O` contains `[min ρ, max ρ]`. -/
theorem instance_optimal [Nonempty O] (hT : (Tplus (Phi a)).Nonempty) {l u : ℝ}
    (hsound : ∀ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) → 0 < Zsys a e →
      Zsys b e / Zsys a e ∈ Set.Icc l u) :
    Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT) ⊆ Set.Icc l u := by
  intro lam hlam
  obtain ⟨e, he, -, hA, hq⟩ := sharp_closed (O := O) hT hlam
  exact hq ▸ hsound e he (by rw [hA]; exact one_pos)

/-- **Corollary 2.2′ (procedure form).** A procedure whose output depends only on the local
data and which is sound for every exterior outputs an interval containing `[min ρ, max ρ]`. -/
theorem local_procedure_contains [Nonempty O]
    (proc : (IA → T → ℝ) → (IB → T → ℝ) → ℝ × ℝ) (hT : (Tplus (Phi a)).Nonempty)
    (hsound : ∀ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) → 0 < Zsys a e →
      Zsys b e / Zsys a e ∈ Set.Icc (proc a b).1 (proc a b).2) :
    Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT) ⊆
      Set.Icc (proc a b).1 (proc a b).2 :=
  instance_optimal hT hsound

/-- The separator enclosure is the tightest sound local interval: it is sound (Theorem 2.1)
and contained in every sound interval (Corollary 2.2′). -/
theorem separator_enclosure_tightest [Nonempty O] (ha : ∀ i τ, 0 ≤ a i τ)
    (hsupp : ∀ τ, Phi a τ = 0 → Phi b τ = 0) (hT : (Tplus (Phi a)).Nonempty) :
    (∀ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) → 0 < Zsys a e →
      Zsys b e / Zsys a e ∈ Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT)) ∧
    ∀ l u : ℝ, (∀ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) → 0 < Zsys a e →
      Zsys b e / Zsys a e ∈ Set.Icc l u) →
      Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT) ⊆ Set.Icc l u :=
  ⟨fun _ he hZ => local_enclosure ha he hsupp hZ hT, fun _ _ h => instance_optimal hT h⟩

/-! ### Theorem 2.3: separation of absolute and relative information -/

/-- **Theorem 2.3.** For every `z > 0` and `λ ∈ [min ρ, max ρ]` some nonnegative exterior has
`Z_A = z` and `Z_B / Z_A = λ`. -/
theorem separation [Nonempty O] (hT : (Tplus (Phi a)).Nonempty) {z lam : ℝ} (hz : 0 < z)
    (hlam : lam ∈ Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT)) :
    ∃ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) ∧ Zsys a e = z ∧ Zsys b e / Zsys a e = lam := by
  obtain ⟨τ₁, τ₂, ψ, hψ, -, hA, hB⟩ := exists_weight_of_mem_Icc hT hlam
  refine ⟨spread (fun τ => z * ψ τ), spread_nonneg fun τ => mul_nonneg hz.le (hψ τ), ?_, ?_⟩
  · rw [Zsys_spread]
    simp only [mul_left_comm _ z, ← mul_sum, hA, mul_one]
  · rw [Zsys_spread, Zsys_spread]
    simp only [mul_left_comm _ z, ← mul_sum, hA, hB, mul_one]
    field_simp

/-- **Theorem 2.3, exact joint image.** Under Theorem 2.1's hypotheses on the local data, the
set of pairs `(Z_A, Z_B/Z_A)` realised by nonnegative exteriors with `Z_A > 0` is exactly
`(0, ∞) × [min ρ, max ρ]`. -/
theorem achievable_eq [Nonempty O] (ha : ∀ i τ, 0 ≤ a i τ)
    (hsupp : ∀ τ, Phi a τ = 0 → Phi b τ = 0) (hT : (Tplus (Phi a)).Nonempty) :
    {p : ℝ × ℝ | ∃ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) ∧ 0 < Zsys a e ∧
        p = (Zsys a e, Zsys b e / Zsys a e)} =
      Set.Ioi 0 ×ˢ Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT) := by
  ext ⟨z, lam⟩
  constructor
  · rintro ⟨e, he, hZ, hp⟩
    simp only [Prod.mk.injEq] at hp
    obtain ⟨rfl, rfl⟩ := hp
    exact ⟨hZ, local_enclosure ha he hsupp hZ hT⟩
  · rintro ⟨hz, hlam⟩
    obtain ⟨e, he, hA, hq⟩ := separation (O := O) hT (Set.mem_Ioi.mp hz) hlam
    refine ⟨e, he, by rw [hA]; exact hz, ?_⟩
    simp only [Prod.mk.injEq]
    exact ⟨hA.symm, hq.symm⟩

/-- **Theorem 2.3 (i).** No finite upper bound and no positive lower bound on `Z_A` holds for
all exteriors: local data cannot bound `Z_A`. -/
theorem Z_unbounded [Nonempty O] (hT : (Tplus (Phi a)).Nonempty) (z : ℝ) (hz : 0 < z) :
    ∃ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) ∧ Zsys a e = z := by
  obtain ⟨e, he, hA, -⟩ := separation (b := a) (O := O) hT hz
    ⟨le_refl _, rhoMin_le_rhoMax hT⟩
  exact ⟨e, he, hA⟩

/-- **Theorem 2.3 (i), free energy.** `F_A = -log Z_A` is unbounded in both directions over
the exteriors compatible with fixed local data. -/
theorem freeEnergy_unbounded [Nonempty O] (hT : (Tplus (Phi a)).Nonempty) (F : ℝ) :
    ∃ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) ∧ 0 < Zsys a e ∧ -Real.log (Zsys a e) = F := by
  obtain ⟨e, he, hA⟩ := Z_unbounded (O := O) hT (Real.exp (-F)) (Real.exp_pos _)
  exact ⟨e, he, hA ▸ Real.exp_pos _, by rw [hA, Real.log_exp, neg_neg]⟩

/-! ### Theorem 2.4: nested hierarchy -/

section Hierarchy

variable {T' : Type*} [Fintype T']

/-- Pushforward `Φ'(τ') = Σ_τ Φ(τ) K(τ, τ')` through a shared kernel. -/
noncomputable def push (φ : T → ℝ) (K : T → T' → ℝ) (τ' : T') : ℝ := ∑ τ, φ τ * K τ τ'

/-- Mixing weights `κ(τ'; τ) = Φ_A(τ) K(τ,τ') / Φ'_A(τ')`. -/
noncomputable def kappa (φA : T → ℝ) (K : T → T' → ℝ) (τ' : T') (τ : T) : ℝ :=
  φA τ * K τ τ' / push φA K τ'

omit [Fintype T'] in
/-- Absorbing an annulus with a shared factor `k` produces the kernel form. -/
theorem Phi_absorb {I W : Type*} [Fintype I] [Fintype W] (x : I → T → ℝ)
    (k : T → W → T' → ℝ) (τ' : T') :
    Phi (fun (p : I × T × W) (τ' : T') => x p.1 p.2.1 * k p.2.1 p.2.2 τ') τ' =
      push (Phi x) (fun τ τ' => ∑ w, k τ w τ') τ' :=
  Zsys_eq_sum x (fun τ w => k τ w τ')

variable {φA φB : T → ℝ} {K : T → T' → ℝ}

omit [Fintype T'] in
theorem kappa_nonneg (hφA : ∀ τ, 0 ≤ φA τ) (hK : ∀ τ τ', 0 ≤ K τ τ') (τ' : T') (τ : T) :
    0 ≤ kappa φA K τ' τ :=
  div_nonneg (mul_nonneg (hφA τ) (hK τ τ')) (sum_nonneg fun τ _ => mul_nonneg (hφA τ) (hK τ τ'))

omit [Fintype T'] in
theorem kappa_sum_eq_one {τ' : T'} (hpos : 0 < push φA K τ') : ∑ τ, kappa φA K τ' τ = 1 := by
  simp only [kappa]
  rw [← sum_div]
  exact div_self hpos.ne'

omit [Fintype T'] in
theorem kappa_eq_zero_of_not_pos (hφA : ∀ τ, 0 ≤ φA τ) {τ : T} (h : ¬ 0 < φA τ) (τ' : T') :
    kappa φA K τ' τ = 0 := by
  rw [kappa, le_antisymm (not_lt.mp h) (hφA τ), zero_mul, zero_div]

omit [Fintype T'] in
/-- The support condition propagates to the larger separator. -/
theorem push_supp (hφA : ∀ τ, 0 ≤ φA τ) (hK : ∀ τ τ', 0 ≤ K τ τ')
    (hsupp : ∀ τ, φA τ = 0 → φB τ = 0) (τ' : T') (h0 : push φA K τ' = 0) :
    push φB K τ' = 0 := by
  have hz := (sum_eq_zero_iff_of_nonneg fun τ _ => mul_nonneg (hφA τ) (hK τ τ')).mp h0
  refine sum_eq_zero fun τ hτ => ?_
  rcases mul_eq_zero.mp (hz τ hτ) with h | h
  · rw [hsupp τ h, zero_mul]
  · rw [h, mul_zero]

omit [Fintype T'] in
/-- **Theorem 2.4 (mixture identity).** `ρ'(τ') = Σ_τ κ(τ'; τ) ρ(τ)`. -/
theorem rho_push_eq (hsupp : ∀ τ, φA τ = 0 → φB τ = 0) {τ' : T'} (hpos : 0 < push φA K τ') :
    rho (push φA K) (push φB K) τ' = ∑ τ, kappa φA K τ' τ * rho φA φB τ := by
  rw [rho, push, sum_div]
  refine sum_congr rfl fun τ _ => ?_
  rw [kappa, rho]
  by_cases h : φA τ = 0
  · rw [hsupp τ h, h]; simp
  · field_simp

omit [Fintype T'] in
/-- **Theorem 2.4 (enclosure).** Every ratio of the larger separator lies in any interval that
bounds the ratios of the smaller one on its support. -/
theorem push_mem_Icc (hφA : ∀ τ, 0 ≤ φA τ) (hK : ∀ τ τ', 0 ≤ K τ τ')
    (hsupp : ∀ τ, φA τ = 0 → φB τ = 0) {m M : ℝ} (hm : ∀ τ, 0 < φA τ → m ≤ rho φA φB τ)
    (hM : ∀ τ, 0 < φA τ → rho φA φB τ ≤ M) {τ' : T'} (hpos : 0 < push φA K τ') :
    rho (push φA K) (push φB K) τ' ∈ Set.Icc m M :=
  ratio_mem_Icc_of_local hφA (fun τ => hK τ τ') hsupp hpos hm hM

/-- **Theorem 2.4 (nesting).** `[min ρ', max ρ'] ⊆ [min ρ, max ρ]`. -/
theorem hierarchy_nested (hφA : ∀ τ, 0 ≤ φA τ) (hK : ∀ τ τ', 0 ≤ K τ τ')
    (hsupp : ∀ τ, φA τ = 0 → φB τ = 0) (hT : (Tplus φA).Nonempty)
    (hT' : (Tplus (push φA K)).Nonempty) :
    Set.Icc (rhoMin (push φA K) (push φB K) hT') (rhoMax (push φA K) (push φB K) hT') ⊆
      Set.Icc (rhoMin φA φB hT) (rhoMax φA φB hT) := by
  have hmem := fun τ' (h : τ' ∈ Tplus (push φA K)) => push_mem_Icc hφA hK hsupp
    (fun _ h => rhoMin_le hT h) (fun _ h => le_rhoMax hT h) (mem_Tplus.mp h)
  exact Set.Icc_subset_Icc (le_inf' _ _ fun τ' h => (hmem τ' h).1)
    (sup'_le _ _ fun τ' h => (hmem τ' h).2)

/-- **Spec discrepancy for Theorem 2.4.** Without the support condition the nesting fails:
`Φ_A = (1, 0)`, `Φ_B = (1, 1)`, `K ≡ 1` onto a single new state gives `ρ' = 2` while
`ρ = 1` on `T⁺`. -/
theorem hierarchy_needs_support :
    ∃ (φA φB : Bool → ℝ) (K : Bool → Unit → ℝ), (∀ τ, 0 ≤ φA τ) ∧ (∀ τ, 0 ≤ φB τ) ∧
      (∀ τ τ', 0 ≤ K τ τ') ∧ (∀ τ, 0 < φA τ → rho φA φB τ = 1) ∧
      0 < push φA K () ∧ rho (push φA K) (push φB K) () = 2 := by
  refine ⟨fun x => if x then 1 else 0, fun _ => 1, fun _ _ => 1, fun τ => by cases τ <;> simp,
    fun _ => zero_le_one, fun _ _ => zero_le_one, fun τ h => ?_, ?_, ?_⟩
  · cases τ <;> simp_all [rho]
  · simp [push]
  · norm_num [rho, push]

end Hierarchy

/-! ### Theorem 2.5: finite termination -/

/-- **Theorem 2.5.** When the separator is a single boundary state, the ratio is exact for
every exterior with `Z_A > 0`, and the enclosure collapses to it.  No sign hypothesis is
needed. -/
theorem finite_termination [Unique T] {e : T → O → ℝ} (hZ : 0 < Zsys a e)
    (hT : (Tplus (Phi a)).Nonempty) :
    Zsys b e / Zsys a e = rho (Phi a) (Phi b) default ∧
      rhoMin (Phi a) (Phi b) hT = Zsys b e / Zsys a e ∧
      rhoMax (Phi a) (Phi b) hT = Zsys b e / Zsys a e := by
  have hA : Zsys a e = Phi a default * Psi e default := by
    rw [Zsys_eq_sum, Fintype.sum_unique]
  have hB : Zsys b e = Phi b default * Psi e default := by
    rw [Zsys_eq_sum, Fintype.sum_unique]
  have hψ : Psi e default ≠ 0 := fun h => by rw [hA, h, mul_zero] at hZ; exact lt_irrefl 0 hZ
  have hq : Zsys b e / Zsys a e = rho (Phi a) (Phi b) default := by
    rw [hA, hB, rho, mul_div_mul_right _ _ hψ]
  have hall : ∀ τ : T, rho (Phi a) (Phi b) τ = rho (Phi a) (Phi b) default := fun τ => by
    rw [Unique.eq_default τ]
  obtain ⟨τ₀, hτ₀⟩ := hT
  refine ⟨hq, ?_, ?_⟩
  · rw [hq]
    exact le_antisymm ((inf'_le _ hτ₀).trans_eq (hall τ₀))
      (le_inf' _ _ fun τ _ => (hall τ).ge)
  · rw [hq]
    exact le_antisymm (sup'_le _ _ fun τ _ => (hall τ).le)
      ((hall τ₀).symm.trans_le (le_sup' _ hτ₀))

/-! ### Theorem 2.6: a-posteriori contraction -/

omit [Fintype T] in
/-- **Theorem 2.6, key inequality.** For weight vectors with equal total mass (in particular
probability vectors) and `m ≤ f ≤ M` on `s`:
`|p·f - p'·f| ≤ TV(p, p') · (M - m)` with `TV(p,p') = ½ Σ |p - p'|`. -/
theorem abs_sum_sub_le_tv_mul {s : Finset T} (p p' f : T → ℝ)
    (hsum : ∑ τ ∈ s, p τ = ∑ τ ∈ s, p' τ) {m M : ℝ} (hm : ∀ τ ∈ s, m ≤ f τ)
    (hM : ∀ τ ∈ s, f τ ≤ M) :
    |∑ τ ∈ s, p τ * f τ - ∑ τ ∈ s, p' τ * f τ| ≤ (1 / 2 * ∑ τ ∈ s, |p τ - p' τ|) * (M - m) := by
  have key : ∑ τ ∈ s, (p τ - p' τ) * (f τ - (m + M) / 2) =
      (∑ τ ∈ s, p τ * f τ - ∑ τ ∈ s, p' τ * f τ) -
        (m + M) / 2 * (∑ τ ∈ s, p τ - ∑ τ ∈ s, p' τ) := by
    rw [mul_sub, mul_sum, mul_sum, ← sum_sub_distrib, ← sum_sub_distrib, ← sum_sub_distrib]
    exact sum_congr rfl fun τ _ => by ring
  rw [hsum, sub_self, mul_zero, sub_zero] at key
  rw [← key]
  calc |∑ τ ∈ s, (p τ - p' τ) * (f τ - (m + M) / 2)|
      ≤ ∑ τ ∈ s, |(p τ - p' τ) * (f τ - (m + M) / 2)| := abs_sum_le_sum_abs _ _
    _ = ∑ τ ∈ s, |p τ - p' τ| * |f τ - (m + M) / 2| := by simp only [abs_mul]
    _ ≤ ∑ τ ∈ s, |p τ - p' τ| * ((M - m) / 2) := sum_le_sum fun τ hτ =>
        mul_le_mul_of_nonneg_left
          (abs_le.mpr ⟨by linarith [hm τ hτ], by linarith [hM τ hτ]⟩) (abs_nonneg _)
    _ = (1 / 2 * ∑ τ ∈ s, |p τ - p' τ|) * (M - m) := by rw [← sum_mul]; ring

section Contraction

variable {T' : Type*} [Fintype T'] {φA φB : T → ℝ} {K : T → T' → ℝ}

omit [Fintype T'] in
/-- **Theorem 2.6 (pairwise Dobrushin bound).** For two states of the larger separator,
`|ρ'(τ'₁) - ρ'(τ'₂)| ≤ TV(κ(τ'₁;·), κ(τ'₂;·)) · (M - m)` whenever `m ≤ ρ ≤ M` on `T⁺`. -/
theorem hierarchy_pair_contraction (hφA : ∀ τ, 0 ≤ φA τ) (hsupp : ∀ τ, φA τ = 0 → φB τ = 0)
    {m M : ℝ} (hm : ∀ τ, 0 < φA τ → m ≤ rho φA φB τ) (hM : ∀ τ, 0 < φA τ → rho φA φB τ ≤ M)
    {τ₁ τ₂ : T'} (h₁ : 0 < push φA K τ₁) (h₂ : 0 < push φA K τ₂) :
    |rho (push φA K) (push φB K) τ₁ - rho (push φA K) (push φB K) τ₂| ≤
      (1 / 2 * ∑ τ, |kappa φA K τ₁ τ - kappa φA K τ₂ τ|) * (M - m) := by
  classical
  -- replace the junk values of `ρ` off the support by `m`
  set f : T → ℝ := fun τ => if 0 < φA τ then rho φA φB τ else m with hf
  have hcongr : ∀ τ', ∑ τ, kappa φA K τ' τ * rho φA φB τ = ∑ τ, kappa φA K τ' τ * f τ :=
    fun τ' => sum_congr rfl fun τ _ => by
      by_cases h : 0 < φA τ
      · simp [hf, h]
      · rw [kappa_eq_zero_of_not_pos hφA h, zero_mul, zero_mul]
  rw [rho_push_eq hsupp h₁, rho_push_eq hsupp h₂, hcongr, hcongr]
  refine abs_sum_sub_le_tv_mul _ _ f
    ((kappa_sum_eq_one h₁).trans (kappa_sum_eq_one h₂).symm) (fun τ _ => ?_) (fun τ _ => ?_)
  · by_cases h : 0 < φA τ
    · simp only [hf, if_pos h]; exact hm τ h
    · simp only [hf, if_neg h]; exact le_rfl
  · have hmM : m ≤ M := by
      obtain ⟨τ₀, -, hτ₀⟩ := CoupledRatio.exists_pos_of_sum_pos h₁
      have hτ0 : 0 < φA τ₀ := by
        rcases (hφA τ₀).lt_or_eq with h | h
        · exact h
        · rw [← h, zero_mul] at hτ₀; exact absurd hτ₀ (lt_irrefl 0)
      exact (hm τ₀ hτ0).trans (hM τ₀ hτ0)
    by_cases h : 0 < φA τ
    · simp only [hf, if_pos h]; exact hM τ h
    · simp only [hf, if_neg h]; exact hmM

/-- **Theorem 2.6 (oscillation contraction).** If `δ` bounds the total-variation distance
between all pairs of mixing rows on the new support (the Dobrushin coefficient), then
`osc(ρ') ≤ δ · osc(ρ)`. -/
theorem osc_contraction (hφA : ∀ τ, 0 ≤ φA τ) (hsupp : ∀ τ, φA τ = 0 → φB τ = 0)
    (hT : (Tplus φA).Nonempty) (hT' : (Tplus (push φA K)).Nonempty) {δ : ℝ}
    (hδ : ∀ τ₁ ∈ Tplus (push φA K), ∀ τ₂ ∈ Tplus (push φA K),
      1 / 2 * ∑ τ, |kappa φA K τ₁ τ - kappa φA K τ₂ τ| ≤ δ) :
    rhoMax (push φA K) (push φB K) hT' - rhoMin (push φA K) (push φB K) hT' ≤
      δ * (rhoMax φA φB hT - rhoMin φA φB hT) := by
  obtain ⟨τm, hm, hem⟩ := exists_mem_eq_inf' hT' (rho (push φA K) (push φB K))
  obtain ⟨τp, hp, hep⟩ := exists_mem_eq_sup' hT' (rho (push φA K) (push φB K))
  rw [rhoMin, rhoMax, hem, hep]
  have h := hierarchy_pair_contraction (K := K) hφA hsupp (fun _ h => rhoMin_le hT h)
    (fun _ h => le_rhoMax hT h) (mem_Tplus.mp hp) (mem_Tplus.mp hm)
  calc rho (push φA K) (push φB K) τp - rho (push φA K) (push φB K) τm
      ≤ |rho (push φA K) (push φB K) τp - rho (push φA K) (push φB K) τm| := le_abs_self _
    _ ≤ _ := h
    _ ≤ δ * (rhoMax φA φB hT - rhoMin φA φB hT) :=
        mul_le_mul_of_nonneg_right (hδ τp hp τm hm) (sub_nonneg.mpr (rhoMin_le_rhoMax hT))

end Contraction

/-! ### Executable separator checker -/

section Executable

/-- Insert a ratio into an optional interval. -/
def includeRatio (q : ℚ) : Option (ℚ × ℚ) → Option (ℚ × ℚ)
  | none => some (q, q)
  | some (lo, hi) => some (min q lo, max q hi)

/-- `[min ρ, max ρ]` over the entries `(Φ_A(τ), Φ_B(τ))` with `Φ_A(τ) > 0`. -/
def ratioBounds : List (ℚ × ℚ) → Option (ℚ × ℚ)
  | [] => none
  | p :: ps => if 0 < p.1 then includeRatio (p.2 / p.1) (ratioBounds ps) else ratioBounds ps

/-- Input gate: `Φ_A ≥ 0`, `Φ_B ≥ 0` and the support condition `Φ_A = 0 → Φ_B = 0`. -/
def separatorWF (phi : List (ℚ × ℚ)) : Bool :=
  phi.all fun p => decide (0 ≤ p.1 ∧ 0 ≤ p.2 ∧ (p.1 = 0 → p.2 = 0))

/-- **Executable separator checker.** Given the list of `(Φ_A(τ), Φ_B(τ))` over boundary
states, returns `[min ρ, max ρ]`, or `none` on malformed input, on a support violation
(the sharp upper bound would be `+∞`), or when `T⁺ = ∅`. -/
def checkSeparator (phi : List (ℚ × ℚ)) : Option (ℚ × ℚ) :=
  if separatorWF phi then ratioBounds phi else none

theorem ratioBounds_none : ∀ (phi : List (ℚ × ℚ)), ratioBounds phi = none →
    ∀ p ∈ phi, ¬ 0 < p.1
  | [], _ => by simp
  | p :: ps, h => by
    simp only [ratioBounds] at h
    split_ifs at h with hp
    · cases hr : ratioBounds ps <;> simp [hr, includeRatio] at h
    · intro q hq
      rcases List.mem_cons.mp hq with rfl | hq
      · exact hp
      · exact ratioBounds_none ps h q hq

/-- Specification of `ratioBounds`: the endpoints bound every ratio on the support and are
attained. -/
theorem ratioBounds_spec : ∀ (phi : List (ℚ × ℚ)) {lo hi : ℚ},
    ratioBounds phi = some (lo, hi) →
      (∀ p ∈ phi, 0 < p.1 → lo ≤ p.2 / p.1 ∧ p.2 / p.1 ≤ hi) ∧
      (∃ p ∈ phi, 0 < p.1 ∧ p.2 / p.1 = lo) ∧ (∃ p ∈ phi, 0 < p.1 ∧ p.2 / p.1 = hi)
  | [], _, _, h => by simp [ratioBounds] at h
  | p :: ps, lo, hi, h => by
    simp only [ratioBounds] at h
    split_ifs at h with hp
    · rcases hr : ratioBounds ps with _ | ⟨lo', hi'⟩
      · rw [hr] at h
        simp only [includeRatio, Option.some.injEq, Prod.mk.injEq] at h
        obtain ⟨rfl, rfl⟩ := h
        refine ⟨fun q hq hq1 => ?_, ⟨p, by simp, hp, rfl⟩, ⟨p, by simp, hp, rfl⟩⟩
        rcases List.mem_cons.mp hq with rfl | hq
        · exact ⟨le_rfl, le_rfl⟩
        · exact absurd hq1 (ratioBounds_none ps hr q hq)
      · rw [hr] at h
        simp only [includeRatio, Option.some.injEq, Prod.mk.injEq] at h
        obtain ⟨rfl, rfl⟩ := h
        obtain ⟨hb, ⟨q₁, hq₁, hq₁p, hq₁e⟩, ⟨q₂, hq₂, hq₂p, hq₂e⟩⟩ := ratioBounds_spec ps hr
        refine ⟨fun q hq hq1 => ?_, ?_, ?_⟩
        · rcases List.mem_cons.mp hq with rfl | hq
          · exact ⟨min_le_left _ _, le_max_left _ _⟩
          · exact ⟨(min_le_right _ _).trans (hb q hq hq1).1,
              (hb q hq hq1).2.trans (le_max_right _ _)⟩
        · rcases min_choice (p.2 / p.1) lo' with h' | h'
          · exact ⟨p, by simp, hp, h'.symm⟩
          · exact ⟨q₁, List.mem_cons_of_mem _ hq₁, hq₁p, hq₁e.trans h'.symm⟩
        · rcases max_choice (p.2 / p.1) hi' with h' | h'
          · exact ⟨p, by simp, hp, h'.symm⟩
          · exact ⟨q₂, List.mem_cons_of_mem _ hq₂, hq₂p, hq₂e.trans h'.symm⟩
    · obtain ⟨hb, h₁, h₂⟩ := ratioBounds_spec ps h
      refine ⟨fun q hq hq1 => ?_, ?_, ?_⟩
      · rcases List.mem_cons.mp hq with rfl | hq
        · exact absurd hq1 hp
        · exact hb q hq hq1
      · obtain ⟨q, hq, h⟩ := h₁; exact ⟨q, List.mem_cons_of_mem _ hq, h⟩
      · obtain ⟨q, hq, h⟩ := h₂; exact ⟨q, List.mem_cons_of_mem _ hq, h⟩

variable {phi : List (ℚ × ℚ)}

theorem checkSeparator_wf {bounds : ℚ × ℚ} (h : checkSeparator phi = some bounds) :
    separatorWF phi = true ∧ ratioBounds phi = some bounds := by
  unfold checkSeparator at h
  split_ifs at h with hwf
  exact ⟨hwf, h⟩

/-- **Soundness of the separator checker.** If `checkSeparator` returns `(lo, hi)`, then for
every exterior marginal `Ψ ≥ 0` over the listed boundary states with `Z_A > 0`,
`lo ≤ Z_B / Z_A ≤ hi`. -/
theorem checkSeparator_sound {lo hi : ℚ} (h : checkSeparator phi = some (lo, hi))
    (ψ : Fin phi.length → ℝ) (hψ : ∀ τ, 0 ≤ ψ τ)
    (hZ : 0 < ∑ τ, (phi[τ].1 : ℝ) * ψ τ) :
    (∑ τ, (phi[τ].2 : ℝ) * ψ τ) / (∑ τ, (phi[τ].1 : ℝ) * ψ τ) ∈
      Set.Icc (lo : ℝ) (hi : ℝ) := by
  obtain ⟨hwf, hb⟩ := checkSeparator_wf h
  obtain ⟨hbounds, -, -⟩ := ratioBounds_spec phi hb
  have hent : ∀ τ : Fin phi.length,
      0 ≤ phi[τ].1 ∧ 0 ≤ phi[τ].2 ∧ (phi[τ].1 = 0 → phi[τ].2 = 0) := fun τ => by
    have := List.all_eq_true.mp hwf phi[τ] (List.getElem_mem τ.isLt)
    simpa only [decide_eq_true_eq] using this
  refine ratio_mem_Icc_of_local (fun τ => by exact_mod_cast (hent τ).1) hψ
    (fun τ h0 => by exact_mod_cast (hent τ).2.2 (by exact_mod_cast h0)) hZ
    (fun τ hτ => ?_) (fun τ hτ => ?_)
  · have := (hbounds phi[τ] (List.getElem_mem τ.isLt) (by exact_mod_cast hτ)).1
    exact_mod_cast this
  · have := (hbounds phi[τ] (List.getElem_mem τ.isLt) (by exact_mod_cast hτ)).2
    exact_mod_cast this

/-- **End-to-end soundness for factorised systems.** For any interiors whose marginals are the
listed values and any nonnegative exterior with `Z_A > 0`, `Z_B/Z_A ∈ [lo, hi]`. -/
theorem checkSeparator_sound_sys {lo hi : ℚ} (h : checkSeparator phi = some (lo, hi))
    {IA IB : Type*} [Fintype IA] [Fintype IB] {a : IA → Fin phi.length → ℝ}
    {b : IB → Fin phi.length → ℝ} (hA : ∀ τ, Phi a τ = (phi[τ].1 : ℝ))
    (hB : ∀ τ, Phi b τ = (phi[τ].2 : ℝ)) {e : Fin phi.length → O → ℝ}
    (he : ∀ τ o, 0 ≤ e τ o) (hZ : 0 < Zsys a e) :
    Zsys b e / Zsys a e ∈ Set.Icc (lo : ℝ) (hi : ℝ) := by
  rw [Zsys_eq_sum] at hZ ⊢
  rw [Zsys_eq_sum]
  simp only [hA, hB] at hZ ⊢
  exact checkSeparator_sound h _ (Psi_nonneg he) hZ

/-- **Exactness of the separator checker.** Both returned endpoints are local ratios of
boundary states in the support, so (Theorem 2.2) they are attained by exteriors. -/
theorem checkSeparator_attained {lo hi : ℚ} (h : checkSeparator phi = some (lo, hi)) :
    (∃ p ∈ phi, 0 < p.1 ∧ p.2 / p.1 = lo) ∧ (∃ p ∈ phi, 0 < p.1 ∧ p.2 / p.1 = hi) :=
  (ratioBounds_spec phi (checkSeparator_wf h).2).2

end Executable

end CoupledSeparator
end FoldBounds
