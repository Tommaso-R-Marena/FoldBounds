module

public import FoldBounds.CoupledSeparator

/-!
# FoldBounds.CoupledBinding — certified affinity and selectivity (CEC §3, Theorems 4.1–4.3)

Pocket states `τ ∈ T` (finite).  `ΨR τ ≥ 0` is the apo receptor weight with the pocket in
state `τ`, summed over everything else; `ΦL τ` is the ligand's summed weight in the pocket in
state `τ`; `ZL` is the ligand's free weight.  Under the locality assumption (A-loc) the complex
partition function is `Z_RL = Σ_τ ΨR(τ) ΦL(τ)`, so the binding constant is the ratio of finite
sums `K_L = Z_RL / (Z_R Z_L)` (`Kbind`).  `Kbind_of_sys` derives this form from a microscopic
factorised model through Lemma 2.0.

* `Kbind_eq_sum`: the conformational-selection identity `K_L = Σ_τ p_apo(τ) k_L(τ)`;
* Theorem 4.1 `affinity_mem_Icc`, `affinity_sharp`;
* Theorem 4.2 `selectivity_mem_Icc`, `selectivity_sharp`, `L2_tighter_of_one_lt_min`
  (requires `k₁(τ) > 0` for every `τ`, which the spec leaves implicit;
  `selectivity_needs_positive_k1` is a counterexample without it);
* Theorem 4.3 `nonidentifiable` (two explicit point-mass receptor ensembles),
  `nonidentifiable_of_min_lt_one_lt_max`, `nonidentifiable_pos` (strictly positive
  ensembles).
-/

@[expose] public section

open Finset

namespace FoldBounds
namespace CoupledBinding

open CoupledSeparator

variable {T : Type*} [Fintype T]

/-- Binding constant `K_L = Z_RL / (Z_R Z_L)` with `Z_R = Σ_τ ΨR(τ)` and
`Z_RL = Σ_τ ΨR(τ) ΦL(τ)`. -/
noncomputable def Kbind (ΨR ΦL : T → ℝ) (ZL : ℝ) : ℝ :=
  (∑ τ, ΨR τ * ΦL τ) / ((∑ τ, ΨR τ) * ZL)

/-- Local affinity `k_L(τ) = Φ_L(τ) / Z_L`. -/
noncomputable def kloc (ΦL : T → ℝ) (ZL : ℝ) (τ : T) : ℝ := ΦL τ / ZL

/-- Apo pocket population `p_apo(τ) = Ψ_R(τ) / Z_R`. -/
noncomputable def papo (ΨR : T → ℝ) (τ : T) : ℝ := ΨR τ / ∑ τ, ΨR τ

/-- **Conformational selection** (classical): `K_L = Σ_τ p_apo(τ) k_L(τ)`. -/
theorem Kbind_eq_sum (ΨR ΦL : T → ℝ) (ZL : ℝ) :
    Kbind ΨR ΦL ZL = ∑ τ, papo ΨR τ * kloc ΦL ZL τ := by
  rw [Kbind, sum_div]
  exact sum_congr rfl fun τ _ => by rw [papo, kloc, div_mul_div_comm]

/-- The microscopic model: receptor exterior `eR : T → O → ℝ` (apo weights) and bound-ligand
pose weights `wL : P → T → ℝ` touching only the pocket (A-loc).  Then
`Z_RL / (Z_R Z_L) = Kbind (Ψ eR) (Φ wL) Z_L` by Lemma 2.0. -/
theorem Kbind_of_sys {O P : Type*} [Fintype O] [Fintype P] (eR : T → O → ℝ) (wL : P → T → ℝ)
    (ZL : ℝ) :
    Zsys wL eR / (Zsys (fun (_ : Unit) (_ : T) => (1 : ℝ)) eR * ZL) =
      Kbind (Psi eR) (Phi wL) ZL := by
  rw [Zsys_eq_sum, Zsys_eq_sum, Kbind]
  simp only [Phi, Finset.univ_unique, sum_singleton, one_mul]
  congr 1
  exact sum_congr rfl fun τ _ => mul_comm _ _

theorem Kbind_eq_ratio (ΨR ΦL : T → ℝ) (ZL : ℝ) :
    Kbind ΨR ΦL ZL = (∑ τ, kloc ΦL ZL τ * ΨR τ) / ∑ τ, (fun _ => (1 : ℝ)) τ * ΨR τ := by
  rw [Kbind_eq_sum, sum_div]
  exact sum_congr rfl fun τ _ => by simp only [papo, one_mul]; ring

/-! ### Theorem 4.1: certified affinity -/

/-- **Theorem 4.1 (enclosure).** For every receptor ensemble (`Ψ_R ≥ 0`, `Z_R > 0`),
`min_τ k_L(τ) ≤ K_L ≤ max_τ k_L(τ)`. -/
theorem affinity_mem_Icc [Nonempty T] {ΨR : T → ℝ} (hΨ : ∀ τ, 0 ≤ ΨR τ)
    (hZR : 0 < ∑ τ, ΨR τ) (ΦL : T → ℝ) (ZL : ℝ) :
    Kbind ΨR ΦL ZL ∈
      Set.Icc (univ.inf' univ_nonempty (kloc ΦL ZL)) (univ.sup' univ_nonempty (kloc ΦL ZL)) := by
  rw [Kbind_eq_ratio]
  exact ratio_mem_Icc_of_local (fun _ => zero_le_one) hΨ (fun _ h => absurd h one_ne_zero)
    (by simpa using hZR) (fun τ _ => by simpa using inf'_le (kloc ΦL ZL) (mem_univ τ))
    (fun τ _ => by simpa using le_sup' (kloc ΦL ZL) (mem_univ τ))

/-- **Theorem 4.1 (sharpness).** Every value of `[min k_L, max k_L]` is the binding constant of
some receptor ensemble. -/
theorem affinity_sharp [Nonempty T] (ΦL : T → ℝ) (ZL : ℝ) {κ : ℝ}
    (hκ : κ ∈ Set.Icc (univ.inf' univ_nonempty (kloc ΦL ZL))
      (univ.sup' univ_nonempty (kloc ΦL ZL))) :
    ∃ ΨR : T → ℝ, (∀ τ, 0 ≤ ΨR τ) ∧ 0 < ∑ τ, ΨR τ ∧ Kbind ΨR ΦL ZL = κ := by
  obtain ⟨τ₁, -, he₁⟩ := exists_mem_eq_inf' univ_nonempty (kloc ΦL ZL)
  obtain ⟨τ₂, -, he₂⟩ := exists_mem_eq_sup' univ_nonempty (kloc ΦL ZL)
  have hl : kloc ΦL ZL τ₁ ≤ κ := by rw [← he₁]; exact hκ.1
  have hu : κ ≤ kloc ΦL ZL τ₂ := by rw [← he₂]; exact hκ.2
  have hseg : κ ∈ segment ℝ (kloc ΦL ZL τ₁ / (fun _ => (1 : ℝ)) τ₁)
      (kloc ΦL ZL τ₂ / (fun _ => (1 : ℝ)) τ₂) := by
    simp only [div_one]
    rw [segment_eq_Icc (hl.trans hu)]
    exact ⟨hl, hu⟩
  obtain ⟨ψ, hψ, -, hA, hB⟩ := exists_two_point (φA := fun _ => (1 : ℝ)) one_pos one_pos hseg
  have hA' : ∑ τ, ψ τ = 1 := by simpa using hA
  refine ⟨ψ, hψ, by rw [hA']; exact one_pos, ?_⟩
  rw [Kbind_eq_ratio, hA]
  simpa [mul_comm] using hB

/-! ### Theorem 4.2: certified selectivity -/

theorem kloc_weighted_pos [Nonempty T] {ΨR k : T → ℝ} (hΨ : ∀ τ, 0 ≤ ΨR τ)
    (hZR : 0 < ∑ τ, ΨR τ) (hk : ∀ τ, 0 < k τ) : 0 < ∑ τ, k τ * ΨR τ := by
  obtain ⟨τ, -, hτ⟩ := CoupledRatio.exists_pos_of_sum_pos hZR
  exact sum_pos' (fun τ _ => mul_nonneg (hk τ).le (hΨ τ)) ⟨τ, mem_univ _, mul_pos (hk τ) hτ⟩

/-- The selectivity `S = K_{L2} / K_{L1}` is the exterior-weighted ratio of the local
affinities; the receptor partition function cancels. -/
theorem selectivity_eq_ratio {ΨR : T → ℝ} (hZR : ∑ τ, ΨR τ ≠ 0) (Φ1 Φ2 : T → ℝ) (Z1 Z2 : ℝ) :
    Kbind ΨR Φ2 Z2 / Kbind ΨR Φ1 Z1 =
      (∑ τ, kloc Φ2 Z2 τ * ΨR τ) / ∑ τ, kloc Φ1 Z1 τ * ΨR τ := by
  rw [Kbind_eq_ratio, Kbind_eq_ratio]
  exact div_div_div_cancel_right₀ (by simpa using hZR) _ _

variable {Φ1 Φ2 : T → ℝ} {Z1 Z2 : ℝ}

/-- **Theorem 4.2 (enclosure).** If `k₁(τ) > 0` for every pocket state, then for every receptor
ensemble `min_τ k₂/k₁ ≤ K_{L2}/K_{L1} ≤ max_τ k₂/k₁`. -/
theorem selectivity_mem_Icc [Nonempty T] {ΨR : T → ℝ} (hΨ : ∀ τ, 0 ≤ ΨR τ)
    (hZR : 0 < ∑ τ, ΨR τ) (hk1 : ∀ τ, 0 < kloc Φ1 Z1 τ) :
    Kbind ΨR Φ2 Z2 / Kbind ΨR Φ1 Z1 ∈
      Set.Icc (univ.inf' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ)
        (univ.sup' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ) := by
  rw [selectivity_eq_ratio hZR.ne']
  exact ratio_mem_Icc_of_local (fun τ => (hk1 τ).le) hΨ (fun τ h => absurd h (hk1 τ).ne')
    (kloc_weighted_pos hΨ hZR hk1)
    (fun τ _ => inf'_le (fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ) (mem_univ τ))
    (fun τ _ => le_sup' (fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ) (mem_univ τ))

/-- **Theorem 4.2 (sharpness).** Every value of `[min k₂/k₁, max k₂/k₁]` is the selectivity of
some receptor ensemble. -/
theorem selectivity_sharp [Nonempty T] (hk1 : ∀ τ, 0 < kloc Φ1 Z1 τ) {σ : ℝ}
    (hσ : σ ∈ Set.Icc (univ.inf' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ)
      (univ.sup' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ)) :
    ∃ ΨR : T → ℝ, (∀ τ, 0 ≤ ΨR τ) ∧ 0 < ∑ τ, ΨR τ ∧ Kbind ΨR Φ2 Z2 / Kbind ΨR Φ1 Z1 = σ := by
  obtain ⟨τ₁, -, he₁⟩ := exists_mem_eq_inf' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ
  obtain ⟨τ₂, -, he₂⟩ := exists_mem_eq_sup' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ
  have hl : kloc Φ2 Z2 τ₁ / kloc Φ1 Z1 τ₁ ≤ σ := by rw [← he₁]; exact hσ.1
  have hu : σ ≤ kloc Φ2 Z2 τ₂ / kloc Φ1 Z1 τ₂ := by rw [← he₂]; exact hσ.2
  have hseg : σ ∈ segment ℝ (kloc Φ2 Z2 τ₁ / kloc Φ1 Z1 τ₁) (kloc Φ2 Z2 τ₂ / kloc Φ1 Z1 τ₂) := by
    rw [segment_eq_Icc (hl.trans hu)]
    exact ⟨hl, hu⟩
  obtain ⟨ψ, hψ, -, hA, hB⟩ := exists_two_point (hk1 τ₁) (hk1 τ₂) hseg
  have hpos : 0 < ∑ τ, ψ τ := by
    by_contra h0
    have hz := (sum_eq_zero_iff_of_nonneg fun τ _ => hψ τ).mp
      (le_antisymm (not_lt.mp h0) (sum_nonneg fun τ _ => hψ τ))
    rw [sum_eq_zero fun τ hτ => by rw [hz τ hτ, mul_zero]] at hA
    exact zero_ne_one hA
  refine ⟨ψ, hψ, hpos, ?_⟩
  rw [selectivity_eq_ratio hpos.ne', hA, hB, div_one]

/-- **Theorem 4.2 (certified ordering).** If `min_τ k₂/k₁ > 1`, L2 binds more tightly than L1
for every receptor ensemble. -/
theorem L2_tighter_of_one_lt_min [Nonempty T] {ΨR : T → ℝ} (hΨ : ∀ τ, 0 ≤ ΨR τ)
    (hZR : 0 < ∑ τ, ΨR τ) (hk1 : ∀ τ, 0 < kloc Φ1 Z1 τ)
    (h1 : 1 < univ.inf' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ) :
    Kbind ΨR Φ1 Z1 < Kbind ΨR Φ2 Z2 := by
  have hK1 : 0 < Kbind ΨR Φ1 Z1 := by
    rw [Kbind_eq_ratio]
    exact div_pos (kloc_weighted_pos hΨ hZR hk1) (by simpa using hZR)
  exact (one_lt_div hK1).mp (h1.trans_le (selectivity_mem_Icc hΨ hZR hk1).1)

/-- **Spec discrepancy for Theorem 4.2.** Without `k₁ > 0` everywhere the enclosure fails:
`k₁ = (1, 0)`, `k₂ = (1, 1)` and the uniform receptor ensemble give `S = 2`, while
`k₂/k₁ ≤ 1` at every pocket state (Lean's `x / 0 = 0`; the ratio is `+∞` mathematically). -/
theorem selectivity_needs_positive_k1 :
    ∃ (ΨR Φ1 Φ2 : Bool → ℝ), (∀ τ, 0 ≤ ΨR τ) ∧ 0 < ∑ τ, ΨR τ ∧ (∀ τ, 0 ≤ Φ1 τ) ∧
      (∀ τ, kloc Φ2 1 τ / kloc Φ1 1 τ ≤ 1) ∧ Kbind ΨR Φ2 1 / Kbind ΨR Φ1 1 = 2 := by
  refine ⟨fun _ => 1, fun x => if x then 1 else 0, fun _ => 1, fun _ => zero_le_one, ?_,
    fun τ => by cases τ <;> simp, fun τ => by cases τ <;> simp [kloc], ?_⟩
  · simp
  · norm_num [Kbind]

/-! ### Theorem 4.3: non-identifiability -/

/-- Point-mass receptor ensemble concentrated on the pocket state `τ₀`. -/
def pointEnsemble [DecidableEq T] (τ₀ : T) : T → ℝ := fun τ => if τ = τ₀ then 1 else 0

theorem Kbind_point [DecidableEq T] (τ₀ : T) (ΦL : T → ℝ) (ZL : ℝ) :
    Kbind (pointEnsemble τ₀) ΦL ZL = kloc ΦL ZL τ₀ := by
  simp [Kbind, pointEnsemble, kloc]

omit [Fintype T] in
theorem pointEnsemble_nonneg [DecidableEq T] (τ₀ τ : T) : 0 ≤ pointEnsemble τ₀ τ := by
  unfold pointEnsemble; split_ifs <;> norm_num

theorem pointEnsemble_sum [DecidableEq T] (τ₀ : T) : ∑ τ, pointEnsemble τ₀ τ = 1 := by
  simp [pointEnsemble]

/-- **Theorem 4.3 (non-identifiability), explicit ensembles.** If `k₂/k₁ < 1` at `τm` and
`k₂/k₁ > 1` at `τp` (with `k₁ > 0` there), the receptor ensemble concentrated on `τp` makes L2
the tighter binder and the one concentrated on `τm` makes L1 the tighter binder. -/
theorem nonidentifiable [DecidableEq T] {τm τp : T} (hm1 : 0 < kloc Φ1 Z1 τm)
    (hp1 : 0 < kloc Φ1 Z1 τp) (hm : kloc Φ2 Z2 τm / kloc Φ1 Z1 τm < 1)
    (hp : 1 < kloc Φ2 Z2 τp / kloc Φ1 Z1 τp) :
    Kbind (pointEnsemble τp) Φ1 Z1 < Kbind (pointEnsemble τp) Φ2 Z2 ∧
      Kbind (pointEnsemble τm) Φ2 Z2 < Kbind (pointEnsemble τm) Φ1 Z1 := by
  rw [Kbind_point, Kbind_point, Kbind_point, Kbind_point]
  exact ⟨(one_lt_div hp1).mp hp, (div_lt_one hm1).mp hm⟩

/-- **Theorem 4.3 (min/max form).** If `min_τ k₂/k₁ < 1 < max_τ k₂/k₁` (and `k₁ > 0`), some
receptor ensemble makes L2 the tighter binder and another makes L1 the tighter binder. -/
theorem nonidentifiable_of_min_lt_one_lt_max [Nonempty T] [DecidableEq T]
    (hk1 : ∀ τ, 0 < kloc Φ1 Z1 τ)
    (hlo : (univ.inf' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ) < 1)
    (hhi : 1 < univ.sup' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ) :
    ∃ τm τp : T,
      Kbind (pointEnsemble τp) Φ1 Z1 < Kbind (pointEnsemble τp) Φ2 Z2 ∧
      Kbind (pointEnsemble τm) Φ2 Z2 < Kbind (pointEnsemble τm) Φ1 Z1 := by
  obtain ⟨τm, -, hem⟩ := exists_mem_eq_inf' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ
  obtain ⟨τp, -, hep⟩ := exists_mem_eq_sup' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ
  exact ⟨τm, τp, nonidentifiable (hk1 τm) (hk1 τp) (hem ▸ hlo) (hep ▸ hhi)⟩

/-- **Theorem 4.3 with strictly positive ensembles.** Under the same hypotheses there are
receptor ensembles with `Ψ_R > 0` everywhere on both sides of the ordering. -/
theorem nonidentifiable_pos [Nonempty T] (hk1 : ∀ τ, 0 < kloc Φ1 Z1 τ)
    (hlo : (univ.inf' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ) < 1)
    (hhi : 1 < univ.sup' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ) :
    ∃ Ψ₁ Ψ₂ : T → ℝ, (∀ τ, 0 < Ψ₁ τ) ∧ (∀ τ, 0 < Ψ₂ τ) ∧
      Kbind Ψ₁ Φ1 Z1 < Kbind Ψ₁ Φ2 Z2 ∧ Kbind Ψ₂ Φ2 Z2 < Kbind Ψ₂ Φ1 Z1 := by
  obtain ⟨τm, -, hem⟩ := exists_mem_eq_inf' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ
  obtain ⟨τp, -, hep⟩ := exists_mem_eq_sup' univ_nonempty fun τ => kloc Φ2 Z2 τ / kloc Φ1 Z1 τ
  have hlo' : kloc Φ2 Z2 τm / kloc Φ1 Z1 τm < 1 := hem ▸ hlo
  have hhi' : 1 < kloc Φ2 Z2 τp / kloc Φ1 Z1 τp := hep ▸ hhi
  have hk1' : ∀ τ, 0 ≤ kloc Φ1 Z1 τ := fun τ => (hk1 τ).le
  obtain ⟨Ψ₁, hΨ₁, hZ₁, hq₁⟩ := exists_pos_weight (φB := kloc Φ2 Z2) hk1' (hk1 τm) (hk1 τp)
    (lam := (1 + kloc Φ2 Z2 τp / kloc Φ1 Z1 τp) / 2) (by linarith) (by linarith)
  obtain ⟨Ψ₂, hΨ₂, hZ₂, hq₂⟩ := exists_pos_weight (φB := kloc Φ2 Z2) hk1' (hk1 τm) (hk1 τp)
    (lam := (kloc Φ2 Z2 τm / kloc Φ1 Z1 τm + 1) / 2) (by linarith) (by linarith)
  have hpos : ∀ Ψ : T → ℝ, (∀ τ, 0 < Ψ τ) → 0 < ∑ τ, Ψ τ := fun Ψ h =>
    sum_pos (fun τ _ => h τ) univ_nonempty
  have hK : ∀ Ψ : T → ℝ, (∀ τ, 0 < Ψ τ) → 0 < Kbind Ψ Φ1 Z1 := fun Ψ h => by
    rw [Kbind_eq_ratio]
    exact div_pos (kloc_weighted_pos (fun τ => (h τ).le) (hpos Ψ h) hk1) (by simpa using hpos Ψ h)
  have hS : ∀ Ψ : T → ℝ, (∀ τ, 0 < Ψ τ) →
      Kbind Ψ Φ2 Z2 / Kbind Ψ Φ1 Z1 =
        (∑ τ, kloc Φ2 Z2 τ * Ψ τ) / ∑ τ, kloc Φ1 Z1 τ * Ψ τ := fun Ψ h =>
    selectivity_eq_ratio (hpos Ψ h).ne' Φ1 Φ2 Z1 Z2
  refine ⟨Ψ₁, Ψ₂, hΨ₁, hΨ₂, ?_, ?_⟩
  · have h := hS Ψ₁ hΨ₁
    rw [hq₁] at h
    exact (one_lt_div (hK Ψ₁ hΨ₁)).mp (by rw [h]; linarith)
  · have h := hS Ψ₂ hΨ₂
    rw [hq₂] at h
    exact (div_lt_one (hK Ψ₂ hΨ₂)).mp (by rw [h]; linarith)

end CoupledBinding
end FoldBounds
