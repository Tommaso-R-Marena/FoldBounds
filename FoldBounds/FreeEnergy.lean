module

public import FoldBounds.Probability

/-!
# FoldBounds.FreeEnergy — certified Helmholtz free-energy bounds

`F = -(1/β) log Z` (energies in units where `k_B T = 1/β`).  The temperature assumption
is `0 < β` (positive, finite temperature); positivity of `Z` is automatic for nonempty `C`
(`Z_pos`) but the *lower end* `Zlo` of a certified enclosure must be shown positive
separately — `certifiedZLower_pos` does this for block certificates.

Note: with Lean's convention `x / 0 = 0`, `freeEnergy 0 E = 0`; all theorems below assume
`0 < β`, so this convention is never used.
-/

@[expose] public section

open Finset Real

namespace FoldBounds

variable {C : Type*} [Fintype C]

/-- Helmholtz free energy `F = -(1/β) log Z`. -/
noncomputable def freeEnergy (β : ℝ) (E : C → ℝ) : ℝ := -Real.log (Z β E) / β

/-- Free energy of a macrostate `A`: `F_A = -(1/β) log Z_A`. -/
noncomputable def freeEnergyOn (β : ℝ) (E : C → ℝ) (A : Finset C) : ℝ :=
  -Real.log (ZOn β E A) / β

variable {β : ℝ} (E : C → ℝ)

/-- Monotone transfer of an enclosure of a positive quantity through `x ↦ -(log x)/β`. -/
theorem neg_log_div_mem_Icc (hβ : 0 < β) {x lo hi : ℝ} (hlo : 0 < lo)
    (hx : x ∈ Set.Icc lo hi) :
    -Real.log x / β ∈ Set.Icc (-Real.log hi / β) (-Real.log lo / β) := by
  have hxpos : 0 < x := hlo.trans_le hx.1
  have h1 := Real.log_le_log hlo hx.1
  have h2 := Real.log_le_log hxpos hx.2
  exact ⟨div_le_div_of_nonneg_right (by linarith) hβ.le,
    div_le_div_of_nonneg_right (by linarith) hβ.le⟩

/-- **Theorem 6.** A certified enclosure `0 < Zlo ≤ Z ≤ Zhi` gives the certified
enclosure `-(log Zhi)/β ≤ F ≤ -(log Zlo)/β` for `0 < β`. -/
theorem freeEnergy_mem_Icc (hβ : 0 < β) {Zlo Zhi : ℝ} (hZlo : 0 < Zlo)
    (hZ : Z β E ∈ Set.Icc Zlo Zhi) :
    freeEnergy β E ∈ Set.Icc (-Real.log Zhi / β) (-Real.log Zlo / β) :=
  neg_log_div_mem_Icc hβ hZlo hZ

/-- The free energy never exceeds the energy of any conformation (`F ≤ E_min`). -/
theorem freeEnergy_le_energy [Nonempty C] (hβ : 0 < β) (c : C) :
    freeEnergy β E ≤ E c := by
  have hw : weight β E c ≤ Z β E := by
    simpa [ZOn] using ZOn_le_Z β E {c}
  have hlog := Real.log_le_log (weight_pos β E c) hw
  rw [weight, Real.log_exp] at hlog
  rw [freeEnergy, div_le_iff₀ hβ]
  linarith

/-- Lower bound by the ground-state energy and the entropy of the full space:
`Emin - (log |C|)/β ≤ F` whenever `Emin ≤ E c` for all `c`. -/
theorem le_freeEnergy_of_energy_ge [Nonempty C] (hβ : 0 < β) (Emin : ℝ)
    (h : ∀ c, Emin ≤ E c) :
    Emin - Real.log (Fintype.card C) / β ≤ freeEnergy β E := by
  have hZ : Z β E ≤ (Fintype.card C : ℝ) * Real.exp (-β * Emin) := by
    simpa using ZOn_le_card_mul_exp E hβ.le Finset.univ Emin (fun c _ => h c)
  have hcard : (0 : ℝ) < Fintype.card C := by exact_mod_cast Fintype.card_pos
  have hlog := Real.log_le_log (Z_pos β E) hZ
  rw [Real.log_mul hcard.ne' (Real.exp_pos _).ne', Real.log_exp] at hlog
  have e : Emin - Real.log (Fintype.card C) / β =
      (β * Emin - Real.log (Fintype.card C)) / β := by
    field_simp
  rw [freeEnergy, e]
  exact div_le_div_of_nonneg_right (by linarith) hβ.le

/-- Variational (enumeration) bound: for any nonempty explored set `S`, `F ≤ F_S`. -/
theorem freeEnergy_le_freeEnergyOn (hβ : 0 < β) {S : Finset C} (hS : S.Nonempty) :
    freeEnergy β E ≤ freeEnergyOn β E S := by
  have h := Real.log_le_log (ZOn_pos β E hS) (ZOn_le_Z β E S)
  exact div_le_div_of_nonneg_right (by linarith) hβ.le

/-- Free-energy representation of probabilities: `P(A) = exp (-β (F_A - F))`. -/
theorem prob_eq_exp_freeEnergy [Nonempty C] (hβ : 0 < β) {A : Finset C}
    (hA : A.Nonempty) :
    prob β E A = Real.exp (-β * (freeEnergyOn β E A - freeEnergy β E)) := by
  have hb : β ≠ 0 := hβ.ne'
  have e1 : -β * (freeEnergyOn β E A - freeEnergy β E) =
      Real.log (ZOn β E A) - Real.log (Z β E) := by
    simp only [freeEnergyOn, freeEnergy]; field_simp; ring
  rw [e1, Real.exp_sub, Real.exp_log (ZOn_pos β E hA), Real.exp_log (Z_pos β E), prob]

/-- **Certified free-energy differences.** A certified probability enclosure
`0 < pL ≤ P(A) ≤ pU` yields `-(log pU)/β ≤ F_A - F ≤ -(log pL)/β`. -/
theorem freeEnergy_diff_mem_Icc [Nonempty C] (hβ : 0 < β) {A : Finset C}
    (hA : A.Nonempty) {pL pU : ℝ} (hpL : 0 < pL) (hp : prob β E A ∈ Set.Icc pL pU) :
    freeEnergyOn β E A - freeEnergy β E ∈
      Set.Icc (-Real.log pU / β) (-Real.log pL / β) := by
  have key : freeEnergyOn β E A - freeEnergy β E = -Real.log (prob β E A) / β := by
    rw [prob_eq_exp_freeEnergy E hβ hA, Real.log_exp]
    field_simp
  rw [key]
  exact neg_log_div_mem_Icc hβ hpL hp

/-- The certified lower end of the block enclosure of `Z` is strictly positive. -/
theorem certifiedZLower_pos [DecidableEq C] [Nonempty C] {ι : Type*} (blk : C → ι)
    (u : ι → ℝ) (S : Finset C) :
    0 < ZOn β E S + blockBound β blk u Sᶜ := by
  rcases S.eq_empty_or_nonempty with rfl | hS
  · simp only [ZOn_empty, Finset.compl_empty, zero_add, blockBound]
    exact Finset.sum_pos (fun _ _ => Real.exp_pos _) Finset.univ_nonempty
  · exact add_pos_of_pos_of_nonneg (ZOn_pos β E hS) (blockBound_nonneg β blk u _)

/-- **Theorem 6 (from partial enumeration).** Exact enumeration of `S` and block energy
certificates on `Sᶜ` give a certified free-energy enclosure, for `0 < β`. -/
theorem freeEnergy_mem_Icc_of_blocks [DecidableEq C] [Nonempty C] {ι : Type*}
    (hβ : 0 < β) (blk : C → ι) (ℓ u : ι → ℝ) (S : Finset C)
    (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    freeEnergy β E ∈
      Set.Icc (-Real.log (ZOn β E S + blockBound β blk ℓ Sᶜ) / β)
        (-Real.log (ZOn β E S + blockBound β blk u Sᶜ) / β) :=
  freeEnergy_mem_Icc E hβ (certifiedZLower_pos E blk u S)
    (Z_mem_Icc_of_blocks E blk hβ.le S ℓ u hℓ hu)

end FoldBounds
