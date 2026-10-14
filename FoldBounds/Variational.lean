module

public import FoldBounds.Compositional

/-!
# FoldBounds.Variational — Gibbs–Bogoliubov certificates

For a split `E = E₀ + V` into a reference landscape `E₀` (whose ensemble is tractable,
e.g. decoupled domains) and a perturbation `V`, the Gibbs–Bogoliubov inequality gives

`Z(E₀ + V) ≥ Z(E₀) · exp (-β ⟨V⟩₀)`   (every real `β`), and
`F(E₀ + V) ≤ F(E₀) + ⟨V⟩₀`             (`0 < β`),

where `⟨V⟩₀ = ∑_c w₀(c) V(c) / Z(E₀)` is the reference-ensemble average.  Since
`⟨V⟩₀ ≤ max V`, this sharpens the one-sided compositional bound
`Z ≥ Z₁ Z₂ e^{-β Jhi}` of `FoldBounds.Compositional` (see `Z_coupledEnergy_ge_gibbs`).
The proof uses only the tangent-line inequality `exp x ≥ exp m · (1 + x - m)`.
-/

@[expose] public section

open Finset Real

namespace FoldBounds

variable {C : Type*} [Fintype C]

/-- Reference-ensemble average `⟨V⟩₀ = ∑ w₀ V / Z₀`. -/
noncomputable def refAverage (β : ℝ) (E₀ V : C → ℝ) : ℝ :=
  (∑ c, weight β E₀ c * V c) / Z β E₀

/-- **Gibbs–Bogoliubov lower bound on `Z`** (any real `β`). -/
theorem Z_ge_gibbs_bogoliubov [Nonempty C] (β : ℝ) (E₀ V : C → ℝ) :
    Z β E₀ * Real.exp (-β * refAverage β E₀ V) ≤ Z β (fun c => E₀ c + V c) := by
  set m := -β * refAverage β E₀ V
  have hZ := Z_pos β E₀
  have hsum : ∑ c, weight β E₀ c * (-β * V c) = m * Z β E₀ := by
    simp only [m, refAverage]
    rw [show (∑ c, weight β E₀ c * (-β * V c)) = -β * ∑ c, weight β E₀ c * V c by
      rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun c _ => by ring]
    field_simp
  have hpt : ∀ c, weight β E₀ c * (Real.exp m * (1 + (-β * V c) - m)) ≤
      weight β (fun c => E₀ c + V c) c := by
    intro c
    have h1 := Real.add_one_le_exp (-β * V c - m)
    have hsplit : weight β (fun c => E₀ c + V c) c =
        weight β E₀ c * (Real.exp m * Real.exp (-β * V c - m)) := by
      simp only [weight, ← Real.exp_add]; ring_nf
    rw [hsplit]
    refine mul_le_mul_of_nonneg_left ?_ (weight_pos β E₀ c).le
    exact mul_le_mul_of_nonneg_left (by linarith) (Real.exp_pos m).le
  have hlhs : Z β E₀ * Real.exp m =
      ∑ c, weight β E₀ c * (Real.exp m * (1 + (-β * V c) - m)) := by
    have : ∑ c, weight β E₀ c * (Real.exp m * (1 + (-β * V c) - m)) =
        Real.exp m * (Z β E₀ + ∑ c, weight β E₀ c * (-β * V c) - m * Z β E₀) := by
      simp only [Z, Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun c _ => by ring
    rw [this, hsum]; ring
  rw [hlhs]
  exact Finset.sum_le_sum fun c _ => hpt c

/-- **Gibbs–Bogoliubov upper bound on `F`**: `F(E₀ + V) ≤ F(E₀) + ⟨V⟩₀` for `0 < β`. -/
theorem freeEnergy_le_gibbs_bogoliubov [Nonempty C] {β : ℝ} (hβ : 0 < β) (E₀ V : C → ℝ) :
    freeEnergy β (fun c => E₀ c + V c) ≤ freeEnergy β E₀ + refAverage β E₀ V := by
  have h := Z_ge_gibbs_bogoliubov β E₀ V
  have hpos : 0 < Z β E₀ * Real.exp (-β * refAverage β E₀ V) :=
    mul_pos (Z_pos β E₀) (Real.exp_pos _)
  have hlog := Real.log_le_log hpos h
  rw [Real.log_mul (Z_pos β E₀).ne' (Real.exp_pos _).ne', Real.log_exp] at hlog
  have e : freeEnergy β E₀ + refAverage β E₀ V =
      -(Real.log (Z β E₀) + -β * refAverage β E₀ V) / β := by
    rw [freeEnergy]; field_simp; ring
  rw [e, freeEnergy]
  exact div_le_div_of_nonneg_right (by linarith) hβ.le

/-- The reference average never exceeds a uniform upper bound on the perturbation, so the
Gibbs–Bogoliubov bound is at least as strong as the worst-case bound. -/
theorem refAverage_le [Nonempty C] (β : ℝ) (E₀ V : C → ℝ) {Vhi : ℝ} (h : ∀ c, V c ≤ Vhi) :
    refAverage β E₀ V ≤ Vhi := by
  rw [refAverage, div_le_iff₀ (Z_pos β E₀), Z, Finset.mul_sum]
  exact Finset.sum_le_sum fun c _ => by
    rw [mul_comm Vhi]; exact mul_le_mul_of_nonneg_left (h c) (weight_pos β E₀ c).le

/-- **Compositional Gibbs–Bogoliubov bound.** With the decoupled two-domain landscape as
reference, `Z ≥ Z₁ Z₂ · exp (-β ⟨J⟩₁₂)` where `⟨J⟩₁₂` is the average coupling in the
product ensemble.  By `refAverage_le` this is never weaker than `Z ≥ Z₁ Z₂ e^{-β Jhi}`. -/
theorem Z_coupledEnergy_ge_gibbs {C₁ C₂ : Type*} [Fintype C₁] [Fintype C₂] [Nonempty C₁]
    [Nonempty C₂] (β : ℝ) (E₁ : C₁ → ℝ) (E₂ : C₂ → ℝ) (J : C₁ → C₂ → ℝ) :
    Z β E₁ * Z β E₂ *
        Real.exp (-β * refAverage β (coupledEnergy E₁ E₂ fun _ _ => 0) (fun p => J p.1 p.2)) ≤
      Z β (coupledEnergy E₁ E₂ J) := by
  have h := Z_ge_gibbs_bogoliubov β (coupledEnergy E₁ E₂ fun _ _ => 0) (fun p => J p.1 p.2)
  rw [Z_coupledEnergy_zero] at h
  convert h using 2
  funext p
  simp [coupledEnergy]

end FoldBounds
