module

public import FoldBounds.FreeEnergy

/-!
# FoldBounds.Compositional — compositional (domain-decomposition) bounds

Conformations of a two-domain system are pairs `(x, y) ∈ C₁ × C₂` with energy
`E(x, y) = E₁ x + E₂ y + J x y`, where `J` is the inter-domain coupling.
If the coupling is certified to lie in `[Jlo, Jhi]` then, for `0 ≤ β`,

`Z₁ Z₂ e^{-β Jhi} ≤ Z ≤ Z₁ Z₂ e^{-β Jlo}`   and, for `0 < β`,
`F₁ + F₂ + Jlo ≤ F ≤ F₁ + F₂ + Jhi`.

Only `|C₁| + |C₂|` energy evaluations are needed instead of `|C₁| · |C₂|`; the price is a
multiplicative uncertainty `e^{β (Jhi - Jlo)}` on `Z`.  For `J = 0` the bound is the exact
factorisation `Z = Z₁ Z₂`.
-/

@[expose] public section

open Finset Real

namespace FoldBounds

variable {C₁ C₂ : Type*} [Fintype C₁] [Fintype C₂]

/-- Two-domain energy `E(x, y) = E₁ x + E₂ y + J x y`. -/
def coupledEnergy (E₁ : C₁ → ℝ) (E₂ : C₂ → ℝ) (J : C₁ → C₂ → ℝ) : C₁ × C₂ → ℝ :=
  fun p => E₁ p.1 + E₂ p.2 + J p.1 p.2

variable (β : ℝ) (E₁ : C₁ → ℝ) (E₂ : C₂ → ℝ) (J : C₁ → C₂ → ℝ)

theorem Z_coupledEnergy :
    Z β (coupledEnergy E₁ E₂ J) =
      ∑ x, ∑ y, weight β E₁ x * weight β E₂ y * Real.exp (-β * J x y) := by
  unfold Z
  rw [Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun y _ => ?_
  simp only [weight, coupledEnergy, ← Real.exp_add]
  ring_nf

private theorem Z_mul_Z_mul (k : ℝ) :
    Z β E₁ * Z β E₂ * k = ∑ x, ∑ y, weight β E₁ x * weight β E₂ y * k := by
  unfold Z
  rw [Finset.sum_mul_sum, Finset.sum_mul]
  exact Finset.sum_congr rfl fun x _ => Finset.sum_mul _ _ _

/-- **Exact factorisation** for non-interacting domains (any real `β`). -/
theorem Z_coupledEnergy_zero :
    Z β (coupledEnergy E₁ E₂ fun _ _ => 0) = Z β E₁ * Z β E₂ := by
  rw [Z_coupledEnergy, ← mul_one (Z β E₁ * Z β E₂), Z_mul_Z_mul]
  simp

/-- **Compositional partition-function bound.** -/
theorem Z_coupledEnergy_mem_Icc {β : ℝ} (hβ : 0 ≤ β) {Jlo Jhi : ℝ}
    (hJ : ∀ x y, J x y ∈ Set.Icc Jlo Jhi) :
    Z β (coupledEnergy E₁ E₂ J) ∈
      Set.Icc (Z β E₁ * Z β E₂ * Real.exp (-β * Jhi))
        (Z β E₁ * Z β E₂ * Real.exp (-β * Jlo)) := by
  rw [Z_coupledEnergy, Z_mul_Z_mul, Z_mul_Z_mul]
  have hw : ∀ x y, 0 ≤ weight β E₁ x * weight β E₂ y :=
    fun x y => (mul_pos (weight_pos β E₁ x) (weight_pos β E₂ y)).le
  constructor
  · refine Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ => ?_
    refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (hw x y)
    nlinarith [(hJ x y).2]
  · refine Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ => ?_
    refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (hw x y)
    nlinarith [(hJ x y).1]

/-- **Compositional free-energy bound**: `F₁ + F₂ + Jlo ≤ F ≤ F₁ + F₂ + Jhi` for `0 < β`. -/
theorem freeEnergy_coupledEnergy_mem_Icc [Nonempty C₁] [Nonempty C₂] {β : ℝ} (hβ : 0 < β)
    {Jlo Jhi : ℝ} (hJ : ∀ x y, J x y ∈ Set.Icc Jlo Jhi) :
    freeEnergy β (coupledEnergy E₁ E₂ J) ∈
      Set.Icc (freeEnergy β E₁ + freeEnergy β E₂ + Jlo)
        (freeEnergy β E₁ + freeEnergy β E₂ + Jhi) := by
  have hZ := Z_coupledEnergy_mem_Icc E₁ E₂ J hβ.le hJ
  have h1 := Z_pos β E₁
  have h2 := Z_pos β E₂
  have hlo : 0 < Z β E₁ * Z β E₂ * Real.exp (-β * Jhi) := by positivity
  have hI := freeEnergy_mem_Icc (coupledEnergy E₁ E₂ J) hβ hlo hZ
  have hlog : ∀ t, -Real.log (Z β E₁ * Z β E₂ * Real.exp (-β * t)) / β =
      freeEnergy β E₁ + freeEnergy β E₂ + t := by
    intro t
    rw [Real.log_mul (by positivity) (by positivity), Real.log_mul h1.ne' h2.ne',
      Real.log_exp, freeEnergy, freeEnergy]
    field_simp
    ring
  rw [hlog, hlog] at hI
  exact hI

end FoldBounds
