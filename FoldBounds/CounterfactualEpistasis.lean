import Mathlib

/-!
# An exact identity underlying ensemble-mediated epistasis

This is a generic rational algebra lemma, NOT a proof of any protein-specific
free-energy model, physical conformer distribution, or the Python microstate
enumerator. It identifies the covariance cross term when two mutation factors
are multiplicative in *each* conformation (no direct pair interaction).
-/
namespace FoldBounds.CounterfactualEpistasis

/-- Even when individual mutation factors multiply within each conformer,
    their *ensemble averages* need not multiply. The difference is an exact
    weighted covariance. For two microstates and arbitrary rational weights,
    the interaction can be positive or negative. -/
theorem two_state_mutation_covariance
    (w₁ w₂ r₁ r₂ s₁ s₂ : ℚ) :
    (w₁ + w₂) * (w₁ * r₁ * s₁ + w₂ * r₂ * s₂) -
      (w₁ * r₁ + w₂ * r₂) * (w₁ * s₁ + w₂ * s₂) =
      w₁ * w₂ * (r₁ - r₂) * (s₁ - s₂) := by
  ring

/-- Conditional covariance may be strictly positive if the factors are
    ordered consistently in two conformations of positive weights. -/
theorem two_state_positive_covariance
    (w₁ w₂ r₁ r₂ s₁ s₂ : ℚ)
    (h₁ : 0 < w₁) (h₂ : 0 < w₂)
    (hr : r₂ < r₁) (hs : s₂ < s₁) :
    0 < (w₁ + w₂) * (w₁ * r₁ * s₁ + w₂ * r₂ * s₂) -
      (w₁ * r₁ + w₂ * r₂) * (w₁ * s₁ + w₂ * s₂) := by
  rw [two_state_mutation_covariance]
  have h_r : 0 < r₁ - r₂ := sub_pos.mpr hr
  have h_s : 0 < s₁ - s₂ := sub_pos.mpr hs
  exact mul_pos (mul_pos (mul_pos h₁ h₂) h_r) h_s

/-- Anti-correlated mutation reweighting can give strictly negative
    population-level epistasis with no pairwise state energy term. -/
theorem two_state_negative_covariance
    (w₁ w₂ r₁ r₂ s₁ s₂ : ℚ)
    (h₁ : 0 < w₁) (h₂ : 0 < w₂)
    (hr : r₂ < r₁) (hs : s₁ < s₂) :
    (w₁ + w₂) * (w₁ * r₁ * s₁ + w₂ * r₂ * s₂) -
      (w₁ * r₁ + w₂ * r₂) * (w₁ * s₁ + w₂ * s₂) < 0 := by
  rw [two_state_mutation_covariance]
  have h_r : 0 < r₁ - r₂ := sub_pos.mpr hr
  have ha : 0 < w₁ * w₂ * (r₁ - r₂) :=
    mul_pos (mul_pos h₁ h₂) h_r
  have hb : s₁ - s₂ < 0 := sub_neg.mpr hs
  exact mul_neg_of_pos_of_neg ha hb

/-- An energy model that receives only a coarse amino-acid encoding
    cannot distinguish sequences with identical encodings, in any state.
    This is an identifiability constraint, not a theorem about real protein
    conformations or a particular biochemical classification. -/
theorem encoded_energy_indistinguishable
    {α σ : Type} (hp : α → Bool) (energy : List Bool → σ → ℚ)
    (xs ys : List α) (conf : σ)
    (h : xs.map hp = ys.map hp) :
    energy (xs.map hp) conf = energy (ys.map hp) conf := by
  rw [h]

end FoldBounds.CounterfactualEpistasis