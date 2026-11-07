import Mathlib

/-!
# Generic stationarity of proof-gated reversible kinetic mixtures

These theorems establish the abstract finite-state identity used by the
Python exact Markov witnesses. This is **not** a Lean refinement of the
Python transition decoder, pivot geometry, max-flow search, or kinetics.
-/

namespace FoldBounds.JepaKinetics

open scoped BigOperators

variable {α : Type*} [Fintype α]

/-- Detailed balance and stochastic rows imply stationarity of `w`.
    Unlike a convergence theorem, this does not assume irreducibility. -/
theorem stationary_of_detailedBalance
    (w : α → ℚ) (P : α → α → ℚ)
    (rowSum : ∀ x, (∑ y, P x y) = 1)
    (balance : ∀ x y, w x * P x y = w y * P y x) :
    ∀ y, (∑ x, w x * P x y) = w y := by
  intro y
  calc
    (∑ x, w x * P x y) = ∑ x, w y * P y x := by
      apply Finset.sum_congr rfl
      intro x _
      exact balance x y
    _ = w y * (∑ x, P y x) := by rw [Finset.mul_sum]
    _ = w y := by rw [rowSum y]; ring

/-- Convex combinations of rational row-stochastic matrices remain
    row-stochastic (algebraically for any rational mixing weight). -/
theorem mix_rowSum
    (P Q : α → α → ℚ) (a : ℚ)
    (rowP : ∀ x, (∑ y, P x y) = 1)
    (rowQ : ∀ x, (∑ y, Q x y) = 1) :
    ∀ x, (∑ y, ((1-a)*P x y+a*Q x y)) = 1 := by
  intro x
  calc
    (∑ y, ((1-a)*P x y+a*Q x y))
      = (1-a)*(∑ y, P x y)+a*(∑ y, Q x y) := by
        rw [Finset.sum_add_distrib, Finset.mul_sum, Finset.mul_sum]
    _ = 1 := by rw [rowP x, rowQ x]; ring

/-- Detailed balance is preserved by mixing any two kernels targeting
    the same stationary weight function. -/
theorem mix_detailedBalance
    (w : α → ℚ) (P Q : α → α → ℚ) (a : ℚ)
    (balanceP : ∀ x y, w x * P x y = w y * P y x)
    (balanceQ : ∀ x y, w x * Q x y = w y * Q y x) :
    ∀ x y,
      w x * ((1-a)*P x y+a*Q x y)
      = w y * ((1-a)*P y x+a*Q y x) := by
  intro x y
  calc
    w x * ((1-a)*P x y+a*Q x y)
      = (1-a)*(w x*P x y)+a*(w x*Q x y) := by ring
    _ = (1-a)*(w y*P y x)+a*(w y*Q y x) := by
      rw [balanceP x y, balanceQ x y]
    _ = w y * ((1-a)*P y x+a*Q y x) := by ring

/-- The checked mixture used by FoldBounds therefore preserves `w`. -/
theorem mix_stationary
    (w : α → ℚ) (P Q : α → α → ℚ) (a : ℚ)
    (rowP : ∀ x, (∑ y, P x y) = 1)
    (rowQ : ∀ x, (∑ y, Q x y) = 1)
    (balanceP : ∀ x y, w x * P x y = w y * P y x)
    (balanceQ : ∀ x y, w x * Q x y = w y * Q y x) :
    ∀ y, (∑ x, w x*((1-a)*P x y+a*Q x y)) = w y :=
  stationary_of_detailedBalance w
    (fun x y => (1-a)*P x y+a*Q x y)
    (mix_rowSum P Q a rowP rowQ)
    (mix_detailedBalance w P Q a balanceP balanceQ)

end FoldBounds.JepaKinetics
