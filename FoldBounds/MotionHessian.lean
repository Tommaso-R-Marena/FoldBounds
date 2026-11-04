module

public import FoldBounds.MotionInverseFourth

/-!
# FoldBounds.MotionHessian

Finite rational Hessian of a ridge-regularized rank-one elastic network.
For contact incidences b, spring coefficients w, and ridge τ, the explicit
matrix Hᵢⱼ = ∑ₑ wₑ bₑᵢ bₑⱼ + τ δᵢⱼ
is symmetric and has the expected sum-of-squares quadratic form.

This establishes matrix/source energy equivalence for rank-one incidences.
Connection to the protein Cartesian incidence construction follows separately.
-/

@[expose] public section

namespace FoldBounds.MotionHessian

variable {I : Type*} [Fintype I] [DecidableEq I]

def dotQ (v u : I → ℚ) : ℚ :=
  ∑ i : I, v i * u i

def matVecQ (H : Matrix I I ℚ) (u : I → ℚ) (i : I) : ℚ :=
  ∑ j : I, H i j * u j

def hessianQ {m : ℕ} (b : Fin m → I → ℚ)
    (w : Fin m → ℚ) (τ : ℚ) : Matrix I I ℚ :=
  fun i j => (∑ e : Fin m, w e * b e i * b e j) +
    if i = j then τ else 0

def springEnergyQ {m : ℕ} (b : Fin m → I → ℚ)
    (w : Fin m → ℚ) (τ : ℚ) (u : I → ℚ) : ℚ :=
  (∑ e : Fin m, w e * (dotQ (b e) u) ^ 2) + τ * dotQ u u

theorem hessianQ_symmetric {m : ℕ} (b : Fin m → I → ℚ)
    (w : Fin m → ℚ) (τ : ℚ) (i j : I) :
    hessianQ b w τ i j = hessianQ b w τ j i := by
  have hs : (∑ e : Fin m, w e * b e i * b e j) =
      (∑ e : Fin m, w e * b e j * b e i) := by
    apply Finset.sum_congr rfl
    intro e _
    ring
  by_cases heq : i = j
  · subst j
    rfl
  · have hne : j ≠ i := Ne.symm heq
    simp [hessianQ, heq, hne, hs]

theorem matVecQ_hessian {m : ℕ} (b : Fin m → I → ℚ)
    (w : Fin m → ℚ) (τ : ℚ) (u : I → ℚ) (i : I) :
    matVecQ (hessianQ b w τ) u i =
      (∑ e : Fin m, (w e * b e i) * dotQ (b e) u) +
        τ * u i := by
  have hs :
      (∑ j : I, (∑ e : Fin m, w e * b e i * b e j) * u j) =
        ∑ e : Fin m, (w e * b e i) * dotQ (b e) u := by
    calc
      (∑ j : I, (∑ e : Fin m, w e * b e i * b e j) * u j) =
          ∑ j : I, ∑ e : Fin m, (w e * b e i * b e j) * u j := by
            apply Finset.sum_congr rfl
            intro j _
            rw [Finset.sum_mul]
      _ = ∑ e : Fin m, ∑ j : I, (w e * b e i * b e j) * u j := by
            rw [Finset.sum_comm]
      _ = ∑ e : Fin m, (w e * b e i) * dotQ (b e) u := by
            apply Finset.sum_congr rfl
            intro e _
            rw [Finset.mul_sum]
            apply Finset.sum_congr rfl
            intro j _
            ring
  have hid : (∑ j : I, (if i = j then τ else 0) * u j) =
      τ * u i := by
    simp
  unfold matVecQ hessianQ
  simp only [add_mul, Finset.sum_add_distrib]
  rw [hs, hid]

theorem dot_matVec_hessian {m : ℕ} (b : Fin m → I → ℚ)
    (w : Fin m → ℚ) (τ : ℚ) (u : I → ℚ) :
    dotQ u (matVecQ (hessianQ b w τ) u) =
      springEnergyQ b w τ u := by
  have hcontact :
      (∑ i : I, u i *
        (∑ e : Fin m, (w e * b e i) * dotQ (b e) u)) =
      ∑ e : Fin m, w e * (dotQ (b e) u) ^ 2 := by
    calc
      (∑ i : I, u i *
        (∑ e : Fin m, (w e * b e i) * dotQ (b e) u)) =
          ∑ i : I, ∑ e : Fin m,
            u i * ((w e * b e i) * dotQ (b e) u) := by
              apply Finset.sum_congr rfl
              intro i _
              rw [Finset.mul_sum]
      _ = ∑ e : Fin m, ∑ i : I,
           u i * ((w e * b e i) * dotQ (b e) u) := by
             rw [Finset.sum_comm]
      _ = ∑ e : Fin m, w e * (dotQ (b e) u) ^ 2 := by
            apply Finset.sum_congr rfl
            intro e _
            calc
              (∑ i : I, u i * ((w e * b e i) * dotQ (b e) u)) =
                  (w e * dotQ (b e) u) * dotQ (b e) u := by
                    rw [Finset.mul_sum]
                    apply Finset.sum_congr rfl
                    intro i _
                    ring
              _ = w e * (dotQ (b e) u) ^ 2 := by ring
  have hridge : (∑ i : I, u i * (τ * u i)) = τ * dotQ u u := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    ring
  change (∑ i : I, u i * matVecQ (hessianQ b w τ) u i) =
      (∑ e : Fin m, w e * (dotQ (b e) u) ^ 2) + τ * dotQ u u
  simp_rw [matVecQ_hessian]
  simp only [mul_add, Finset.sum_add_distrib]
  rw [hcontact, hridge]

theorem hessianQ_coercive {m : ℕ}
    (b : Fin m → I → ℚ) (w : Fin m → ℚ)
    (τ : ℚ) (u : I → ℚ) (hw : ∀ e, 0 ≤ w e) :
    τ * dotQ u u ≤ dotQ u (matVecQ (hessianQ b w τ) u) := by
  rw [dot_matVec_hessian]
  unfold springEnergyQ
  have hc : 0 ≤ ∑ e : Fin m, w e * (dotQ (b e) u) ^ 2 := by
    apply Finset.sum_nonneg
    intro e _
    exact mul_nonneg (hw e) (sq_nonneg _)
  linarith

end FoldBounds.MotionHessian
