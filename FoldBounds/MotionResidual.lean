module

public import FoldBounds.MotionProteinHessian

/-!
# Residual certificates for exact finite rational Hessians

This closes the *linear-algebraic* residual gap for the explicit
ridge-regularized rank-one Hessian built in MotionHessian.

The output checker uses only rational arithmetic and rejects nonpositive
ridge coefficients and negative springs. If an exact solution y exists
to H y = z, a successful check proves that z.dot(y) lies between the
returned variational lower and residual-based upper bound.

Inverting an actual floating-point matrix is NOT trusted here. Solvability
is an explicit theorem premise; proving it from ridge positivity and
producing verified physical input bytes remain subsequent tasks.
-/

@[expose] public section
namespace FoldBounds.MotionResidual

open FoldBounds.MotionHessian

variable {I : Type*} [Fintype I] [DecidableEq I]

def vsub (u v : I → ℚ) : I → ℚ := fun i => u i - v i

theorem dotQ_comm (u v : I → ℚ) :
    dotQ u v = dotQ v u := by
  unfold dotQ
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem dotQ_vsub_left (u v t : I → ℚ) :
    dotQ (vsub u v) t = dotQ u t - dotQ v t := by
  unfold dotQ vsub
  simp only [sub_mul, Finset.sum_sub_distrib]

theorem dotQ_vsub_right (t u v : I → ℚ) :
    dotQ t (vsub u v) = dotQ t u - dotQ t v := by
  unfold dotQ vsub
  simp only [mul_sub, Finset.sum_sub_distrib]

theorem matVecQ_vsub (H : Matrix I I ℚ) (u v : I → ℚ) :
    matVecQ H (vsub u v) =
      vsub (matVecQ H u) (matVecQ H v) := by
  funext i
  simp only [matVecQ, vsub, mul_sub, Finset.sum_sub_distrib]

theorem dotQ_self_nonneg (u : I → ℚ) : 0 ≤ dotQ u u := by
  unfold dotQ
  apply Finset.sum_nonneg
  intro i _
  nlinarith [sq_nonneg (u i)]

theorem dotQ_matVec_symmetric (H : Matrix I I ℚ)
    (hH : ∀ i j, H i j = H j i) (u v : I → ℚ) :
    dotQ u (matVecQ H v) = dotQ v (matVecQ H u) := by
  calc
    dotQ u (matVecQ H v) =
        ∑ i : I, ∑ j : I, u i * (H i j * v j) := by
          change (∑ i : I, u i * (∑ j : I, H i j * v j)) = _
          apply Finset.sum_congr rfl
          intro i _
          rw [Finset.mul_sum]
    _ = ∑ j : I, ∑ i : I, u i * (H i j * v j) := by
          rw [Finset.sum_comm]
    _ = ∑ j : I, v j * (∑ i : I, H j i * u i) := by
          apply Finset.sum_congr rfl
          intro j _
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro i _
          rw [hH i j]
          ring
    _ = dotQ v (matVecQ H u) := rfl

/-- Exact quadratic-completion identity for symmetric rational matrices. -/
theorem dual_gap_identity (H : Matrix I I ℚ)
    (hH : ∀ i j, H i j = H j i) (x y : I → ℚ) :
    dotQ (matVecQ H y) y -
      (2 * dotQ (matVecQ H y) x -
        dotQ x (matVecQ H x)) =
    dotQ (vsub y x) (matVecQ H (vsub y x)) := by
  have hxy := dotQ_matVec_symmetric H hH y x
  calc
    dotQ (matVecQ H y) y -
        (2 * dotQ (matVecQ H y) x -
          dotQ x (matVecQ H x)) =
      dotQ y (matVecQ H y) -
        2 * dotQ x (matVecQ H y) +
          dotQ x (matVecQ H x) := by
            rw [dotQ_comm (matVecQ H y) y]
            rw [dotQ_comm (matVecQ H y) x]
            ring
    _ = dotQ (vsub y x) (matVecQ H (vsub y x)) := by
          rw [matVecQ_vsub]
          simp only [dotQ_vsub_left, dotQ_vsub_right]
          rw [hxy]
          ring

/-- Nonnegativity of a whole squared residual, expanded exactly. -/
theorem residual_square_ineq (r d : I → ℚ) (τ : ℚ) :
    0 ≤ dotQ r r - 2 * τ * dotQ r d +
      τ ^ 2 * dotQ d d := by
  have hn : 0 ≤ ∑ i : I, (r i - τ * d i) ^ 2 := by
    apply Finset.sum_nonneg
    intro i _
    exact sq_nonneg _
  have hs : (∑ i : I, (r i - τ * d i) ^ 2) =
      dotQ r r - 2 * τ * dotQ r d + τ ^ 2 * dotQ d d := by
    calc
      (∑ i : I, (r i - τ * d i) ^ 2) =
          ∑ i : I, (r i * r i - 2 * τ * (r i * d i) +
            τ ^ 2 * (d i * d i)) := by
              apply Finset.sum_congr rfl
              intro i _
              ring
      _ = dotQ r r - 2 * τ * dotQ r d +
            τ ^ 2 * dotQ d d := by
              simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib]
              change (∑ i : I, r i * r i) -
                (∑ i : I, 2 * τ * (r i * d i)) +
                  (∑ i : I, τ ^ 2 * (d i * d i)) =
                (∑ i : I, r i * r i) -
                   2 * τ * (∑ i : I, r i * d i) +
                     τ ^ 2 * (∑ i : I, d i * d i)
              rw [Finset.mul_sum, Finset.mul_sum]
  rw [hs] at hn
  exact hn

/-- Residual response bounds for a symmetric τ-coercive rational matrix.
The exact solution appears only in the mathematical theorem, not in the
data required by the numerical checker. -/
theorem symmetric_residual_enclosure (H : Matrix I I ℚ)
    (τ : ℚ) (z x y : I → ℚ)
    (hτ : 0 < τ)
    (hH : ∀ i j, H i j = H j i)
    (hcoercive : ∀ d : I → ℚ,
      τ * dotQ d d ≤ dotQ d (matVecQ H d))
    (hy : ∀ i, matVecQ H y i = z i) :
    (2 * dotQ z x - dotQ x (matVecQ H x)) ≤ dotQ z y ∧
      dotQ z y ≤
        (2 * dotQ z x - dotQ x (matVecQ H x)) +
          dotQ (vsub z (matVecQ H x))
            (vsub z (matVecQ H x)) / τ := by
  have hz : z = matVecQ H y := by
    funext i
    exact (hy i).symm
  subst z
  let d := vsub y x
  let r := vsub (matVecQ H y) (matVecQ H x)
  let E := dotQ d (matVecQ H d)
  let N := dotQ d d
  let R := dotQ r r
  have hgap : dotQ (matVecQ H y) y -
      (2 * dotQ (matVecQ H y) x -
        dotQ x (matVecQ H x)) = E := by
    simpa [d, E] using dual_gap_identity H hH x y
  have hr : r = matVecQ H d := by
    exact matVecQ_vsub H y x |>.symm
  have hrd : dotQ r d = E := by
    rw [dotQ_comm r d, hr]
  have hn : 0 ≤ N := dotQ_self_nonneg d
  have hc : τ * N ≤ E := hcoercive d
  have hE : 0 ≤ E := le_trans (mul_nonneg hτ.le hn) hc
  have hsq : 0 ≤ R - 2 * τ * E + τ ^ 2 * N := by
    have h := residual_square_ineq r d τ
    rw [hrd] at h
    exact h
  have hmul : τ * (τ * N) ≤ τ * E :=
    mul_le_mul_of_nonneg_left hc hτ.le
  have hR : τ * E ≤ R := by
    nlinarith [hsq, hmul]
  have hER : E ≤ R / τ := (le_div_iff₀ hτ).2 (by nlinarith [hR])
  constructor
  · dsimp [E] at hE
    linarith [hgap, hE]
  · linarith [hgap, hER]

/-- The exact finite rank-one Hessian satisfies all premises above. -/
theorem rankOne_residual_enclosure {m : ℕ}
    (b : Fin m → I → ℚ) (w : Fin m → ℚ) (τ : ℚ)
    (z x y : I → ℚ)
    (hτ : 0 < τ) (hw : ∀ e, 0 ≤ w e)
    (hy : ∀ i, matVecQ (hessianQ b w τ) y i = z i) :
    (2 * dotQ z x - dotQ x (matVecQ (hessianQ b w τ) x)) ≤
        dotQ z y ∧
    dotQ z y ≤
      (2 * dotQ z x - dotQ x (matVecQ (hessianQ b w τ) x)) +
        dotQ (vsub z (matVecQ (hessianQ b w τ) x))
          (vsub z (matVecQ (hessianQ b w τ) x)) / τ := by
  exact symmetric_residual_enclosure (hessianQ b w τ) τ z x y hτ
    (hessianQ_symmetric b w τ)
    (fun d => hessianQ_coercive b w τ d hw) hy

/-- Fail-closed executable exact rational residual certificate: proposals
may choose an approximate solve x, but cannot choose the endpoints. -/
def checkResponse {m : ℕ} (b : Fin m → I → ℚ) (w : Fin m → ℚ)
    (τ : ℚ) (z x : I → ℚ) : Option (ℚ × ℚ) :=
  if decide (0 < τ ∧ ∀ e, 0 ≤ w e) then
    let H := hessianQ b w τ
    let L := 2 * dotQ z x - dotQ x (matVecQ H x)
    let r := vsub z (matVecQ H x)
    some (L, L + dotQ r r / τ)
  else none

/-- A successful executable check is sound for every exact solution of
the specified rational system; the floating-point proposal is not trusted. -/
theorem checkResponse_sound {m : ℕ}
    (b : Fin m → I → ℚ) (w : Fin m → ℚ)
    (τ : ℚ) (z x y : I → ℚ) (bounds : ℚ × ℚ)
    (h : checkResponse b w τ z x = some bounds)
    (hy : ∀ i, matVecQ (hessianQ b w τ) y i = z i) :
    bounds.1 ≤ dotQ z y ∧ dotQ z y ≤ bounds.2 := by
  have hvalid : decide (0 < τ ∧ ∀ e, 0 ≤ w e) = true := by
    cases hh : decide (0 < τ ∧ ∀ e, 0 ≤ w e) with
    | false => simp [checkResponse, hh] at h
    | true => rfl
  have hs : 0 < τ ∧ ∀ e, 0 ≤ w e := of_decide_eq_true hvalid
  have heq :
      (2 * dotQ z x - dotQ x (matVecQ (hessianQ b w τ) x),
       2 * dotQ z x - dotQ x (matVecQ (hessianQ b w τ) x) +
         dotQ (vsub z (matVecQ (hessianQ b w τ) x))
           (vsub z (matVecQ (hessianQ b w τ) x)) / τ) =
      bounds := by
    simpa [checkResponse, hvalid] using h
  cases heq
  exact rankOne_residual_enclosure b w τ z x y hs.1 hs.2 hy

end FoldBounds.MotionResidual
