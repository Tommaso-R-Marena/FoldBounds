module

public import FoldBounds.MotionVerifiedInverse

/-!
# FoldBounds.MotionExactSolvability

Coercivity of the explicit rational rank-one Hessian implies that its
finite-dimensional matrix action is injective and hence bijective.
Consequently every rational forcing vector admits a unique exact
rational solution of H*y=z.

This removes the existence premise from the residual/score theorem
*without* requiring a cubic-size submitted inverse witness. The
external numerical solver remains untrusted; the physical inputs
are still rational and fixed.
-/

@[expose] public section
namespace FoldBounds.MotionExactSolvability

open FoldBounds.MotionHessian
open FoldBounds.MotionResidual
open FoldBounds.MotionUncertainty
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse

variable {I : Type*} [Fintype I] [DecidableEq I]

theorem dotQ_self_eq_zero_iff (u : I → ℚ) :
    dotQ u u = 0 ↔ ∀ i, u i = 0 := by
  constructor
  · intro hz i
    have hs : (u i) ^ 2 ≤ (∑ j : I, (u j) ^ 2) := by
      apply Finset.single_le_sum
      · intro j _
        exact sq_nonneg (u j)
      · exact Finset.mem_univ i
    have hsum : (∑ j : I, (u j) ^ 2) = 0 := by
      calc
        (∑ j : I, (u j) ^ 2) =
            (∑ j : I, u j * u j) := by
              apply Finset.sum_congr rfl
              intro j _
              ring
        _ = 0 := hz
    have hsq : (u i) ^ 2 ≤ 0 := by
      simpa [hsum] using hs
    nlinarith [sq_nonneg (u i)]
  · intro h
    unfold dotQ
    simp [h]

/-- Every symmetric τ-coercive rational square matrix has a trivial kernel. -/
theorem coercive_matrix_injective (H : Matrix I I ℚ) (τ : ℚ)
    (hτ : 0 < τ)
    (hcoercive : ∀ d : I → ℚ,
      τ * dotQ d d ≤ dotQ d (matVecQ H d)) :
    Function.Injective H.mulVec := by
  intro a b hab
  have heq : matVecQ H a = matVecQ H b := hab
  let d := vsub a b
  have hz : ∀ i, matVecQ H d i = 0 := by
    intro i
    have hdiff := congrFun (matVecQ_vsub H a b) i
    rw [hdiff]
    exact sub_eq_zero.mpr (congrFun heq i)
  have hzero : dotQ d (matVecQ H d) = 0 := by
    unfold dotQ
    simp [hz]
  have hnorm : 0 ≤ dotQ d d := dotQ_self_nonneg d
  have hnormzero : dotQ d d = 0 := by
    have hc := hcoercive d
    nlinarith
  have hd : ∀ i, d i = 0 :=
    (dotQ_self_eq_zero_iff d).mp hnormzero
  funext i
  have hsub : a i - b i = 0 := hd i
  exact sub_eq_zero.mp hsub

/-- Solvability follows from coercivity using the standard finite matrix
equivalence between injectivity and surjectivity of its action. -/
theorem coercive_matrix_solvable (H : Matrix I I ℚ) (τ : ℚ)
    (hτ : 0 < τ)
    (hcoercive : ∀ d : I → ℚ,
      τ * dotQ d d ≤ dotQ d (matVecQ H d))
    (z : I → ℚ) :
    ∃! y : I → ℚ, matVecQ H y = z := by
  have hi : Function.Injective H.mulVec :=
    coercive_matrix_injective H τ hτ hcoercive
  have hu : IsUnit H := (Matrix.mulVec_injective_iff_isUnit).mp hi
  have hs : Function.Surjective H.mulVec :=
    (Matrix.mulVec_surjective_iff_isUnit).mpr hu
  obtain ⟨y, hy⟩ := hs z
  refine ⟨y, hy, ?_⟩
  intro y' hy'
  exact hi (hy'.trans hy.symm)

/-- A strictly positive rational ridge and nonnegative finite spring
weights guarantee a unique exact solution for every rational force. -/
theorem rankOne_solvable {m : ℕ}
    (b : Fin m → I → ℚ) (w : Fin m → ℚ) (τ : ℚ)
    (hτ : 0 < τ) (hw : ∀ e, 0 ≤ w e)
    (z : I → ℚ) :
    ∃! y : I → ℚ, matVecQ (hessianQ b w τ) y = z := by
  exact coercive_matrix_solvable (hessianQ b w τ) τ hτ
    (fun d => hessianQ_coercive b w τ d hw) z

/-- An accepted residual certificate bounds a provably existing, unique
exact model response. No inverse witness or external existence assumption. -/
theorem checkResponse_exists_unique_sound {m : ℕ}
    (b : Fin m → I → ℚ) (w : Fin m → ℚ)
    (τ : ℚ) (z candidate : I → ℚ) (bounds : ℚ × ℚ)
    (h : checkResponse b w τ z candidate = some bounds) :
    ∃! y : I → ℚ,
      matVecQ (hessianQ b w τ) y = z ∧
      bounds.1 ≤ dotQ z y ∧ dotQ z y ≤ bounds.2 := by
  have hg : decide (0 < τ ∧ ∀ e, 0 ≤ w e) = true := by
    cases ht : decide (0 < τ ∧ ∀ e, 0 ≤ w e) with
    | false => simp [checkResponse, ht] at h
    | true => rfl
  obtain ⟨ht, hw⟩ := of_decide_eq_true hg
  obtain ⟨y, hy, huniq⟩ := rankOne_solvable b w τ ht hw z
  have hb := checkResponse_sound b w τ z candidate y bounds h
    (fun i => congrFun hy i)
  refine ⟨y, ⟨hy, hb.1, hb.2⟩, ?_⟩
  intro y' hy'
  exact huniq y' hy'.1

/-- At fixed rational coordinates, a successful protein pair certificate
encloses the score of a uniquely existing exact rational solution. -/
theorem checkPairScore_exists_unique_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate : (Fin n × Fin 3) → ℚ) (bounds : ℚ × ℚ)
    (h : checkPairScore edges coords cutoffSq τ
      pair candidate = some bounds) :
    ∃! y : (Fin n × Fin 3) → ℚ,
      matVecQ (proteinHessian edges coords cutoffSq τ) y =
        pairForce pair coords ∧
      bounds.1 ≤ pairScore pair coords y ∧
      pairScore pair coords y ≤ bounds.2 := by
  have hgate : decide (0 < cutoffSq ∧ 0 < distSq pair coords ∧
      ∀ e : Fin m, 0 < distSq (edges e) coords) = true := by
    cases hg : decide (0 < cutoffSq ∧ 0 < distSq pair coords ∧
        ∀ e : Fin m, 0 < distSq (edges e) coords) with
    | false => simp [checkPairScore, hg] at h
    | true => rfl
  let b : Fin m → (Fin n × Fin 3) → ℚ :=
    fun e => incidence (edges e) coords
  let w : Fin m → ℚ :=
    fun e => contactWeight (edges e) coords cutoffSq
  have hr : checkResponse b w τ (pairForce pair coords)
      candidate ≠ none := by
    intro hn
    simp [checkPairScore, hgate, b, w, hn] at h
  have hτ : 0 < τ := by
    cases hval : checkResponse b w τ (pairForce pair coords)
        candidate with
    | none => exact (hr hval).elim
    | some raw =>
      have hvalid : decide (0 < τ ∧ ∀ e, 0 ≤ w e) = true := by
        cases hs : decide (0 < τ ∧ ∀ e, 0 ≤ w e) with
        | false => simp [checkResponse, hs] at hval
        | true => rfl
      exact (of_decide_eq_true hvalid).1
  have hw : ∀ e, 0 ≤ w e := by
    intro e
    exact contactWeight_nonneg (edges e) coords cutoffSq
  obtain ⟨y, hy, huniq⟩ :=
    rankOne_solvable b w τ hτ hw (pairForce pair coords)
  have hb := checkPairScore_sound edges coords cutoffSq τ pair
    candidate y bounds h (fun i => congrFun hy i)
  refine ⟨y, ⟨hy, hb.1, hb.2⟩, ?_⟩
  intro y' hy'
  exact huniq y' hy'.1

end FoldBounds.MotionExactSolvability
