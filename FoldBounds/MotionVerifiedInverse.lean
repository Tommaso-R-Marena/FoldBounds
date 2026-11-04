module

public import FoldBounds.MotionPairOrdering

/-!
# FoldBounds.MotionVerifiedInverse

A finite rational *witness-checking* construction that discharges the
remaining existence premise of the model-relative motion score theorem.

The prover may supply any proposed inverse B; the executable gate verifies
H * B = I by exact rational arithmetic. This implies that y = B * z is an
exact solution of H * y = z for EVERY force vector z. The witness is
untrusted. A bad inverse fails closed.

Unlike taking existence as an external premise, the successful Option
certificate provides an actual exact solution and a score enclosure.
Checking the full inverse is cubic in matrix dimension and is intended
as an initial trustworthy baseline, not the final scalable implementation.
-/

@[expose] public section

namespace FoldBounds.MotionVerifiedInverse

open FoldBounds.MotionHessian
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionResidual
open FoldBounds.MotionUncertainty

variable {I : Type*} [Fintype I] [DecidableEq I]

/-- Finite exact rational right-inverse gate. -/
def checkRightInverse (H B : Matrix I I ℚ) : Bool :=
  decide (∀ i j, (∑ k : I, H i k * B k j) =
    if i = j then (1 : ℚ) else 0)

/-- Passing the executable gate means the untrusted matrix is genuinely
a right inverse on every coordinate. -/
theorem checkRightInverse_valid (H B : Matrix I I ℚ)
    (h : checkRightInverse H B = true) :
    ∀ i j, (∑ k : I, H i k * B k j) =
      if i = j then (1 : ℚ) else 0 := by
  exact of_decide_eq_true h

/-- Construct the exact solution for any force vector, with no existence
premise, once the rational inverse witness passes the gate. -/
theorem checked_inverse_solves (H B : Matrix I I ℚ)
    (h : checkRightInverse H B = true) (z : I → ℚ) :
    matVecQ H (matVecQ B z) = z := by
  have hh := checkRightInverse_valid H B h
  funext i
  calc
    matVecQ H (matVecQ B z) i =
        ∑ k : I, H i k * (∑ j : I, B k j * z j) := rfl
    _ = ∑ k : I, ∑ j : I, H i k * (B k j * z j) := by
          apply Finset.sum_congr rfl
          intro k _
          rw [Finset.mul_sum]
    _ = ∑ j : I, ∑ k : I, H i k * (B k j * z j) := by
          rw [Finset.sum_comm]
    _ = ∑ j : I, (∑ k : I, H i k * B k j) * z j := by
          apply Finset.sum_congr rfl
          intro j _
          rw [Finset.sum_mul]
          apply Finset.sum_congr rfl
          intro k _
          ring
    _ = ∑ j : I, (if i = j then (1 : ℚ) else 0) * z j := by
          apply Finset.sum_congr rfl
          intro j _
          rw [hh i j]
    _ = z i := by simp

/-- All submitted inverse entries are checked, not trusted. The subsequent
score bounds are independently constructed from the approximate solve. -/
def checkPairScoreWithInverse {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate : (Fin n × Fin 3) → ℚ)
    (inverseWitness :
      Matrix (Fin n × Fin 3) (Fin n × Fin 3) ℚ) :
    Option (ℚ × ℚ) :=
  if checkRightInverse (proteinHessian edges coords cutoffSq τ)
      inverseWitness then
    checkPairScore edges coords cutoffSq τ pair candidate
  else none

/-- An accepted certificate gives the exact response bounds for the
*constructed* rational exact solution B*z: no assumed solvability,
no trust in the numerical inverse witness, and no external oracle. -/
theorem checkPairScoreWithInverse_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate : (Fin n × Fin 3) → ℚ)
    (B : Matrix (Fin n × Fin 3) (Fin n × Fin 3) ℚ)
    (bounds : ℚ × ℚ)
    (h : checkPairScoreWithInverse edges coords cutoffSq τ
      pair candidate B = some bounds) :
    let exactY := matVecQ B (pairForce pair coords)
    bounds.1 ≤ pairScore pair coords exactY ∧
      pairScore pair coords exactY ≤ bounds.2 := by
  have hInv : checkRightInverse
      (proteinHessian edges coords cutoffSq τ) B = true := by
    cases hi : checkRightInverse
        (proteinHessian edges coords cutoffSq τ) B with
    | false => simp [checkPairScoreWithInverse, hi] at h
    | true => rfl
  have hScore : checkPairScore edges coords cutoffSq τ
      pair candidate = some bounds := by
    simpa [checkPairScoreWithInverse, hInv] using h
  have hy : ∀ i,
      matVecQ (proteinHessian edges coords cutoffSq τ)
        (matVecQ B (pairForce pair coords)) i = pairForce pair coords i := by
    intro i
    exact congrFun (checked_inverse_solves
      (proteinHessian edges coords cutoffSq τ)
      B hInv (pairForce pair coords)) i
  exact checkPairScore_sound edges coords cutoffSq τ pair candidate
    (matVecQ B (pairForce pair coords)) bounds hScore hy

end FoldBounds.MotionVerifiedInverse
