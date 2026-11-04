import FoldBounds.MotionExactSolvability

open FoldBounds.MotionHessian
open FoldBounds.MotionResidual
open FoldBounds.MotionUncertainty
open FoldBounds.MotionVariableDirections
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionExactSolvability

private def b : Fin 1 → Fin 2 → ℚ := ![![1, -1]]
private def w : Fin 1 → ℚ := ![2]
private def z : Fin 2 → ℚ := ![1, 0]
private def candidate : Fin 2 → ℚ := ![1 / 2, 0]

example : ∃! y : Fin 2 → ℚ, matVecQ (hessianQ b w 1) y = z := by
  apply rankOne_solvable b w 1
  · norm_num
  · intro i
    fin_cases i
    norm_num

example : ∃! y : Fin 2 → ℚ,
    matVecQ (hessianQ b w 1) y = z ∧
      (1 / 4 : ℚ) ≤ dotQ z y ∧ dotQ z y ≤ 3 / 2 := by
  have hc : checkResponse b w 1 z candidate = some (1 / 4, 3 / 2) := by
    native_decide
  exact checkResponse_exists_unique_sound b w 1 z candidate
    (1 / 4, 3 / 2) hc

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 2 := ![e]
private def coords : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
private def zero : (Fin 2 × Fin 3) → ℚ := fun _ => 0

/-- Proof of an exact finite protein-motion response without a claimed
inverse, external exact solution or trusted eigensolver. -/
example : ∃! y : (Fin 2 × Fin 3) → ℚ,
    matVecQ (proteinHessian edges coords 100 (1 / 5)) y =
      pairForce e coords ∧
    (0 : ℚ) ≤ pairScore e coords y ∧ pairScore e coords y ≤ 10 := by
  have hc : checkPairScore edges coords 100 (1 / 5) e zero =
    some (0, 10) := by native_decide
  exact checkPairScore_exists_unique_sound edges coords
    100 (1 / 5) e zero (0, 10) hc
