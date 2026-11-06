import FoldBounds.MotionProteinResponse

open FoldBounds.MotionUncertainty
open FoldBounds.MotionVariableDirections
open FoldBounds.MotionHessian
open FoldBounds.MotionInverseFourth
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 2 := ![e]
private def x : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
private def zero : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![0, 0, 0]]
-- τ = 1/5 and the spring effective difference mode eigenvalue = 277/160.
private def exactY : Fin 2 → Fin 3 → ℚ :=
  ![![-1280 / 277, 0, 0], ![1280 / 277, 0, 0]]
private def brokenX : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![0, 0, 0]]

example : distSq e x = 64 := by native_decide
example : pairScore e x (flat exactY) = 320 / 277 := by native_decide
example : checkPairScore edges x 100 (1 / 5) e (flat exactY) =
    some (320 / 277, 320 / 277) := by native_decide

-- A zero approximate solution is safe but yields a wider enclosure.
example : checkPairScore edges x 100 (1 / 5) e (flat zero) =
    some (0, 10) := by native_decide

-- Fail closed on denominator zero, ridge zero, invalid contact cutoff.
example : checkPairScore edges brokenX 100 (1 / 5) e (flat zero) =
    none := by native_decide
example : checkPairScore edges x 100 0 e (flat zero) =
    none := by native_decide
example : checkPairScore edges x 0 (1 / 5) e (flat zero) =
    none := by native_decide

-- The computed end-to-end certificate proves the exact response in a
-- nontrivial finite 3D test, rather than trusting the proposed solution.
example :
    (320 / 277 : ℚ) ≤ pairScore e x (flat exactY) ∧
      pairScore e x (flat exactY) ≤ (320 / 277 : ℚ) := by
  have hy : ∀ i, matVecQ (proteinHessian edges x 100 (1 / 5))
      (flat exactY) i = pairForce e x i := by
    intro i
    have heq : matVecQ (proteinHessian edges x 100 (1 / 5))
        (flat exactY) = pairForce e x := by native_decide
    exact congrFun heq i
  have hc : checkPairScore edges x 100 (1 / 5) e (flat exactY) =
      some (320 / 277, 320 / 277) := by native_decide
  exact checkPairScore_sound edges x 100 (1 / 5) e
    (flat exactY) (flat exactY) (320 / 277, 320 / 277) hc hy
