import FoldBounds.MotionProteinHessian

open FoldBounds.MotionHessian
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionUncertainty
open FoldBounds.MotionInverseFourth

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 2 := ![e]
private def x : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
private def u : Fin 2 → Fin 3 → ℚ :=
  ![![1, 0, 0], ![0, 0, 0]]

-- The diagonal spring coefficient is exactly 49/64 at an 8 A contact.
example : proteinHessian edges x 100 1 (0, 0) (0, 0) =
    1 + 49 / 64 := by native_decide
example : proteinHessian edges x 100 1 (0, 0) (1, 0) =
    -(49 / 64) := by native_decide
example : proteinHessian edges x 100 1 (0, 1) (0, 1) = 1 := by native_decide

example : weightedRidgeEnergy edges x u 100 1 = 1 + 49 / 64 :=
  by native_decide
example : dotQ (flat u)
    (matVecQ (proteinHessian edges x 100 1) (flat u)) =
    1 + 49 / 64 := by native_decide

example (v : Fin 2 → Fin 3 → ℚ) :
    dotQ (flat v)
      (matVecQ (proteinHessian edges x 100 1) (flat v)) =
      weightedRidgeEnergy edges x v 100 1 :=
  proteinHessian_energy edges x v 100 1

example (v : Fin 2 → Fin 3 → ℚ) :
    dotQ (flat v) (flat v) ≤
      dotQ (flat v)
        (matVecQ (proteinHessian edges x 100 1) (flat v)) :=
  proteinHessian_coercive edges x v 100 1
