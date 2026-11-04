import FoldBounds.MotionPairOrdering

open FoldBounds.MotionPairOrdering
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionUncertainty

example : strictlyOrdered (0, 1) (2, 3) = true := by native_decide
example : strictlyOrdered (0, 1) (1, 3) = false := by native_decide
example : strictlyOrdered (0, 3) (2, 4) = false := by native_decide

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 2 := ![e]
private def x : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
private def zero : (Fin 2 × Fin 3) → ℚ := fun _ => 0

-- The exact same physical pair cannot be certified strictly below itself.
example : checkPairOrder edges x 100 (1 / 5) e e zero zero = false :=
  by native_decide
-- Invalid ridge must not produce a comparison.
example : checkPairOrder edges x 100 0 e e zero zero = false :=
  by native_decide
