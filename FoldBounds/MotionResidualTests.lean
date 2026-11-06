import FoldBounds.MotionResidual

open FoldBounds.MotionHessian
open FoldBounds.MotionResidual

private def b : Fin 1 → Fin 2 → ℚ := ![![1, -1]]
private def w : Fin 1 → ℚ := ![2]
private def y : Fin 2 → ℚ := ![3 / 5, 2 / 5]
private def z : Fin 2 → ℚ := ![1, 0]
private def x : Fin 2 → ℚ := ![1 / 2, 0]

-- Exact 2x2 linear system H = [[3,-2],[-2,3]].
example : matVecQ (hessianQ b w 1) y = z := by native_decide
example : dotQ z y = 3 / 5 := by native_decide
example : checkResponse b w 1 z y = some (3 / 5, 3 / 5) := by native_decide
example : checkResponse b w 1 z x = some (1 / 4, 3 / 2) := by native_decide

-- Negative stiffness, zero ridge and negative ridge must be rejected.
example : checkResponse b w 0 z x = none := by native_decide
example : checkResponse b w (-1) z x = none := by native_decide
example : checkResponse b (![-1]) 1 z x = none := by native_decide

-- A successful certificate really encloses the fixed exact solution.
example : (1 / 4 : ℚ) ≤ dotQ z y ∧ dotQ z y ≤ 3 / 2 := by
  have hy : ∀ i, matVecQ (hessianQ b w 1) y i = z i := by
    intro i
    have hh : matVecQ (hessianQ b w 1) y = z := by native_decide
    exact congrFun hh i
  have hc : checkResponse b w 1 z x = some (1 / 4, 3 / 2) := by native_decide
  exact checkResponse_sound b w 1 z x y (1 / 4, 3 / 2) hc hy
