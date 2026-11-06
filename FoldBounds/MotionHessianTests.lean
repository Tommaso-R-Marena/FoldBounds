import FoldBounds.MotionHessian

open FoldBounds.MotionHessian

private def b : Fin 1 → Fin 2 → ℚ := ![![1, -1]]
private def w : Fin 1 → ℚ := ![2]
private def u : Fin 2 → ℚ := ![1, 0]

-- Exact analytic 2x2: H = [[3,-2],[-2,3]], uHu = 3.
example : hessianQ b w 1 0 0 = 3 := by native_decide
example : hessianQ b w 1 0 1 = -2 := by native_decide
example : hessianQ b w 1 1 0 = -2 := by native_decide
example : hessianQ b w 1 1 1 = 3 := by native_decide
example : matVecQ (hessianQ b w 1) u 0 = 3 := by native_decide
example : matVecQ (hessianQ b w 1) u 1 = -2 := by native_decide
example : springEnergyQ b w 1 u = 3 := by native_decide
example : dotQ u (matVecQ (hessianQ b w 1) u) = 3 := by native_decide

-- Exact universal mathematical results, not just numerical tests.
example (v : Fin 2 → ℚ) :
    hessianQ b w 1 0 1 = hessianQ b w 1 1 0 :=
  hessianQ_symmetric b w 1 0 1

example (v : Fin 2 → ℚ) :
    dotQ v (matVecQ (hessianQ b w 1) v) = springEnergyQ b w 1 v :=
  dot_matVec_hessian b w 1 v
