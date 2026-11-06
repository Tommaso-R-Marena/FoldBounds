import FoldBounds.MotionVerifiedInverse

open FoldBounds.MotionHessian
open FoldBounds.MotionUncertainty
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionVerifiedInverse

private def b : Fin 1 → Fin 2 → ℚ := ![![1, -1]]
private def w : Fin 1 → ℚ := ![2]
private def inv2 : Matrix (Fin 2) (Fin 2) ℚ :=
  ![![3 / 5, 2 / 5], ![2 / 5, 3 / 5]]

example : checkRightInverse (hessianQ b w 1) inv2 = true := by native_decide
example : matVecQ (hessianQ b w 1) (matVecQ inv2 ![1, 0]) =
    ![1, 0] := by native_decide
example : checkRightInverse (hessianQ b w 1)
    (fun _ _ => 0) = false := by native_decide

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 2 := ![e]
private def x : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
private def exactY : Fin 2 → Fin 3 → ℚ :=
  ![![-1280 / 277, 0, 0], ![1280 / 277, 0, 0]]

/-- The full inverse of the exact six-dimensional physical Hessian at
ridge 1/5, not a claim of invertibility from numerical eigensolvers. -/
private def inverse6 : Matrix (Fin 2 × Fin 3) (Fin 2 × Fin 3) ℚ :=
  fun p q =>
    if p.2 = 0 ∧ q.2 = 0 then
      if p.1 = q.1 then 1545 / 554 else 1225 / 554
    else if p = q then 5 else 0

private def zeroInverse : Matrix (Fin 2 × Fin 3) (Fin 2 × Fin 3) ℚ :=
  fun _ _ => 0

example : checkRightInverse (proteinHessian edges x 100 (1 / 5))
    inverse6 = true := by native_decide
example : checkRightInverse (proteinHessian edges x 100 (1 / 5))
    zeroInverse = false := by native_decide

example : checkPairScoreWithInverse edges x 100 (1 / 5)
    e (flat exactY) inverse6 =
    some (320 / 277, 320 / 277) := by native_decide
example : checkPairScoreWithInverse edges x 100 (1 / 5)
    e (flat exactY) zeroInverse = none := by native_decide

-- This does not take a separately supplied exact solution hypothesis.
example :
    let y := matVecQ inverse6 (pairForce e x)
    (320 / 277 : ℚ) ≤ pairScore e x y ∧
      pairScore e x y ≤ (320 / 277 : ℚ) := by
  have h : checkPairScoreWithInverse edges x 100 (1 / 5)
      e (flat exactY) inverse6 = some (320 / 277, 320 / 277) :=
    by native_decide
  exact checkPairScoreWithInverse_sound edges x 100 (1 / 5)
    e (flat exactY) inverse6 (320 / 277, 320 / 277) h
