import FoldBounds.MotionGNMSpectralSOS

open FoldBounds.MotionHessian
open FoldBounds.MotionGNMSpectralSOS

-- A = I + [1,-1]ᵀ[1,-1], i.e. 2-vertex connected grounded matrix.
private def testA : Matrix (Fin 2) (Fin 2) ℚ :=
  fun i j => if i = j then 2 else -1
private def testFactor : Fin 1 → Fin 2 → ℚ :=
  fun _ i => if i = 0 then 1 else -1
private def testWeights : Fin 1 → ℚ := fun _ => 1
private def force : Fin 2 → ℚ := fun i => if i = 0 then 1 else -1
private def zero : Fin 2 → ℚ := fun _ => 0

example : checkSpectralSOS testA 1 testFactor testWeights = true := by native_decide
example : checkSpectralSOS testA 2 testFactor testWeights = false := by native_decide
example : checkSpectralSOS testA 1 testFactor (fun _ => -1) = false := by native_decide
example : checkSpectralResponse testA 1 testFactor testWeights force zero = some (0,2) := by
  native_decide
example : ∃! y : Fin 2 → ℚ,
  matVecQ testA y = force ∧ (0:ℚ) ≤ dotQ force y ∧ dotQ force y ≤ 2 := by
  have h : checkSpectralResponse testA 1 testFactor testWeights force zero = some (0,2) := by
    native_decide
  exact checkSpectralResponse_sound testA 1 testFactor testWeights force zero (0,2) h

#print axioms FoldBounds.MotionGNMSpectralSOS.checkSpectralSOS_sound
#print axioms FoldBounds.MotionGNMSpectralSOS.checkSpectralSOS_coercive
#print axioms FoldBounds.MotionGNMSpectralSOS.checkSpectralSOS_solvable
#print axioms FoldBounds.MotionGNMSpectralSOS.checkSpectralResponse_sound