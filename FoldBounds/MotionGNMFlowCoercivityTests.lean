import FoldBounds.MotionGNMFlowCoercivity

open FoldBounds.MotionHessian
open FoldBounds.MotionGNMFlowCoercivity

-- Two grounded vertices and two edges to the anchored (deleted) reference.
private def simpleIncidence : Fin 2 → Fin 2 → ℚ :=
  fun e i => if e = i then 1 else 0
private def exactFlow : Fin 2 → Fin 2 → ℚ :=
  fun i e => if i = e then 1 else 0
private def falseFlow : Fin 2 → Fin 2 → ℚ := fun _ _ => 0
private def force : Fin 2 → ℚ := fun i => if i = 0 then 1 else -1
private def zeroApprox : Fin 2 → ℚ := fun _ => 0

example : checkFlowCertificate simpleIncidence exactFlow = true := by native_decide
example : checkFlowCertificate simpleIncidence falseFlow = false := by native_decide
example : flowBound exactFlow = 2 := by norm_num [flowBound, exactFlow]
example : checkFlowResponse simpleIncidence exactFlow force zeroApprox = some (0,4) := by
  native_decide
example : checkFlowResponse simpleIncidence falseFlow force zeroApprox = none := by
  native_decide
example : ∃! y : Fin 2 → ℚ,
    matVecQ (hessianQ simpleIncidence (fun _ => 1) 0) y = force ∧
    (0 : ℚ) ≤ dotQ force y ∧ dotQ force y ≤ 4 := by
  have hc : checkFlowResponse simpleIncidence exactFlow force zeroApprox = some (0,4) := by
    native_decide
  exact checkFlowResponse_sound simpleIncidence exactFlow force zeroApprox (0,4) hc

#print axioms FoldBounds.MotionGNMFlowCoercivity.checkFlowIdentity_sound
#print axioms FoldBounds.MotionGNMFlowCoercivity.checkFlowCertificate_sound
#print axioms FoldBounds.MotionGNMFlowCoercivity.flow_reconstruct
#print axioms FoldBounds.MotionGNMFlowCoercivity.flow_coercivity
#print axioms FoldBounds.MotionGNMFlowCoercivity.flow_response_exists_unique
#print axioms FoldBounds.MotionGNMFlowCoercivity.checkFlowResponse_sound