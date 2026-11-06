module

public import FoldBounds.MotionGNMGeometricNogoods

@[expose] public section
namespace FoldBounds.MotionGNMGeometricNogoods

/-- A small, executable contradictory two-row contact-box example. -/
example : checkLinearNoGood (1/20 : ℚ)
    (fun _ : Fin 1 => (0 : ℚ)) (-1 : ℚ) = true := by
  norm_num [checkLinearNoGood]

/-- An exact-rational example proving a necessary contact-on linear gate. -/
example :
    2 * (∑ j : Fin 1, (1 : ℚ) * (0 : ℚ)) ≤
      (2 : ℚ) - (∑ _j : Fin 1, (1 : ℚ) ^ 2) := by
  norm_num

end FoldBounds.MotionGNMGeometricNogoods