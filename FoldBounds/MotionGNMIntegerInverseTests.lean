import FoldBounds.MotionGNMIntegerInverse

open FoldBounds.MotionGNMIntegerInverse

private def chainA : Matrix (Fin 2) (Fin 2) ℤ :=
  ![![1, -1], ![-1, 2]]
private def chainInverse : Matrix (Fin 2) (Fin 2) ℤ :=
  ![![2, 1], ![1, 1]]
private def wrongInverse : Matrix (Fin 2) (Fin 2) ℤ :=
  ![![2, 1], ![1, 2]]
private def triangleA : Matrix (Fin 2) (Fin 2) ℤ :=
  ![![2, -1], ![-1, 2]]

example : checkScaledInverse chainA chainInverse 1 = true := by
  native_decide
example : checkScaledInverse chainA wrongInverse 1 = false := by
  native_decide
example : checkScaledInverse triangleA wrongInverse 3 = true := by
  native_decide
example : checkScaledInverse chainA chainInverse 0 = false := by
  native_decide

#print axioms FoldBounds.MotionGNMIntegerInverse.checkScaledInverse_valid
#print axioms FoldBounds.MotionGNMIntegerInverse.checkScaledInverse_to_rational
#print axioms FoldBounds.MotionGNMIntegerInverse.checkScaledInverse_rightInverse
#print axioms FoldBounds.MotionGNMIntegerInverse.checkScaledInverse_solves