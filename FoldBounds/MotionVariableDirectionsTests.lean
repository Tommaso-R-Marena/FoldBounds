import FoldBounds.MotionVariableDirections

open FoldBounds.MotionVariableDirections
open FoldBounds.MotionCoordinates
open FoldBounds.MotionUncertainty

private def zero : Fin 3 → Box :=
  ![⟨0, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def near : Fin 3 → Box :=
  ![⟨7, 9⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def switch : Fin 3 → Box :=
  ![⟨9, 11⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def far : Fin 3 → Box :=
  ![⟨11, 13⟩, ⟨0, 0⟩, ⟨0, 0⟩]

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def disp : Fin 2 → Fin 3 → ℚ :=
  ![![1, 0, 0], ![0, 0, 0]]

example : lowerTerm e (![zero, near]) disp 100 1 = 49 := by native_decide
example : upperTerm e (![zero, near]) disp 100 1 = 81 := by native_decide
example : lowerTerm e (![zero, switch]) disp 100 1 = 0 := by native_decide
example : upperTerm e (![zero, switch]) disp 100 1 = 121 := by native_decide
example : lowerTerm e (![zero, far]) disp 100 1 = 0 := by native_decide
example : upperTerm e (![zero, far]) disp 100 1 = 0 := by native_decide

private def actual : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
example : energyTerm e actual disp 100 1 = 64 := by native_decide
example : lowerTerm e (![zero, near]) disp 100 1 ≤
    energyTerm e actual disp 100 1 := by native_decide
example : energyTerm e actual disp 100 1 ≤
    upperTerm e (![zero, near]) disp 100 1 := by native_decide
