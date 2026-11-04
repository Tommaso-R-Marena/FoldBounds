import FoldBounds.MotionCoordinates

open FoldBounds.MotionCoordinates
open FoldBounds.MotionContactChecker
open FoldBounds.MotionUncertainty

private def pointBox : Fin 3 → Box :=
  ![⟨0, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def xBox : Fin 3 → Box :=
  ![⟨7, 9⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def farBox : Fin 3 → Box :=
  ![⟨11, 13⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def switchBox : Fin 3 → Box :=
  ![⟨9, 11⟩, ⟨0, 0⟩, ⟨0, 0⟩]

example : (squaredDistance pointBox xBox).lower = 49 := by native_decide
example : (squaredDistance pointBox xBox).upper = 81 := by native_decide
example : (squaredDistance pointBox farBox).lower = 121 := by native_decide
example : (squaredDistance pointBox farBox).upper = 169 := by native_decide
example : (squaredDistance pointBox switchBox).lower = 81 := by native_decide
example : (squaredDistance pointBox switchBox).upper = 121 := by native_decide

private def demoEdges : Fin 3 → Spring 4 :=
  ![⟨0, 1, ![1, 0, 0]⟩, ⟨0, 2, ![1, 0, 0]⟩,
    ⟨0, 3, ![1, 0, 0]⟩]
private def demoBoxes : Fin 4 → Fin 3 → Box :=
  ![pointBox, xBox, switchBox, farBox]
private def demoContacts : RawContacts 3 :=
  contactBoxes demoEdges demoBoxes 100 (![1, 1, 1]) (![2, 2, 2])

example : accepts demoContacts = true := by native_decide
example : lowerWeight demoContacts 0 = 1 := by native_decide
example : upperWeight demoContacts 0 = 2 := by native_decide
example : lowerWeight demoContacts 1 = 0 := by native_decide
example : upperWeight demoContacts 1 = 2 := by native_decide
example : lowerWeight demoContacts 2 = 0 := by native_decide
example : upperWeight demoContacts 2 = 0 := by native_decide

example : (squared (⟨-2, 3⟩ : Box)).lower = 0 := by native_decide
example : (squared (⟨-2, 3⟩ : Box)).upper = 9 := by native_decide
