module

public import FoldBounds.MotionGNMAdaptiveCases

@[expose] public section
namespace FoldBounds.MotionGNMAdaptiveCasesTests

open FoldBounds.MotionGNMAdaptiveCases

-- Every allowed world is sandwiched whether a contact is selected or not.
example : caseLowWeight (fun _ : Fin 2 => false)
      (fun _ => true) (fun _ => false)
      (selectedState (fun _ => false) (fun _ => true)) 0 = 0 := by
  decide

example : caseHighWeight (fun _ : Fin 2 => false)
      (fun _ => true) (fun _ => false)
      (selectedState (fun _ => false) (fun _ => true)) 0 = 1 := by
  decide

-- This case of the theorem is conditional on score-envelope premises; the
-- proof itself is not replacing real dataset interval validation.
example : caseLowWeight (fun _ : Fin 1 => true)
      (fun _ => false) (fun _ => false)
      (selectedState (fun _ => false) (fun _ => false)) 0 = 1 := by
  decide

end FoldBounds.MotionGNMAdaptiveCasesTests