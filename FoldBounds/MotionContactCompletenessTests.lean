import FoldBounds.MotionContactCompleteness

open FoldBounds.MotionUncertainty
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionContactCompleteness

private def e01 : Spring 3 := ⟨0, 1, ![1, 0, 0]⟩
private def e10 : Spring 3 := ⟨1, 0, ![1, 0, 0]⟩
private def e02 : Spring 3 := ⟨0, 2, ![1, 0, 0]⟩
private def complete : Fin 1 → Spring 3 := ![e01]
private def missing : Fin 0 → Spring 3 := ![]
private def duplicated : Fin 2 → Spring 3 := ![e01, e01]
private def reversed : Fin 1 → Spring 3 := ![e10]
private def spurious : Fin 1 → Spring 3 := ![e02]
private def x : Fin 3 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0], ![20, 0, 0]]
private def collision : Fin 3 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![0, 0, 0], ![20, 0, 0]]
private def zero : Fin 3 × Fin 3 → ℚ := fun _ => 0
private def wide : ScaledInterval :=
  ⟨0, 10 * (2 : ℤ) ^ 96⟩

example : checkContactNetwork complete x 100 = true := by native_decide
example : checkContactNetwork missing x 100 = false := by native_decide
example : checkContactNetwork duplicated x 100 = false := by native_decide
example : checkContactNetwork reversed x 100 = false := by native_decide
example : checkContactNetwork spurious x 100 = false := by native_decide
example : checkContactNetwork complete collision 100 = false := by native_decide
example : checkContactNetwork complete x 0 = false := by native_decide

-- The old sparse score formula would still define a model if a contact
-- were silently deleted. The new full-network gate rejects that model.
example : checkCompleteSparseScore missing x 100 (1 / 5)
    e01 zero wide = false := by native_decide
example : checkCompleteSparseScore duplicated x 100 (1 / 5)
    e01 zero wide = false := by native_decide

-- This accepted score refers to exactly the coordinate-implied graph,
-- not merely an arbitrary user-selected spring subset.
example : checkCompleteSparseScore complete x 100 (1 / 5)
    e01 zero wide = true := by native_decide

example :
    ExactContactNetwork complete x 100 ∧
      ∃! y : Fin 3 × Fin 3 → ℚ,
        matVecQ (proteinHessian complete x 100 (1 / 5)) y =
          pairForce e01 x ∧
        lowerQ wide ≤ pairScore e01 x y ∧
        pairScore e01 x y ≤ upperQ wide := by
  have h : checkCompleteSparseScore complete x 100 (1 / 5)
      e01 zero wide = true := by native_decide
  exact checkCompleteSparseScore_sound complete x 100 (1 / 5)
    e01 zero wide h

-- At exact cutoff (distance 10 A), the contact must be present.
private def boundary : Fin 3 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![10, 0, 0], ![30, 0, 0]]
example : checkContactNetwork complete boundary 100 = true := by native_decide
example : checkContactNetwork missing boundary 100 = false := by native_decide
