import FoldBounds.MotionScaledCertificates

open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionUncertainty
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionHessian

-- The independent 2^96 outward-rounded interval for 320/277.
private def good : ScaledInterval :=
  ⟨91527119150052664367993011218, 91527119150052664367993011219⟩
private def badLower : ScaledInterval :=
  ⟨91527119150052664367993011219, 91527119150052664367993011219⟩
private def badUpper : ScaledInterval :=
  ⟨91527119150052664367993011218, 91527119150052664367993011218⟩

example : checkValue good (320 / 277) = true := by native_decide
example : checkValue badLower (320 / 277) = false := by native_decide
example : checkValue badUpper (320 / 277) = false := by native_decide
example : checkEnclosure good (320 / 277) (320 / 277) = true :=
  by native_decide

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 2 := ![e]
private def x : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
private def exactY : Fin 2 → Fin 3 → ℚ :=
  ![![-1280 / 277, 0, 0], ![1280 / 277, 0, 0]]

-- The actual rational 3D ANM model with ridge τ=1/5 is checked, not a
-- free-standing arithmetic example. A malicious tighter interval fails.
example : checkScaledPairScore edges x 100 (1 / 5) e
    (flat exactY) good = true := by native_decide
example : checkScaledPairScore edges x 100 (1 / 5) e
    (flat exactY) badLower = false := by native_decide
example : checkScaledPairScore edges x 100 (1 / 5) e
    (flat exactY) badUpper = false := by native_decide
example : checkScaledPairScore edges x 100 0 e
    (flat exactY) good = false := by native_decide

-- An accepted 2^96-scaled certificate establishes exact existence,
-- uniqueness, and enclosure for the physical score, without
-- trusting either the proposer or an inverse solver.
example :
    ∃! y : (Fin 2 × Fin 3) → ℚ,
      matVecQ (proteinHessian edges x 100 (1 / 5)) y =
        pairForce e x ∧
      lowerQ good ≤ pairScore e x y ∧ pairScore e x y ≤ upperQ good := by
  have hc : checkScaledPairScore edges x 100 (1 / 5)
      e (flat exactY) good = true := by native_decide
  exact checkScaledPairScore_exists_unique_sound edges x
    100 (1 / 5) e (flat exactY) good hc
