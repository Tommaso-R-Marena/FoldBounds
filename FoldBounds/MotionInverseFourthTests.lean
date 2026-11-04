import FoldBounds.MotionInverseFourth

open FoldBounds.MotionInverseFourth
open FoldBounds.MotionCoordinates
open FoldBounds.MotionVariableDirections
open FoldBounds.MotionUncertainty

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 2 := ![e]
private def zero : Fin 3 → Box :=
  ![⟨0, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def near : Fin 3 → Box :=
  ![⟨7, 9⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def exactEight : Fin 3 → Box :=
  ![⟨8, 8⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def switching : Fin 3 → Box :=
  ![⟨9, 11⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def distant : Fin 3 → Box :=
  ![⟨11, 13⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def collision : Fin 3 → Box :=
  ![⟨-1, 1⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def reversed : Fin 3 → Box :=
  ![⟨9, 7⟩, ⟨0, 0⟩, ⟨0, 0⟩]
private def u : Fin 2 → Fin 3 → ℚ :=
  ![![1, 0, 0], ![0, 0, 0]]
private def x : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]

-- An exact 8 A contact has d²=64, stretch²=64, E=49/64.
example : weightedEnergyTerm e x u 100 = 49 / 64 := by native_decide
example : weightedLowerTerm e (![zero, exactEight]) u 100 = 49 / 64 := by native_decide
example : weightedUpperTerm e (![zero, exactEight]) u 100 = 49 / 64 := by native_decide
example : safeNetwork edges (![zero, exactEight]) 100 = true := by native_decide

-- A range encompassing the contact cutoff is admitted, with a zero
-- lower-energy endpoint and a positive upper-energy endpoint.
example : safeNetwork edges (![zero, switching]) 100 = true := by native_decide
example : weightedLowerTerm e (![zero, switching]) u 100 = 0 := by native_decide
example : 0 < weightedUpperTerm e (![zero, switching]) u 100 := by native_decide
example : weightedUpperTerm e (![zero, distant]) u 100 = 0 := by native_decide

-- All physical energy cases are enclosed by the frozen coordinate boxes.
example : weightedLowerTerm e (![zero, near]) u 100 ≤
    weightedEnergyTerm e x u 100 := by native_decide
example : weightedEnergyTerm e x u 100 ≤
    weightedUpperTerm e (![zero, near]) u 100 := by native_decide

-- Fail closed on possible collisions, malformed coordinate boxes,
-- nonpositive cutoffs, and contacts with zero separation.
example : safeNetwork edges (![zero, collision]) 100 = false := by native_decide
example : safeNetwork edges (![zero, reversed]) 100 = false := by native_decide
example : safeNetwork edges (![zero, exactEight]) 0 = false := by native_decide
example : safeNetwork edges (![zero, zero]) 100 = false := by native_decide


-- Exact model-level PSD/ridge regression for the same 8 A example.
example : weightedNetworkEnergy edges x u 100 = 49 / 64 := by native_decide
example : weightedRidgeEnergy edges x u 100 1 = 49 / 64 + 1 := by native_decide

-- Universal nonnegative-energy and coercivity proofs are kernel-checked,
-- not merely verified on the numerical example above.
example (xx uu : Fin 2 → Fin 3 → ℚ) (τ : ℚ) :
    τ * displacementSqQ uu ≤
      weightedRidgeEnergy edges xx uu 100 τ := by
  exact weightedRidgeEnergy_coercive edges xx uu 100 τ
