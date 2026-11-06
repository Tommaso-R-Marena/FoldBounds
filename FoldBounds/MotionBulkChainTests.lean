import FoldBounds.MotionBulkChain

open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionBulkChain

private def a : ScaledInterval := ⟨10, 20⟩
private def b : ScaledInterval := ⟨21, 30⟩
private def c : ScaledInterval := ⟨100, 200⟩
private def tie : ScaledInterval := ⟨20, 30⟩
private def overlap : ScaledInterval := ⟨19, 31⟩

example : checkAdjacent ([] : List ScaledInterval) = true := by native_decide
example : checkAdjacent [a] = true := by native_decide
example : checkAdjacent [a,b,c] = true := by native_decide
example : checkAdjacent [a,tie] = false := by native_decide
example : checkAdjacent [a,overlap] = false := by native_decide
example : checkAdjacent [c,b,a] = false := by native_decide

/-- A three-element chain has no hidden unverified comparison. -/
example : StrictAdjacent [a,b,c] := by
  have h : checkAdjacent [a,b,c] = true := by native_decide
  exact checkAdjacent_sound [a,b,c] h

/-- Exact integer endpoint separation implies a certified real order,
assuming each separate physical-model score certificate was accepted. -/
example (x y : ℚ)
    (hx : lowerQ a ≤ x ∧ x ≤ upperQ a)
    (hy : lowerQ b ≤ y ∧ y ≤ upperQ b) : x < y := by
  have hab : separated a b = true := by native_decide
  exact scores_ordered_of_certificates a b x y hx hy hab
