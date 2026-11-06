import FoldBounds.MotionDyadicDivision

open FoldBounds.MotionDyadicDivision
open FoldBounds.MotionScaledCertificates

private def oneThird : ScaledInterval :=
  ⟨((2 : ℤ) ^ 96) / 3, ((2 : ℤ) ^ 96) / 3 + 1⟩

private def minusThird : ScaledInterval :=
  ⟨-(((2 : ℤ) ^ 96) / 3) - 1, -(((2 : ℤ) ^ 96) / 3)⟩

private def tooTight : ScaledInterval :=
  ⟨((2 : ℤ) ^ 96) / 3 + 1, ((2 : ℤ) ^ 96) / 3 + 1⟩

private def exactContact : ScaledInterval :=
  ⟨49 * ((2 : ℤ) ^ 90), 49 * ((2 : ℤ) ^ 90)⟩

example : checkQuotient 49 64 exactContact = true := by native_decide
example : checkQuotient 1 3 oneThird = true := by native_decide
example : checkQuotient (-1) 3 minusThird = true := by native_decide
example : checkQuotient 1 3 tooTight = false := by native_decide
example : checkQuotient 1 0 oneThird = false := by native_decide
example : checkQuotient 1 (-3) oneThird = false := by native_decide

example :
    lowerQ oneThird ≤ (1 : ℚ) / 3 ∧
      (1 : ℚ) / 3 ≤ upperQ oneThird := by
  have hc : checkQuotient 1 3 oneThird = true := by native_decide
  exact checkQuotient_sound 1 3 oneThird hc
