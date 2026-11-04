import FoldBounds.MotionRoundedOps

open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionRoundedOps

private def a : ScaledInterval := ⟨1, 2⟩
private def good : ScaledInterval := ⟨0, 1⟩
private def corrupted : ScaledInterval := ⟨1, 1⟩

-- y=2^96 / 4 => y / 2^96 = 1/4.
-- [1,2]/2^96 multiplied by 1/4 lies inside [0,1]/2^96.
example : checkMul a good (scaleInt / 4) = true := by native_decide
example : checkMul a corrupted (scaleInt / 4) = false := by native_decide

-- Multiply by +1 or -1 represented as signed dyadic integers.
example : checkMul a a scaleInt = true := by native_decide
example : checkMul a (⟨-2, -1⟩) (-scaleInt) = true :=
  by native_decide
example : checkMul a (⟨-1, -1⟩) (-scaleInt) = false :=
  by native_decide

-- Directly instantiate the universal proof of accepted outward rounding:
-- signed values, no floating-point approximation.
example (x : ℚ) (ha : contains a x) :
    contains good (x * (((scaleInt / 4 : ℤ) : ℚ) / binaryScale)) := by
  have hc : checkMul a good (scaleInt / 4) = true := by native_decide
  exact checkMul_sound a good x (scaleInt / 4) ha hc

example (u v : ℚ) (hu : contains a u) (hv : contains a v) :
    contains (addBox a a) (u + v) :=
  addBox_sound a a u v hu hv

example (u : ℚ) (hu : contains a u) :
    contains (negBox a) (-u) :=
  negBox_sound a u hu
