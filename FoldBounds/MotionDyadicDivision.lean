module

public import FoldBounds.MotionRoundedOps

/-!
# FoldBounds.MotionDyadicDivision

A checked integer division gate for the last exact-rational rounding primitive
in the fast frozen-protein fixed-point Hessian construction. An untrusted
proposer supplies signed numerator, positive denominator, and proposed
2^-96 dyadic outward integer endpoints. The checker accepts only if
cross-multiplied INTEGER inequalities establish exact enclosure.

In particular, zero/negative denominators and incorrectly rounded quotients
are rejected. The result does not assume a floating point implementation,
nor trust an external floor or ceil algorithm.
-/

@[expose] public section
namespace FoldBounds.MotionDyadicDivision

open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionRoundedOps

/-- Integer-only gate for the dyadic enclosure of signed num / positive den. -/
def checkQuotient (num den : ℤ) (out : ScaledInterval) : Bool :=
  decide (0 < den ∧
    out.lower * den ≤ num * scaleInt ∧
    num * scaleInt ≤ out.upper * den)

/-- Any accepted integer quotient interval encloses the exact rational
value num/den, even for negative numerators and non-exact dyadic results. -/
theorem checkQuotient_sound (num den : ℤ) (out : ScaledInterval)
    (h : checkQuotient num den out = true) :
    lowerQ out ≤ (num : ℚ) / (den : ℚ) ∧
      (num : ℚ) / (den : ℚ) ≤ upperQ out := by
  have hg : 0 < den ∧
      out.lower * den ≤ num * scaleInt ∧
      num * scaleInt ≤ out.upper * den :=
    of_decide_eq_true h
  have hden : (0 : ℚ) < (den : ℚ) := by exact_mod_cast hg.1
  have ha : ((out.lower * den : ℤ) : ℚ) ≤
      ((num * scaleInt : ℤ) : ℚ) := by exact_mod_cast hg.2.1
  have hb : ((num * scaleInt : ℤ) : ℚ) ≤
      ((out.upper * den : ℤ) : ℚ) := by exact_mod_cast hg.2.2
  have hlo : (out.lower : ℚ) * (den : ℚ) ≤
      (num : ℚ) * binaryScale := by
    simpa [scaleInt_cast] using ha
  have hhi : (num : ℚ) * binaryScale ≤
      (out.upper : ℚ) * (den : ℚ) := by
    simpa [scaleInt_cast] using hb
  constructor
  · unfold lowerQ
    exact (div_le_div_iff₀ binaryScale_pos hden).2 hlo
  · unfold upperQ
    exact (div_le_div_iff₀ hden binaryScale_pos).2 hhi

end FoldBounds.MotionDyadicDivision
