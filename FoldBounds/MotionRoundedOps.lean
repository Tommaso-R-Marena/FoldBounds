module

public import FoldBounds.MotionScaledCertificates

/-!
# FoldBounds.MotionRoundedOps

Exact verification of the two primitive fixed-point interval operations
used by the untrusted bulk motion checker. All interval endpoints are
signed integers representing endpoint / 2^96.

Addition and negation are exact. For a product of an interval-valued
rational x by an integer witness y/2^96, the proposer supplies an
outward-rounded integer interval. A small *integer Boolean gate* checks
the rounding inequalities. The subsequent Lean theorem proves that any
accepted interval contains the true rational product.

No floating-point arithmetic, Python rounding code, or external
numerical solver is in the trusted acceptance proof.
-/

@[expose] public section
namespace FoldBounds.MotionRoundedOps

open FoldBounds.MotionScaledCertificates

def scaleInt : ℤ := 2 ^ 96

theorem scaleInt_cast : (scaleInt : ℚ) = binaryScale := by
  norm_num [scaleInt, binaryScale]

theorem scaleInt_pos : 0 < scaleInt := by
  norm_num [scaleInt]

def contains (a : ScaledInterval) (x : ℚ) : Prop :=
  (a.lower : ℚ) ≤ binaryScale * x ∧
  binaryScale * x ≤ (a.upper : ℚ)

def addBox (a b : ScaledInterval) : ScaledInterval :=
  ⟨a.lower + b.lower, a.upper + b.upper⟩

theorem addBox_sound (a b : ScaledInterval) (x y : ℚ)
    (ha : contains a x) (hb : contains b y) :
    contains (addBox a b) (x + y) := by
  rcases ha with ⟨ha1, ha2⟩
  rcases hb with ⟨hb1, hb2⟩
  change ((a.lower + b.lower : ℤ) : ℚ) ≤ binaryScale * (x + y) ∧
    binaryScale * (x + y) ≤ ((a.upper + b.upper : ℤ) : ℚ)
  constructor <;> push_cast <;> nlinarith

def negBox (a : ScaledInterval) : ScaledInterval :=
  ⟨-a.upper, -a.lower⟩

theorem negBox_sound (a : ScaledInterval) (x : ℚ)
    (ha : contains a x) : contains (negBox a) (-x) := by
  rcases ha with ⟨ha1, ha2⟩
  change ((-a.upper : ℤ) : ℚ) ≤ binaryScale * (-x) ∧
    binaryScale * (-x) ≤ ((-a.lower : ℤ) : ℚ)
  constructor <;> push_cast <;> nlinarith

/-- Check the proposer-supplied result for multiplying an interval by
a signed integer witness whose true rational value is y / 2^96.

Every comparison is on integers; floor/ceil semantics are enforced by
the acceptance inequalities, not assumed from the proposer. -/
def checkMul (a out : ScaledInterval) (y : ℤ) : Bool :=
  if 0 ≤ y then
    decide (out.lower * scaleInt ≤ a.lower * y ∧
      a.upper * y ≤ out.upper * scaleInt)
  else
    decide (out.lower * scaleInt ≤ a.upper * y ∧
      a.lower * y ≤ out.upper * scaleInt)

/-- Exact, checker-enforced signed 2^96 multiplication soundness. -/
theorem checkMul_sound (a out : ScaledInterval) (x : ℚ) (y : ℤ)
    (ha : contains a x) (hc : checkMul a out y = true) :
    contains out (x * ((y : ℚ) / binaryScale)) := by
  have hs : 0 < binaryScale := binaryScale_pos
  have hid : binaryScale * (x * ((y : ℚ) / binaryScale)) =
      x * (y : ℚ) := by
    field_simp [ne_of_gt hs]
    ring
  rcases ha with ⟨ha1, ha2⟩
  have hbounds :
      (out.lower : ℚ) * binaryScale ≤
        (if 0 ≤ y then (a.lower : ℚ) * (y : ℚ)
          else (a.upper : ℚ) * (y : ℚ)) ∧
      (if 0 ≤ y then (a.upper : ℚ) * (y : ℚ)
          else (a.lower : ℚ) * (y : ℚ)) ≤
        (out.upper : ℚ) * binaryScale := by
    by_cases hy : 0 ≤ y
    · have hh : out.lower * scaleInt ≤ a.lower * y ∧
          a.upper * y ≤ out.upper * scaleInt :=
        of_decide_eq_true (by simpa [checkMul, hy] using hc)
      have h1 : (out.lower : ℚ) * binaryScale ≤
          (a.lower : ℚ) * (y : ℚ) := by
        have hv : ((out.lower * scaleInt : ℤ) : ℚ) ≤
            ((a.lower * y : ℤ) : ℚ) := by exact_mod_cast hh.1
        simpa [scaleInt_cast] using hv
      have h2 : (a.upper : ℚ) * (y : ℚ) ≤
          (out.upper : ℚ) * binaryScale := by
        have hv : ((a.upper * y : ℤ) : ℚ) ≤
            ((out.upper * scaleInt : ℤ) : ℚ) := by exact_mod_cast hh.2
        simpa [scaleInt_cast] using hv
      simpa [hy] using And.intro h1 h2
    · have hh : out.lower * scaleInt ≤ a.upper * y ∧
          a.lower * y ≤ out.upper * scaleInt :=
        of_decide_eq_true (by simpa [checkMul, hy] using hc)
      have h1 : (out.lower : ℚ) * binaryScale ≤
          (a.upper : ℚ) * (y : ℚ) := by
        have hv : ((out.lower * scaleInt : ℤ) : ℚ) ≤
            ((a.upper * y : ℤ) : ℚ) := by exact_mod_cast hh.1
        simpa [scaleInt_cast] using hv
      have h2 : (a.lower : ℚ) * (y : ℚ) ≤
          (out.upper : ℚ) * binaryScale := by
        have hv : ((a.lower * y : ℤ) : ℚ) ≤
            ((out.upper * scaleInt : ℤ) : ℚ) := by exact_mod_cast hh.2
        simpa [scaleInt_cast] using hv
      simpa [hy] using And.intro h1 h2
  have hprod :
      (out.lower : ℚ) * binaryScale ≤
        (x * (y : ℚ)) * binaryScale ∧
      (x * (y : ℚ)) * binaryScale ≤
        (out.upper : ℚ) * binaryScale := by
    by_cases hy : 0 ≤ y
    · have hyq : (0 : ℚ) ≤ (y : ℚ) := by exact_mod_cast hy
      have hm1 := mul_le_mul_of_nonneg_right ha1 hyq
      have hm2 := mul_le_mul_of_nonneg_right ha2 hyq
      have hh1 := hbounds.1
      have hh2 := hbounds.2
      simp only [if_pos hy] at hh1 hh2
      constructor <;> nlinarith
    · have hyn : y ≤ 0 := le_of_lt (lt_of_not_ge hy)
      have hyq : (y : ℚ) ≤ 0 := by exact_mod_cast hyn
      have hm1 := mul_le_mul_of_nonpos_right ha2 hyq
      have hm2 := mul_le_mul_of_nonpos_right ha1 hyq
      have hh1 := hbounds.1
      have hh2 := hbounds.2
      simp only [if_neg hy] at hh1 hh2
      constructor <;> nlinarith
  have hlo : (out.lower : ℚ) ≤ x * (y : ℚ) :=
    (mul_le_mul_right hs).mp (by nlinarith [hprod.1])
  have hhi : x * (y : ℚ) ≤ (out.upper : ℚ) :=
    (mul_le_mul_right hs).mp (by nlinarith [hprod.2])
  simpa only [contains, hid] using And.intro hlo hhi

end FoldBounds.MotionRoundedOps
