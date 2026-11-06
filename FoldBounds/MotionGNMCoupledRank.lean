module

public import FoldBounds.MotionGNMIntegerInverse
public import FoldBounds.MotionGNMTopologyRobust

/-!
# Coupled GNM order certificates (two positive-semidefinite changes)

For any graph with A_lo <= A_world <= A_hi in Loewner order, the
inverse difference K=A_lo^-1-A_world^-1 and its complement
Q=A_world^-1-A_hi^-1 are both positive semidefinite.

The *same* K acts on both residue-pair forces. Hence the uncertainty in
R_j-R_i is correlated; bounding R_i and R_j independently loses proof power.
This module formalizes a finite rational, kernel-checkable scalar contraction
and the strict sign test. Matrices and protein file decoding are checked by
an independent integer verifier, not smuggled into the Lean trusted base.

This is a statement about a specified scalar GNM, not biological prediction.
-/

@[expose] public section
namespace FoldBounds.MotionGNMCoupledRank

/-- If a positive-leading-coefficient rational quadratic is nonnegative at
all rationals, its discriminant is nonpositive. -/
theorem quadratic_discriminant_bound
    (a b c : ℚ) (ha : 0 < a)
    (hpoly : ∀ t : ℚ, 0 ≤ a * t * t - 2 * c * t + b) :
    c * c ≤ a * b := by
  have hanz : a ≠ 0 := ne_of_gt ha
  have h := hpoly (c / a)
  have hid : a * (c / a) * (c / a) - 2 * c * (c / a) + b =
      (a * b - c * c) / a := by
    field_simp [hanz]
    ring
  rw [hid] at h
  have hmul := mul_nonneg h (le_of_lt ha)
  have hcancel : ((a * b - c * c) / a) * a = a * b - c * c := by
    field_simp [hanz]
  nlinarith [hmul]

/-- Two PSD quadratic forms, sampled along (t*u-v) and (t*u+v), imply a
*cross-form* contraction. The two sources are K and M-K in the GNM use case.
This is the source of correlations absent from independent intervals. -/
theorem psd_pair_contraction
    (puu puv pvv quu quv qvv : ℚ)
    (ha : 0 < puu + quu)
    (hp : ∀ t : ℚ, 0 ≤ puu * t * t - 2 * puv * t + pvv)
    (hq : ∀ t : ℚ, 0 ≤ quu * t * t + 2 * quv * t + qvv) :
    (puv - quv) * (puv - quv) ≤
      (puu + quu) * (pvv + qvv) := by
  apply quadratic_discriminant_bound (puu + quu) (pvv + qvv)
    (puv - quv) ha
  intro t
  have h1 := hp t
  have h2 := hq t
  nlinarith

/-- A midpoint whose square is strictly larger than the worst-case
quadratic variation enforces a positive true difference. -/
theorem positive_of_coupled_square
    (s actual variance : ℚ)
    (hs : 0 < s)
    (hbound : (2 * actual - s) * (2 * actual - s) ≤ variance)
    (hmargin : variance < s * s) :
    0 < actual := by
  by_contra hnot
  have hle : actual ≤ 0 := le_of_not_gt hnot
  have hprod : 0 ≤ (-actual) * (s - actual) := by
    exact mul_nonneg (by linarith) (by linarith)
  nlinarith [hprod]

/-- Complete *coupled* rank certificate. `di,dj` are the two drop magnitudes
from the sparse to dense graph; `cross` is the exact transfer-current pairing
from the endpoint inverse difference; `s` is the sum of endpoint order gaps.

Hypotheses `hp,hq` are PSD claims about the two inverse-difference pieces.
The data-specific identity of those pieces is not assumed to be checked in
Lean: it must be checked by the source/graph/integer-matrix verifier. -/
theorem coupled_strict_rank_sound
    (puu puv pvv quu quv qvv : ℚ)
    (di dj cross s actual : ℚ)
    (ha : 0 < puu + quu)
    (hp : ∀ t : ℚ, 0 ≤ puu * t * t - 2 * puv * t + pvv)
    (hq : ∀ t : ℚ, 0 ≤ quu * t * t + 2 * quv * t + qvv)
    (hu : puu + quu = di + dj + 2 * cross)
    (hv : pvv + qvv = di + dj - 2 * cross)
    (hc : puv - quv = 2 * actual - s)
    (hs : 0 < s)
    (hcheck : (di + dj) * (di + dj) <
        s * s + 4 * cross * cross) :
    0 < actual := by
  have hcone := psd_pair_contraction puu puv pvv quu quv qvv
    ha hp hq
  have hvariance :
      (2 * actual - s) * (2 * actual - s) ≤
      (di + dj) * (di + dj) - 4 * cross * cross := by
    calc
      (2 * actual - s) * (2 * actual - s) =
          (puv - quv) * (puv - quv) := by rw [hc]
      _ ≤ (puu + quu) * (pvv + qvv) := hcone
      _ = (di + dj) * (di + dj) - 4 * cross * cross := by
        rw [hu, hv]
        ring
  have hmargin :
      (di + dj) * (di + dj) - 4 * cross * cross < s * s := by
    linarith
  exact positive_of_coupled_square s actual
    ((di + dj) * (di + dj) - 4 * cross * cross)
    hs hvariance hmargin

end FoldBounds.MotionGNMCoupledRank