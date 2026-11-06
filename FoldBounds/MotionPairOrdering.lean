module

public import FoldBounds.MotionProteinResponse

/-!
# FoldBounds.MotionPairOrdering

A fail-closed comparison of two normalized motion-score certificates.
The checker returns true only when BOTH independent exact-rational
response enclosures succeed and the first upper endpoint lies strictly
below the second lower endpoint.

Numerical closeness, tied intervals, bad input, failed solves or
overlapping uncertainty lead to abstention. The soundness theorem is
about the executable Boolean result.
-/

@[expose] public section
namespace FoldBounds.MotionPairOrdering

open FoldBounds.MotionUncertainty
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionHessian

def strictlyOrdered (a b : ℚ × ℚ) : Bool := decide (a.2 < b.1)

theorem strictlyOrdered_sound (a b : ℚ × ℚ)
    (ha : strictlyOrdered a b = true)
    (x y : ℚ) (hx : x ≤ a.2) (hy : b.1 ≤ y) : x < y := by
  have ht : a.2 < b.1 := of_decide_eq_true ha
  exact lt_of_le_of_lt hx (lt_of_lt_of_le ht hy)

/-- Two independent rational score certificates, composed into a strict
ordering proof only if both succeed and their intervals are disjoint. -/
def checkPairOrder {n m : ℕ}
    (edges : Fin m → Spring n)
    (x : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ)
    (pairA pairB : Spring n)
    (candidateA candidateB : (Fin n × Fin 3) → ℚ) : Bool :=
  match checkPairScore edges x cutoffSq τ pairA candidateA,
        checkPairScore edges x cutoffSq τ pairB candidateB with
  | some a, some b => strictlyOrdered a b
  | _, _ => false

/-- Successful comparison implies an exact order of the two true finite
model responses. No floating-point sorting or post-hoc tie-breaking. -/
theorem checkPairOrder_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (x : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ)
    (pairA pairB : Spring n)
    (candidateA candidateB yA yB : (Fin n × Fin 3) → ℚ)
    (h : checkPairOrder edges x cutoffSq τ pairA pairB
      candidateA candidateB = true)
    (hyA : ∀ i, matVecQ (proteinHessian edges x cutoffSq τ) yA i =
      pairForce pairA x i)
    (hyB : ∀ i, matVecQ (proteinHessian edges x cutoffSq τ) yB i =
      pairForce pairB x i) :
    pairScore pairA x yA < pairScore pairB x yB := by
  cases hcA : checkPairScore edges x cutoffSq τ pairA candidateA with
  | none =>
      simp [checkPairOrder, hcA] at h
  | some boundsA =>
      cases hcB : checkPairScore edges x cutoffSq τ pairB candidateB with
      | none =>
          simp [checkPairOrder, hcA, hcB] at h
      | some boundsB =>
          have hstrict : strictlyOrdered boundsA boundsB = true := by
            simpa [checkPairOrder, hcA, hcB] using h
          have sa := checkPairScore_sound edges x cutoffSq τ pairA
            candidateA yA boundsA hcA hyA
          have sb := checkPairScore_sound edges x cutoffSq τ pairB
            candidateB yB boundsB hcB hyB
          exact strictlyOrdered_sound boundsA boundsB hstrict
            (pairScore pairA x yA) (pairScore pairB x yB) sa.2 sb.1

end FoldBounds.MotionPairOrdering
