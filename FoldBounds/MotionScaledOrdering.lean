module

public import FoldBounds.MotionScaledCertificates

/-!
# FoldBounds.MotionScaledOrdering

Proof-carrying strict ordering of two fixed-structure protein motion scores
from independently checked 2^-96 integer interval certificates.

A rank comparison is accepted only when both rational model/remainder
certificates succeed and the outer integer endpoints are strictly separated.
Ties, overlapping intervals, invalid structures and bad witnesses abstain.

This is the exact acceptance semantics needed to verify the bulk audit's
1,715 adjacent comparisons, but those source files have NOT yet been
replayed through Lean.
-/

@[expose] public section
namespace FoldBounds.MotionScaledOrdering

open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionHessian
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionUncertainty

def checkScaledOrder {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pairA pairB : Spring n)
    (candidateA candidateB : (Fin n × Fin 3) → ℚ)
    (a b : ScaledInterval) : Bool :=
  (checkScaledPairScore edges coords cutoffSq τ pairA candidateA a &&
    checkScaledPairScore edges coords cutoffSq τ pairB candidateB b) &&
    decide (a.upper < b.lower)

/-- Successful integer interval comparison implies strict physical-model
score ordering for two uniquely solvable finite rational systems.
The returned yA and yB are exact solutions, not untrusted proposals. -/
theorem checkScaledOrder_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pairA pairB : Spring n)
    (candidateA candidateB : (Fin n × Fin 3) → ℚ)
    (a b : ScaledInterval)
    (h : checkScaledOrder edges coords cutoffSq τ
      pairA pairB candidateA candidateB a b = true) :
    ∃ yA yB : (Fin n × Fin 3) → ℚ,
      matVecQ (proteinHessian edges coords cutoffSq τ) yA =
        pairForce pairA coords ∧
      matVecQ (proteinHessian edges coords cutoffSq τ) yB =
        pairForce pairB coords ∧
      pairScore pairA coords yA < pairScore pairB coords yB := by
  have hc : ((checkScaledPairScore edges coords cutoffSq τ pairA candidateA a &&
      checkScaledPairScore edges coords cutoffSq τ pairB candidateB b) = true) ∧
      decide (a.upper < b.lower) = true := by
    exact Bool.and_eq_true.mp h
  have hab := Bool.and_eq_true.mp hc.1
  have ht : a.upper < b.lower := of_decide_eq_true hc.2
  have htq : (a.upper : ℚ) < (b.lower : ℚ) := by exact_mod_cast ht
  have hscaled : upperQ a < lowerQ b := by
    unfold upperQ lowerQ
    apply (div_lt_div_iff₀ binaryScale_pos binaryScale_pos).2
    exact (mul_lt_mul_right binaryScale_pos).2 htq
  obtain ⟨yA, hyA, _⟩ :=
    checkScaledPairScore_exists_unique_sound edges coords cutoffSq τ
      pairA candidateA a hab.1
  obtain ⟨yB, hyB, _⟩ :=
    checkScaledPairScore_exists_unique_sound edges coords cutoffSq τ
      pairB candidateB b hab.2
  refine ⟨yA, yB, hyA.1, hyB.1, ?_⟩
  exact lt_of_lt_of_le (lt_of_le_of_lt hyA.2.2 hscaled) hyB.2.1

end FoldBounds.MotionScaledOrdering
