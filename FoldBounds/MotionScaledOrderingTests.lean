import FoldBounds.MotionScaledOrdering

open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionScaledOrdering
open FoldBounds.MotionUncertainty
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionHessian

private def spring : Spring 3 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 3 := ![spring]
private def pairB : Spring 3 := ⟨1, 2, ![1, 0, 0]⟩
private def coords : Fin 3 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0], ![20, 0, 0]]
private def candidateA : Fin 3 → Fin 3 → ℚ :=
  ![![-1280 / 277, 0, 0], ![1280 / 277, 0, 0], ![0, 0, 0]]
private def candidateB : Fin 3 → Fin 3 → ℚ :=
  ![![-7350 / 277, 0, 0], ![-9270 / 277, 0, 0], ![60, 0, 0]]
private def a : ScaledInterval :=
  ⟨91527119150052664367993011218, 91527119150052664367993011219⟩
private def b : ScaledInterval :=
  ⟨617092998644495698043577880324, 617092998644495698043577880325⟩

-- Both intervals are exact outward rounds of actual model responses
-- 320/277 and 4315/554. A real 3D Hessian underlies each score.
example : pairScore spring coords (flat candidateA) = 320 / 277 :=
  by native_decide
example : pairScore pairB coords (flat candidateB) = 4315 / 554 :=
  by native_decide
example : checkScaledPairScore edges coords 100 (1 / 5)
    spring (flat candidateA) a = true := by native_decide
example : checkScaledPairScore edges coords 100 (1 / 5)
    pairB (flat candidateB) b = true := by native_decide
example : checkScaledOrder edges coords 100 (1 / 5)
    spring pairB (flat candidateA) (flat candidateB) a b = true :=
  by native_decide

-- The reverse comparison cannot pass; nor can a tied interval.
example : checkScaledOrder edges coords 100 (1 / 5)
    pairB spring (flat candidateB) (flat candidateA) b a = false :=
  by native_decide
example : checkScaledOrder edges coords 100 (1 / 5)
    spring spring (flat candidateA) (flat candidateA) a a = false :=
  by native_decide

-- This proof produces actual rational solution witnesses and a strict
-- model score ordering by combining two successful integer certificates.
example :
    ∃ yA yB : (Fin 3 × Fin 3) → ℚ,
      matVecQ (proteinHessian edges coords 100 (1 / 5)) yA =
        pairForce spring coords ∧
      matVecQ (proteinHessian edges coords 100 (1 / 5)) yB =
        pairForce pairB coords ∧
      pairScore spring coords yA < pairScore pairB coords yB := by
  have h : checkScaledOrder edges coords 100 (1 / 5) spring pairB
      (flat candidateA) (flat candidateB) a b = true := by native_decide
  exact checkScaledOrder_sound edges coords 100 (1 / 5)
    spring pairB (flat candidateA) (flat candidateB) a b h
