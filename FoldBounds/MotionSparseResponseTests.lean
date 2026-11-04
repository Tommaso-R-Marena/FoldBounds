import FoldBounds.MotionSparseResponse

open FoldBounds.MotionSparseResponse
open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionHessian
open FoldBounds.MotionUncertainty

private def e : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def edges : Fin 1 → Spring 2 := ![e]
private def x : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
private def u : Fin 2 → Fin 3 → ℚ :=
  ![![-1280 / 277, 0, 0], ![1280 / 277, 0, 0]]
private def scoreBounds : ScaledInterval :=
  ⟨91527119150052664367993011218, 91527119150052664367993011219⟩
private def wrongBounds : ScaledInterval :=
  ⟨91527119150052664367993011219, 91527119150052664367993011219⟩

example : sparseMatVec edges x 100 (1 / 5) (flat u) =
    matVecQ (proteinHessian edges x 100 (1 / 5)) (flat u) := by
  exact sparseMatVec_eq edges x 100 (1 / 5) (flat u)

example : checkSparseScaledPair edges x 100 (1 / 5) e
    (flat u) scoreBounds = true := by native_decide
example : checkSparseScaledPair edges x 100 (1 / 5) e
    (flat u) wrongBounds = false := by native_decide
example : checkSparseScaledPair edges x 100 0 e
    (flat u) scoreBounds = false := by native_decide

example :
    ∃! y : (Fin 2 × Fin 3) → ℚ,
      matVecQ (proteinHessian edges x 100 (1 / 5)) y =
        pairForce e x ∧
      lowerQ scoreBounds ≤ pairScore e x y ∧
      pairScore e x y ≤ upperQ scoreBounds := by
  have hc : checkSparseScaledPair edges x 100 (1 / 5)
      e (flat u) scoreBounds = true := by native_decide
  exact checkSparseScaledPair_sound edges x 100 (1 / 5)
    e (flat u) scoreBounds hc
