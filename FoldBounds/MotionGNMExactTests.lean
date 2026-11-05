module

public import FoldBounds.MotionGNMExact

/-!
# GNM acceptance proof wiring tests

The original attempt to evaluate the complete `checkExactGNM` on a 3-node
example using `native_decide` failed to link the generated native symbol for
`canonicalPair` through a transitive public import under Lean 4.28. This is a
Lean compiler/codegen limitation, NOT a proof that the proposed test passed.

We retain compile-checked theorem integration tests below; the already
existing `MotionContactCompletenessTests` exercises graph completeness via
`native_decide`, and the independent Python tests cover successful and
adversarial GNM inverse/graph cases. The new compound checker is not yet
end-to-end evaluated in Lean on a concrete protein.
-/

@[expose] public section
namespace FoldBounds.MotionGNMExactTests

open FoldBounds.MotionUncertainty
open FoldBounds.MotionContactCompleteness
open FoldBounds.MotionGNMExact
open FoldBounds.MotionHessian

private def lineCoords : Fin 3 → Fin 3 → ℚ :=
  fun i k => if k = 0 then (i.val : ℚ) else 0

private def lineEdges : Fin 2 → Spring 3 :=
  fun e => if e = 0 then canonicalPair 0 1 else canonicalPair 1 2

private def lineInverse : Matrix (Fin 2) (Fin 2) ℚ :=
  fun i j => if i = 0 then (if j = 0 then 2 else 1) else 1

example (v : ℚ)
    (h : checkExactGNM (2 : Fin 3) lineEdges lineCoords
      1 (canonicalPair 0 2) lineInverse = some v) :
    ExactContactNetwork lineEdges lineCoords 1 :=
  (checkExactGNM_sound (2 : Fin 3) lineEdges lineCoords
    1 (canonicalPair 0 2) lineInverse v h).1

example (v : ℚ)
    (h : checkExactGNM (2 : Fin 3) lineEdges lineCoords
      1 (canonicalPair 0 2) lineInverse = some v) :
    ∃ y : Fin 2 → ℚ,
      matVecQ (groundedGNM (2 : Fin 3) lineEdges) y =
        gnmForce (2 : Fin 3) (canonicalPair 0 2) ∧
      dotQ (gnmForce (2 : Fin 3) (canonicalPair 0 2)) y = v :=
  (checkExactGNM_sound (2 : Fin 3) lineEdges lineCoords
    1 (canonicalPair 0 2) lineInverse v h).2

end FoldBounds.MotionGNMExactTests