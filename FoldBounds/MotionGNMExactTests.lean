module

public import FoldBounds.MotionGNMExact

@[expose] public section
namespace FoldBounds.MotionGNMExactTests

open FoldBounds.MotionUncertainty
open FoldBounds.MotionContactCompleteness
open FoldBounds.MotionGNMExact

private def lineCoords : Fin 3 → Fin 3 → ℚ :=
  fun i k => if k = 0 then (i.val : ℚ) else 0

private def lineEdges : Fin 2 → Spring 3 :=
  fun e => if e = 0 then canonicalPair 0 1 else canonicalPair 1 2

private def lineInverse : Matrix (Fin 2) (Fin 2) ℚ :=
  fun i j => if i = 0 then (if j = 0 then 2 else 1) else 1

-- Two unit contacts, 0--1--2.  Ground vertex 2; resistance(0,2)=2.
example : checkExactGNM (2 : Fin 3) lineEdges lineCoords
    1 (canonicalPair 0 2) lineInverse = some 2 := by
  decide

-- Omitted last contact must be rejected even with a plausible inverse.
private def incompleteEdges : Fin 1 → Spring 3 := fun _ => canonicalPair 0 1

example : checkExactGNM (2 : Fin 3) incompleteEdges lineCoords
    1 (canonicalPair 0 2) lineInverse = none := by
  decide

-- Wrong inverse must be rejected even on a complete contact network.
private def wrongInverse : Matrix (Fin 2) (Fin 2) ℚ :=
  fun i j => if i = j then 1 else 0

example : checkExactGNM (2 : Fin 3) lineEdges lineCoords
    1 (canonicalPair 0 2) wrongInverse = none := by
  decide

end FoldBounds.MotionGNMExactTests