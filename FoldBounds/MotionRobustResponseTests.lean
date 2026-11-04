import FoldBounds.MotionRobustResponse

open FoldBounds.MotionRobustResponse
open FoldBounds.MotionContactCompleteness
open FoldBounds.MotionUncertainty
open FoldBounds.MotionHessian
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse

private def pair : Spring 2 := ⟨0, 1, ![1, 0, 0]⟩
private def x8 : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![8, 0, 0]]
private def x12 : Fin 2 → Fin 3 → ℚ :=
  ![![0, 0, 0], ![12, 0, 0]]
private def v : Fin 2 × Fin 3 → ℚ := fun _ => 0

/-- A structural change crosses the 10A contact cutoff, literally changing
the number of springs from 1 to 0. -/
private def worlds : Bool → ProteinWorld 2
  | false => ⟨1, ![pair], x8⟩
  | true => ⟨0, ![], x12⟩

private theorem uniform_two_world_envelope :
    UniformEnvelope worlds 100 (1 / 5) pair v
      64 144 0 0 288 := by
  refine ⟨by norm_num, by norm_num, ?_, ?_, ?_, ?_⟩
  · intro ω
    cases ω <;> native_decide
  · intro ω
    cases ω <;> native_decide
  · intro ω
    cases ω <;> native_decide
  · intro ω
    cases ω <;> native_decide

/-- The exact same rational score bounds hold on both sides of a real
contact-network switch: the theorem cannot assume a fixed adjacency graph. -/
example :
    ∀ ω : Bool, ∃! y : Fin 2 × Fin 3 → ℚ,
      matVecQ
        (proteinHessian (worlds ω).edges (worlds ω).coords 100 (1 / 5)) y =
        pairForce pair (worlds ω).coords ∧
      (0 : ℚ) ≤ worldScore (worlds ω) pair y ∧
      worldScore (worlds ω) pair y ≤ 45 / 2 := by
  have h := uniform_score_sound worlds 100 (1 / 5)
    pair v 64 144 0 0 288 uniform_two_world_envelope
  have hl : (max (0 : ℚ) 0) / 144 = 0 := by norm_num
  have hu : ((0 : ℚ) + 288 / (1 / 5)) / 64 = 45 / 2 := by norm_num
  intro ω
  simpa only [hl, hu] using (h ω)
