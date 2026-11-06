module

public import FoldBounds.MotionRobustResponse

/-!
# FoldBounds.MotionRobustOrdering

Kernel-level deduction of guaranteed protein pair-motion ordering across
a possibly infinite FAMILY OF STRUCTURES and CHANGING CONTACT NETWORKS.

If each pair has a sound uniform response enclosure and the upper bound of
A is strictly below the lower bound of B, then EVERY allowed structure
has pair-motion score A < B. The two exact solutions may differ and are
proved to exist uniquely from positive ridge regularization.

This separates the trusted mathematical theorem from untrusted numerical
coordinate-box interval construction. It does not claim an empirically
validated physical flexibility signal.
-/

@[expose] public section
namespace FoldBounds.MotionRobustOrdering

open FoldBounds.MotionRobustResponse
open FoldBounds.MotionHessian
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionUncertainty

theorem uniform_strict_pair_order {n : ℕ} {Ω : Type*}
    (worlds : Ω → ProteinWorld n)
    (cutoffSq τ : ℚ)
    (pairA pairB : Spring n)
    (candidateA candidateB : Fin n × Fin 3 → ℚ)
    (DloA DhiA LloA LhiA RhiA : ℚ)
    (DloB DhiB LloB LhiB RhiB : ℚ)
    (envA : UniformEnvelope worlds cutoffSq τ pairA candidateA
      DloA DhiA LloA LhiA RhiA)
    (envB : UniformEnvelope worlds cutoffSq τ pairB candidateB
      DloB DhiB LloB LhiB RhiB)
    (hsep : (LhiA + RhiA / τ) / DloA <
      max (0 : ℚ) LloB / DhiB) :
    ∀ ω : Ω,
      ∃ yA yB : Fin n × Fin 3 → ℚ,
        matVecQ (proteinHessian (worlds ω).edges
          (worlds ω).coords cutoffSq τ) yA =
          pairForce pairA (worlds ω).coords ∧
        matVecQ (proteinHessian (worlds ω).edges
          (worlds ω).coords cutoffSq τ) yB =
          pairForce pairB (worlds ω).coords ∧
        worldScore (worlds ω) pairA yA <
          worldScore (worlds ω) pairB yB := by
  intro ω
  obtain ⟨yA, ha, _⟩ := uniform_score_sound worlds cutoffSq τ
    pairA candidateA DloA DhiA LloA LhiA RhiA envA ω
  obtain ⟨yB, hb, _⟩ := uniform_score_sound worlds cutoffSq τ
    pairB candidateB DloB DhiB LloB LhiB RhiB envB ω
  refine ⟨yA, yB, ha.1, hb.1, ?_⟩
  exact lt_of_lt_of_le (lt_of_le_of_lt ha.2.2 hsep) hb.2.1

end FoldBounds.MotionRobustOrdering
