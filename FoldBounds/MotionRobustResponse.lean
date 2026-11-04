module

public import FoldBounds.MotionContactCompleteness

/-!
# FoldBounds.MotionRobustResponse

Uniform exact-rational response bounds for arbitrary families of protein
structures in which BOTH coordinates and the complete 10-angstrom contact
network may change. Each world can have a DIFFERENT number of contacts.

The core certificate uses one untrusted approximate solution for every
possible world; validity is conditional on a uniform enclosure of its
variational lower expression L and squared residual R. These envelope
premises can be discharged by sound interval arithmetic over structural
coordinate boxes (not yet implemented as a Lean executable gate).

The mathematical theorem is uniform in the world family; no assumptions of
contact-set preservation, continuity or invertibility from floating-point
eigensolvers. Exact uniqueness follows from positive ridge coercivity.

This proves a model-computation guarantee. It makes NO claim about
experimental NMR dynamics or physical ensemble coverage.
-/

@[expose] public section
namespace FoldBounds.MotionRobustResponse

open FoldBounds.MotionUncertainty
open FoldBounds.MotionVariableDirections
open FoldBounds.MotionHessian
open FoldBounds.MotionResidual
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionExactSolvability
open FoldBounds.MotionContactCompleteness

/-- A world can supply a DIFFERENT number of contacts. This permits
contact insertion/removal as the underlying protein coordinates vary. -/
structure ProteinWorld (n : ℕ) where
  m : ℕ
  edges : Fin m → Spring n
  coords : Fin n → Fin 3 → ℚ

/-- The trial lower expression (unnormalized). -/
def trialLower {n : ℕ} (w : ProteinWorld n)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (v : Fin n × Fin 3 → ℚ) : ℚ :=
  let z := pairForce pair w.coords
  let hv := matVecQ (proteinHessian w.edges w.coords cutoffSq τ) v
  2 * dotQ z v - dotQ v hv

/-- Squared residual of the same model relative to untrusted witness v. -/
def trialResidualSq {n : ℕ} (w : ProteinWorld n)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (v : Fin n × Fin 3 → ℚ) : ℚ :=
  let z := pairForce pair w.coords
  let hv := matVecQ (proteinHessian w.edges w.coords cutoffSq τ) v
  dotQ (vsub z hv) (vsub z hv)

def pairSeparationSq {n : ℕ} (w : ProteinWorld n)
    (pair : Spring n) : ℚ :=
  distSq pair w.coords

/-- The normalized model score of a variable-contact world and exact solve. -/
def worldScore {n : ℕ} (w : ProteinWorld n)
    (pair : Spring n) (y : Fin n × Fin 3 → ℚ) : ℚ :=
  pairScore pair w.coords y

/-- A uniform structural-uncertainty certificate. Its five interval
premises apply to EVERY possible world, including contact switches.
The graph-completeness premise rejects adversarially omitted springs.
The structure may represent a finite ensemble or infinitely many rational
conformations inside a coordinate box. -/
structure UniformEnvelope {n : ℕ} {Ω : Type*}
    (worlds : Ω → ProteinWorld n) (cutoffSq τ : ℚ)
    (pair : Spring n) (v : Fin n × Fin 3 → ℚ)
    (Dlo Dhi Llo Lhi Rhi : ℚ) : Prop where
  ridge_pos : 0 < τ
  separation_pos : 0 < Dlo
  complete_graph :
    ∀ ω, checkContactNetwork (worlds ω).edges
      (worlds ω).coords cutoffSq = true
  separation_bound :
    ∀ ω, Dlo ≤ pairSeparationSq (worlds ω) pair ∧
      pairSeparationSq (worlds ω) pair ≤ Dhi
  trial_bound :
    ∀ ω, Llo ≤ trialLower (worlds ω) cutoffSq τ pair v ∧
      trialLower (worlds ω) cutoffSq τ pair v ≤ Lhi
  residual_bound :
    ∀ ω, trialResidualSq (worlds ω) cutoffSq τ pair v ≤ Rhi

/-- A uniform bound on a score of a complete, variable-contact model.
Unlike pointwise certificates, the SAME endpoints apply to every world
whose coordinates/contact network satisfy the supplied envelope.

The lower bound is clipped at zero using positive-definite Hessians.
A normalized quadratic response is never negative. -/
theorem uniform_score_sound {n : ℕ} {Ω : Type*}
    (worlds : Ω → ProteinWorld n)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (v : Fin n × Fin 3 → ℚ)
    (Dlo Dhi Llo Lhi Rhi : ℚ)
    (env : UniformEnvelope worlds cutoffSq τ pair v
      Dlo Dhi Llo Lhi Rhi) :
    ∀ ω : Ω, ∃! y : Fin n × Fin 3 → ℚ,
      matVecQ (proteinHessian (worlds ω).edges (worlds ω).coords
        cutoffSq τ) y = pairForce pair (worlds ω).coords ∧
      max 0 Llo / Dhi ≤ worldScore (worlds ω) pair y ∧
      worldScore (worlds ω) pair y ≤ (Lhi + Rhi / τ) / Dlo := by
  intro ω
  let w := worlds ω
  let b : Fin w.m → (Fin n × Fin 3) → ℚ :=
    fun e => incidence (w.edges e) w.coords
  let weights : Fin w.m → ℚ :=
    fun e => contactWeight (w.edges e) w.coords cutoffSq
  let z := pairForce pair w.coords
  let H := proteinHessian w.edges w.coords cutoffSq τ
  have hw : ∀ e : Fin w.m, 0 ≤ weights e := by
    intro e
    exact contactWeight_nonneg (w.edges e) w.coords cutoffSq
  obtain ⟨y, hy, hu⟩ :=
    rankOne_solvable b weights τ env.ridge_pos hw z
  have hh : H = hessianQ b weights τ := rfl
  have hres := rankOne_residual_enclosure b weights τ z v y
    env.ridge_pos hw (fun i => congrFun hy i)
  have ht : trialLower w cutoffSq τ pair v ≤ dotQ z y ∧
      dotQ z y ≤ trialLower w cutoffSq τ pair v +
        trialResidualSq w cutoffSq τ pair v / τ := by
    simpa [trialLower, trialResidualSq, H, hh, z, b, weights]
      using hres
  have hnumlo : Llo ≤ dotQ z y :=
    le_trans (env.trial_bound ω).1 ht.1
  have hnumhi : dotQ z y ≤ Lhi + Rhi / τ := by
    have hr := div_le_div_of_nonneg_right (env.residual_bound ω)
      env.ridge_pos.le
    linarith [(env.trial_bound ω).2, ht.2, hr]
  have hqnonneg : 0 ≤ dotQ z y := by
    have hcoer := hessianQ_coercive b weights τ y hw
    have hn := dotQ_self_nonneg y
    have heq : dotQ z y =
        dotQ y (matVecQ (hessianQ b weights τ) y) := by
      rw [dotQ_comm z y]
      congr 1
      exact hy.symm
    rw [heq]
    exact le_trans (mul_nonneg env.ridge_pos.le hn) hcoer
  have hDlo : 0 < pairSeparationSq w pair :=
    lt_of_lt_of_le env.separation_pos (env.separation_bound ω).1
  have hDhi : 0 < Dhi :=
    lt_of_lt_of_le env.separation_pos
      (le_trans (env.separation_bound ω).1
        (env.separation_bound ω).2)
  have hmaxnonneg : 0 ≤ max (0 : ℚ) Llo := le_max_left _ _
  have hmaxle : max (0 : ℚ) Llo ≤ dotQ z y :=
    max_le hqnonneg hnumlo
  have huppernonneg : 0 ≤ Lhi + Rhi / τ :=
    le_trans hqnonneg hnumhi
  have hlow :
      max (0 : ℚ) Llo / Dhi ≤ dotQ z y / pairSeparationSq w pair := by
    apply (div_le_div_iff₀ hDhi hDlo).2
    have h1 := mul_le_mul_of_nonneg_left
      (env.separation_bound ω).2 hmaxnonneg
    have h2 := mul_le_mul_of_nonneg_right hmaxle hDhi.le
    nlinarith
  have hupp :
      dotQ z y / pairSeparationSq w pair ≤
        (Lhi + Rhi / τ) / Dlo := by
    apply (div_le_div_iff₀ hDlo env.separation_pos).2
    have h1 := mul_le_mul_of_nonneg_right hnumhi env.separation_pos.le
    have h2 := mul_le_mul_of_nonneg_left
      (env.separation_bound ω).1 huppernonneg
    nlinarith
  refine ⟨y, ⟨hy, ?_, ?_⟩, ?_⟩
  · simpa [worldScore, pairScore, pairSeparationSq, w, z] using hlow
  · simpa [worldScore, pairScore, pairSeparationSq, w, z] using hupp
  · intro y' hy'
    exact hu y' hy'.1

end FoldBounds.MotionRobustResponse
