module

public import FoldBounds.MotionExactSolvability

/-!
# FoldBounds.MotionScaledCertificates

A fail-closed kernel-checkable bridge from external 2^96-scaled integer
score enclosures to the Lean-certified exact-rational protein response.

An untrusted bulk numerical proposer supplies only:
  - an approximate rational response witness;
  - two arbitrary signed integers (lo, hi), denoting lo/2^96 and hi/2^96.

The executable checker RECOMPUTES the underlying exact rational ANM
residual enclosure, then compares both proposed endpoints against it.
It returns false for invalid structural inputs, invalid ridge, incorrect
rounding, reversed intervals or a score interval narrower than justified
by the residual. Its proof does not depend on how the integer endpoints
or approximate solves were generated.

This provides a *slow but sound kernel-checkable reference gate*.
It does NOT yet certify the faster integer interval algorithm in
scripts/motion_fixedpoint_core.py, and it does not itself import the
1,717-case frozen NPZ archive.
-/

@[expose] public section

namespace FoldBounds.MotionScaledCertificates

open FoldBounds.MotionHessian
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionExactSolvability
open FoldBounds.MotionUncertainty

/-- Fixed binary denominator used by the independent outward-integer audit. -/
def binaryScale : ℚ := (2 : ℚ) ^ 96

theorem binaryScale_pos : 0 < binaryScale := by
  norm_num [binaryScale]

/-- A proposed closed rational interval in 2^-96 units. -/
structure ScaledInterval where
  lower : ℤ
  upper : ℤ
  deriving Repr, DecidableEq

def lowerQ (s : ScaledInterval) : ℚ := (s.lower : ℚ) / binaryScale
def upperQ (s : ScaledInterval) : ℚ := (s.upper : ℚ) / binaryScale

/-- Exact rational inclusion of a value, checked in integer-scaled units.
The checker does not trust the proposed interval constructor. -/
def checkValue (s : ScaledInterval) (v : ℚ) : Bool :=
  decide ((s.lower : ℚ) ≤ v * binaryScale ∧
    v * binaryScale ≤ (s.upper : ℚ))

theorem checkValue_sound (s : ScaledInterval) (v : ℚ)
    (h : checkValue s v = true) :
    lowerQ s ≤ v ∧ v ≤ upperQ s := by
  have hs : (s.lower : ℚ) ≤ v * binaryScale ∧
    v * binaryScale ≤ (s.upper : ℚ) := of_decide_eq_true h
  constructor
  · exact (div_le_iff₀ binaryScale_pos).2 hs.1
  · exact (le_div_iff₀ binaryScale_pos).2 hs.2

/-- Exact integer-scaled gates for an already-computed rational enclosure. -/
def checkEnclosure (s : ScaledInterval) (lo hi : ℚ) : Bool :=
  decide ((s.lower : ℚ) ≤ lo * binaryScale ∧
    hi * binaryScale ≤ (s.upper : ℚ))

theorem checkEnclosure_sound (s : ScaledInterval) (lo hi : ℚ)
    (h : checkEnclosure s lo hi = true) :
    lowerQ s ≤ lo ∧ hi ≤ upperQ s := by
  have hs : (s.lower : ℚ) ≤ lo * binaryScale ∧
    hi * binaryScale ≤ (s.upper : ℚ) := of_decide_eq_true h
  exact ⟨(div_le_iff₀ binaryScale_pos).2 hs.1,
    (le_div_iff₀ binaryScale_pos).2 hs.2⟩

/-- The proposal is rejected unless Lean's exact model-derived residual
enclosure succeeds *and* the supplied 2^96 integers outwardly enclose it. -/
def checkScaledPairScore {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate : (Fin n × Fin 3) → ℚ)
    (proposal : ScaledInterval) : Bool :=
  match checkPairScore edges coords cutoffSq τ pair candidate with
  | none => false
  | some (lo, hi) => checkEnclosure proposal lo hi

/-- A successful integer certificate guarantees that a UNIQUE exact
rational solution exists, and its normalized model score is inside
the proposed outward-rounded interval. Neither the approximate solver
nor the endpoint integers are trusted. -/
theorem checkScaledPairScore_exists_unique_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate : (Fin n × Fin 3) → ℚ)
    (proposal : ScaledInterval)
    (h : checkScaledPairScore edges coords cutoffSq τ
      pair candidate proposal = true) :
    ∃! y : (Fin n × Fin 3) → ℚ,
      matVecQ (proteinHessian edges coords cutoffSq τ) y =
        pairForce pair coords ∧
      lowerQ proposal ≤ pairScore pair coords y ∧
        pairScore pair coords y ≤ upperQ proposal := by
  cases hc : checkPairScore edges coords cutoffSq τ pair candidate with
  | none =>
      simp [checkScaledPairScore, hc] at h
  | some raw =>
      have hb : checkEnclosure proposal raw.1 raw.2 = true := by
        simpa [checkScaledPairScore, hc] using h
      have hbounds := checkEnclosure_sound proposal raw.1 raw.2 hb
      obtain ⟨y, hy, hu⟩ :=
        checkPairScore_exists_unique_sound edges coords cutoffSq τ
          pair candidate raw hc
      refine ⟨y, ⟨hy.1, ?_, ?_⟩, ?_⟩
      · exact le_trans hbounds.1 hy.2.1
      · exact le_trans hy.2.2 hbounds.2
      · intro y' hy'
        have hv := checkPairScore_sound edges coords cutoffSq τ pair
          candidate y' raw hc (fun i => congrFun hy'.1 i)
        exact hu y' ⟨hy'.1, hv.1, hv.2⟩

end FoldBounds.MotionScaledCertificates
