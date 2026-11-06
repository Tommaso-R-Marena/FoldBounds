module

public import FoldBounds.MotionResidual

/-!
# FoldBounds.MotionProteinResponse

Exact rational certificate for the **normalized pair-motion score** of
the previously implemented inverse-square-weighted protein ANM.

A physical pair with Cartesian difference v has force b_pair with
entries +v at the first residue and -v at the second. If H*y = b_pair,
then q = (b_pair.dot(y)) / (v.dot(v)). This exactly matches the spectral
unit-direction projection definition without introducing irrational
square roots.

The executable checker uses only the frozen reference structure,
the approximate solution, contact topology, spring coefficients and ridge.
It rejects invalid cutoff, collision-prone contacts, zero-length target
pair, nonpositive ridge, or negative spring weights.

The theorem is conditional on an exact solution of the finite rational
linear system. The computed numeric witness itself is untrusted.
-/

@[expose] public section

namespace FoldBounds.MotionProteinResponse

open FoldBounds.MotionUncertainty
open FoldBounds.MotionVariableDirections
open FoldBounds.MotionInverseFourth
open FoldBounds.MotionHessian
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionResidual

/-- Rational pair-force vector, no unit-length square root required. -/
def pairForce {n : ℕ} (pair : Spring n)
    (x : Fin n → Fin 3 → ℚ) : (Fin n × Fin 3) → ℚ :=
  incidence pair x

/-- Model-relative normalized pair score for an exact linear-system solution. -/
def pairScore {n : ℕ} (pair : Spring n)
    (x : Fin n → Fin 3 → ℚ)
    (y : (Fin n × Fin 3) → ℚ) : ℚ :=
  dotQ (pairForce pair x) y / distSq pair x

/-- Exact finite rational gate for a fixed structural reference, pair and
proposed approximate solve. No physical target NMR ensemble is accessed. -/
def checkPairScore {n m : ℕ}
    (edges : Fin m → Spring n)
    (x : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate : (Fin n × Fin 3) → ℚ) : Option (ℚ × ℚ) :=
  if decide (0 < cutoffSq ∧ 0 < distSq pair x ∧
      ∀ e : Fin m, 0 < distSq (edges e) x) then
    (checkResponse
      (fun e => incidence (edges e) x)
      (fun e => contactWeight (edges e) x cutoffSq)
      τ (pairForce pair x) candidate).map
        (fun bounds => (bounds.1 / distSq pair x,
                        bounds.2 / distSq pair x))
  else none

/-- Successful returned score bounds are sound for every exact solution of
the explicit, positive-ridge Cartesian protein Hessian. The checker
independently computes all upper/lower endpoints from rational inputs. -/
theorem checkPairScore_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (x : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate y : (Fin n × Fin 3) → ℚ)
    (bounds : ℚ × ℚ)
    (h : checkPairScore edges x cutoffSq τ pair candidate = some bounds)
    (hy : ∀ i, matVecQ (proteinHessian edges x cutoffSq τ) y i =
      pairForce pair x i) :
    bounds.1 ≤ pairScore pair x y ∧ pairScore pair x y ≤ bounds.2 := by
  have hgate : decide (0 < cutoffSq ∧ 0 < distSq pair x ∧
      ∀ e : Fin m, 0 < distSq (edges e) x) = true := by
    cases hg : decide (0 < cutoffSq ∧ 0 < distSq pair x ∧
        ∀ e : Fin m, 0 < distSq (edges e) x) with
    | false => simp [checkPairScore, hg] at h
    | true => rfl
  have hd : 0 < distSq pair x := (of_decide_eq_true hgate).2.1
  let b : Fin m → (Fin n × Fin 3) → ℚ :=
    fun e => incidence (edges e) x
  let w : Fin m → ℚ := fun e => contactWeight (edges e) x cutoffSq
  have hy' : ∀ i, matVecQ (hessianQ b w τ) y i =
      pairForce pair x i := by
    intro i
    exact hy i
  cases hc : checkResponse b w τ (pairForce pair x) candidate with
  | none =>
      simp [checkPairScore, hgate, b, w, hc] at h
  | some raw =>
      have hEq : (raw.1 / distSq pair x, raw.2 / distSq pair x) =
          bounds := by
        simpa [checkPairScore, hgate, b, w, hc] using h
      cases hEq
      have hs : raw.1 ≤ dotQ (pairForce pair x) y ∧
          dotQ (pairForce pair x) y ≤ raw.2 :=
        checkResponse_sound b w τ (pairForce pair x) candidate y raw hc hy'
      change raw.1 / distSq pair x ≤
          dotQ (pairForce pair x) y / distSq pair x ∧
        dotQ (pairForce pair x) y / distSq pair x ≤
          raw.2 / distSq pair x
      exact ⟨(div_le_div_iff_of_pos_right hd).2 hs.1,
             (div_le_div_iff_of_pos_right hd).2 hs.2⟩

end FoldBounds.MotionProteinResponse
