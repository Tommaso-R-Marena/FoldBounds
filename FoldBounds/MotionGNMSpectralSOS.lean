module

public import FoldBounds.MotionExactSolvability
public import FoldBounds.MotionGNMExact

/-!
# Proof-carrying, rational spectral sum-of-squares lower bounds

A source-only numerical proposal is allowed to suggest a rational spectral
lower bound τ and a family of positive-semidefinite square factors. The
*executable* rational checker verifies **all** matrix entries of

  A = τ I + Σ_e wₑ bₑ bₑᵀ,        τ > 0, wₑ ≥ 0.

The numerical proposal is entirely untrusted. Accepted equality closes the
coercivity, solvability and score interval premises for the unregularized GNM.
The existing rational Hessian theorems supply kernel-checked soundness.

The independent full-protein source verifier uses a compact integer
upper-triangular Cholesky proposal W plus a symmetric diagonally dominant
exact-integer residual S. It reconstructs S as an explicit nonnegative SOS:
  S = diag(d) + Σ_{i<j} |Sᵢⱼ|(eᵢ+sgn(Sᵢⱼ)eⱼ)²,
  dᵢ = Sᵢᵢ - Σ_{j≠i}|Sᵢⱼ| ≥ 0.
That decomposition can be supplied to this finite rational checker.

This module proves soundness of *in-memory* rational SOS certificates. It does
NOT assert that original NPZ/ZIP bytes or the 123x123 integer matrix are
parsed and checked inside the Lean kernel. A separate exact integer
verification program performs that data check.
-/

@[expose] public section
namespace FoldBounds.MotionGNMSpectralSOS

open FoldBounds.MotionHessian
open FoldBounds.MotionResidual
open FoldBounds.MotionExactSolvability

/-- Fail-closed exact spectral SOS gate for a symmetric grounded response
matrix, including all rational entries and nonnegative factor weights. -/
def checkSpectralSOS {n m : ℕ}
    (A : Matrix (Fin n) (Fin n) ℚ) (τ : ℚ)
    (vectors : Fin m → Fin n → ℚ) (weights : Fin m → ℚ) : Bool :=
  decide (0 < τ ∧ (∀ e : Fin m, 0 ≤ weights e) ∧
    (∀ i j : Fin n,
      A i j = hessianQ vectors weights τ i j))

theorem checkSpectralSOS_sound {n m : ℕ}
    (A : Matrix (Fin n) (Fin n) ℚ) (τ : ℚ)
    (vectors : Fin m → Fin n → ℚ) (weights : Fin m → ℚ)
    (h : checkSpectralSOS A τ vectors weights = true) :
    0 < τ ∧ (∀ e : Fin m, 0 ≤ weights e) ∧
      (∀ i j : Fin n,
        A i j = hessianQ vectors weights τ i j) := by
  exact of_decide_eq_true h

/-- An accepted untrusted sum-of-squares certificate establishes the entire
spectral/grounded Poincaré inequality: τ||u||² ≤ uᵀAu. -/
theorem checkSpectralSOS_coercive {n m : ℕ}
    (A : Matrix (Fin n) (Fin n) ℚ) (τ : ℚ)
    (vectors : Fin m → Fin n → ℚ) (weights : Fin m → ℚ)
    (h : checkSpectralSOS A τ vectors weights = true) :
    ∀ u : Fin n → ℚ,
      τ * dotQ u u ≤ dotQ u (matVecQ A u) := by
  obtain ⟨_,hw,hentries⟩ := checkSpectralSOS_sound A τ vectors weights h
  have hA : A = hessianQ vectors weights τ := by
    ext i j
    exact hentries i j
  intro u
  rw [hA]
  exact hessianQ_coercive vectors weights τ u hw

/-- In particular, accepted spectral SOS matrices have a unique exact
rational response for every force, without assuming existence. -/
theorem checkSpectralSOS_solvable {n m : ℕ}
    (A : Matrix (Fin n) (Fin n) ℚ) (τ : ℚ)
    (vectors : Fin m → Fin n → ℚ) (weights : Fin m → ℚ)
    (z : Fin n → ℚ)
    (h : checkSpectralSOS A τ vectors weights = true) :
    ∃! y : Fin n → ℚ, matVecQ A y = z := by
  have ht := (checkSpectralSOS_sound A τ vectors weights h).1
  exact coercive_matrix_solvable A τ ht
    (checkSpectralSOS_coercive A τ vectors weights h) z

/-- Executable residual enclosure from a *checked spectral lower bound*.
Neither the candidate response nor the SOS proposal is trusted. -/
def checkSpectralResponse {n m : ℕ}
    (A : Matrix (Fin n) (Fin n) ℚ) (τ : ℚ)
    (vectors : Fin m → Fin n → ℚ) (weights : Fin m → ℚ)
    (z approx : Fin n → ℚ) : Option (ℚ × ℚ) :=
  if checkSpectralSOS A τ vectors weights then
    let lower := 2 * dotQ z approx - dotQ approx (matVecQ A approx)
    let residual := vsub z (matVecQ A approx)
    some (lower, lower + dotQ residual residual / τ)
  else none

/-- An accepted spectral SOS gate entails a uniquely existing response and
its exact rational residual interval. This is a *model* guarantee, not a
claim about biological predictive accuracy. -/
theorem checkSpectralResponse_sound {n m : ℕ}
    (A : Matrix (Fin n) (Fin n) ℚ) (τ : ℚ)
    (vectors : Fin m → Fin n → ℚ) (weights : Fin m → ℚ)
    (z approx : Fin n → ℚ) (bounds : ℚ × ℚ)
    (h : checkSpectralResponse A τ vectors weights z approx = some bounds) :
    ∃! y : Fin n → ℚ,
      matVecQ A y = z ∧
      bounds.1 ≤ dotQ z y ∧ dotQ z y ≤ bounds.2 := by
  have hgate : checkSpectralSOS A τ vectors weights = true := by
    cases hg : checkSpectralSOS A τ vectors weights with
    | false => simp [checkSpectralResponse, hg] at h
    | true => rfl
  obtain ⟨hτ,hw,hentry⟩ := checkSpectralSOS_sound A τ vectors weights hgate
  have hmatrix : A = hessianQ vectors weights τ := by
    ext i j
    exact hentry i j
  have hbound :
      (2 * dotQ z approx - dotQ approx (matVecQ A approx),
       (2 * dotQ z approx - dotQ approx (matVecQ A approx)) +
         dotQ (vsub z (matVecQ A approx))
           (vsub z (matVecQ A approx)) / τ) = bounds := by
    simpa [checkSpectralResponse, hgate] using h
  obtain ⟨y,hy,hu⟩ := rankOne_solvable vectors weights τ hτ hw z
  have hres := rankOne_residual_enclosure vectors weights τ z approx y
      hτ hw (fun i => congrFun hy i)
  refine ⟨y, ?_, ?_⟩
  · rw [← hmatrix] at hy hres
    rw [← hbound]
    exact ⟨hy,hres.1,hres.2⟩
  · intro y' hy'
    exact hu y' (by simpa only [hmatrix] using hy'.1)

end FoldBounds.MotionGNMSpectralSOS