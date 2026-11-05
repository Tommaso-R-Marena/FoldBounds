module

public import FoldBounds.MotionContactCompleteness
public import FoldBounds.MotionVerifiedInverse

/-!
# FoldBounds.MotionGNMExact

A **scalar** graph-Laplacian model separate from the anisotropic 3D Hessian.
For vertices 0..n, remove one anchor and use the canonical incidence basis
on the remaining n coordinates. The resulting matrix is a grounded graph
Laplacian. This is the exact algebraic effective-resistance baseline of GNM.

An untrusted inverse matrix B is accepted only if a full exact-rational
matrix multiplication confirms A*B = I. An accepted certificate includes
both the *complete physical cutoff graph* and an exact solution / score.

The scalable Python verifier instead checks exact residual intervals with a
spanning-tree coercivity bound; its dataset-specific path argument is not
asserted kernel-verified by this module. No claim about NMR dynamics follows.
-/

@[expose] public section
namespace FoldBounds.MotionGNMExact

open FoldBounds.MotionUncertainty
open FoldBounds.MotionHessian
open FoldBounds.MotionVerifiedInverse
open FoldBounds.MotionContactCompleteness
open FoldBounds.MotionResidual

/-- Ground a full graph on n+1 residues by deleting one anchor coordinate. -/
def scalarIncidence {n m : ℕ}
    (anchor : Fin (n+1)) (edges : Fin m → Spring (n+1))
    (e : Fin m) (i : Fin n) : ℚ :=
  (if (edges e).source = anchor.succAbove i then 1 else 0) -
  (if (edges e).target = anchor.succAbove i then 1 else 0)

/-- Unweighted grounded Kirchhoff matrix, with NO ridge regularization. -/
def groundedGNM {n m : ℕ} (anchor : Fin (n+1))
    (edges : Fin m → Spring (n+1)) : Matrix (Fin n) (Fin n) ℚ :=
  hessianQ (scalarIncidence anchor edges) (fun _ => 1) 0

/-- Scalar differential force for a residue pair. -/
def gnmForce {n : ℕ} (anchor : Fin (n+1))
    (pair : Spring (n+1)) (i : Fin n) : ℚ :=
  (if pair.source = anchor.succAbove i then 1 else 0) -
  (if pair.target = anchor.succAbove i then 1 else 0)

/-- If a supplied exact rational inverse and the complete graph both check,
return the exact effective-resistance score. Every input is rational.
The inverse witness is not trusted; all its entries are checked. -/
def checkExactGNM {n m : ℕ}
    (anchor : Fin (n+1)) (edges : Fin m → Spring (n+1))
    (coords : Fin (n+1) → Fin 3 → ℚ)
    (cutoffSq : ℚ) (pair : Spring (n+1))
    (B : Matrix (Fin n) (Fin n) ℚ) : Option ℚ :=
  if checkContactNetwork edges coords cutoffSq then
    if checkRightInverse (groundedGNM anchor edges) B then
      some (dotQ (gnmForce anchor pair)
        (matVecQ B (gnmForce anchor pair)))
    else none
  else none

/-- Executable acceptance implies an exact solution to the *complete*
reference graph Laplacian and its certified response value. The proof uses
no floating-point eigensolver, custom axioms, or assumed inverse. -/
theorem checkExactGNM_sound {n m : ℕ}
    (anchor : Fin (n+1)) (edges : Fin m → Spring (n+1))
    (coords : Fin (n+1) → Fin 3 → ℚ)
    (cutoffSq : ℚ) (pair : Spring (n+1))
    (B : Matrix (Fin n) (Fin n) ℚ) (value : ℚ)
    (h : checkExactGNM anchor edges coords cutoffSq pair B = some value) :
    ExactContactNetwork edges coords cutoffSq ∧
      ∃ y : Fin n → ℚ,
        matVecQ (groundedGNM anchor edges) y = gnmForce anchor pair ∧
        dotQ (gnmForce anchor pair) y = value := by
  have hgraph : checkContactNetwork edges coords cutoffSq = true := by
    cases hc : checkContactNetwork edges coords cutoffSq with
    | false => simp [checkExactGNM, hc] at h
    | true => rfl
  have hInv : checkRightInverse (groundedGNM anchor edges) B = true := by
    cases hi : checkRightInverse (groundedGNM anchor edges) B with
    | false => simp [checkExactGNM, hgraph, hi] at h
    | true => rfl
  have hscore : dotQ (gnmForce anchor pair)
      (matVecQ B (gnmForce anchor pair)) = value := by
    simpa [checkExactGNM, hgraph, hInv] using h
  refine ⟨checkContactNetwork_sound edges coords cutoffSq hgraph,
    ⟨matVecQ B (gnmForce anchor pair), ?_, hscore⟩⟩
  exact checked_inverse_solves (groundedGNM anchor edges) B hInv
    (gnmForce anchor pair)

/-- General exact residual enclosure under a Poincare/grounded coercivity
certificate: ||u||² <= B*u'A*u for ALL grounded u. A spanning-tree path
argument can provide this premise; the data-level Python checker constructs
and validates that tree separately. -/
theorem gnm_residual_enclosure {I : Type*} [Fintype I] [DecidableEq I]
    (A : Matrix I I ℚ) (B : ℚ) (z x y : I → ℚ)
    (hB : 0 < B) (hA : ∀ i j, A i j = A j i)
    (hcoercive : ∀ u : I → ℚ,
      dotQ u u ≤ B * dotQ u (matVecQ A u))
    (hy : ∀ i, matVecQ A y i = z i) :
    (2 * dotQ z x - dotQ x (matVecQ A x)) ≤ dotQ z y ∧
      dotQ z y ≤ (2 * dotQ z x - dotQ x (matVecQ A x)) +
        B * dotQ (vsub z (matVecQ A x)) (vsub z (matVecQ A x)) := by
  have ht : 0 < (1 / B : ℚ) := one_div_pos.mpr hB
  have hc : ∀ u : I → ℚ,
      (1 / B : ℚ) * dotQ u u ≤ dotQ u (matVecQ A u) := by
    intro u
    have hdiv : dotQ u u / B ≤ dotQ u (matVecQ A u) :=
      (div_le_iff₀ hB).2 (by simpa [mul_comm] using hcoercive u)
    simpa [one_div, div_eq_mul_inv, mul_comm] using hdiv
  obtain ⟨hl, hu⟩ := symmetric_residual_enclosure A (1 / B) z x y
    ht hA hc hy
  constructor
  · exact hl
  · have hs :
      dotQ (vsub z (matVecQ A x)) (vsub z (matVecQ A x)) / (1 / B) =
      B * dotQ (vsub z (matVecQ A x)) (vsub z (matVecQ A x)) := by
        field_simp [ne_of_gt hB]
    rw [hs] at hu
    exact hu

end FoldBounds.MotionGNMExact