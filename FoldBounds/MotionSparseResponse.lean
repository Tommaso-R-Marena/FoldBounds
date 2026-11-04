module

public import FoldBounds.MotionScaledCertificates

/-!
# FoldBounds.MotionSparseResponse

Evaluate a fixed-structure protein Hessian residual in O(number of candidate
contacts incident to a residue) rather than materializing a dense 3N x 3N
matrix. The sparse expression is PROVED equal to the kernel-checked explicit
Cartesian Hessian's action; an untrusted numerical candidate y is not trusted.

A source-to-target spring contributes its signed force to the two endpoint
residues, with all remaining residues receiving an exact zero contribution.
The expression handles inactive contact cutoffs and even self-edges by
matching the signed Cartesian incidence definition.

This is the intended reference checker for actual frozen CspB/crambin
coordinate inputs. It does not itself decode NPZ bytes.
-/

@[expose] public section

namespace FoldBounds.MotionSparseResponse

open FoldBounds.MotionHessian
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionResidual
open FoldBounds.MotionExactSolvability
open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionVariableDirections
open FoldBounds.MotionInverseFourth
open FoldBounds.MotionUncertainty

/-- Contact contribution in axis a for displacement candidate u. -/
def sparseEdgeTerm {n : ℕ} (e : Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ)
    (u : Fin n × Fin 3 → ℚ) (a : Fin 3) : ℚ :=
  contactWeight e coords cutoffSq * coordDelta e coords a *
    stretch e coords (fun j k => u (j,k))

/-- Sparse, incidence-aware spring contribution for one Cartesian row.
Self-edges contribute zero, just as in the explicit Hessian. -/
def sparseRowTerm {n : ℕ} (e : Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ)
    (u : Fin n × Fin 3 → ℚ) (p : Fin n × Fin 3) : ℚ :=
  if p.1 = e.source then
    if p.1 = e.target then 0
    else sparseEdgeTerm e coords cutoffSq u p.2
  else if p.1 = e.target then
    -sparseEdgeTerm e coords cutoffSq u p.2
  else 0

/-- Rowwise matrix-vector product, computing only incident contact forces. -/
def sparseMatVec {n m : ℕ} (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (u : Fin n × Fin 3 → ℚ)
    (p : Fin n × Fin 3) : ℚ :=
  (∑ e : Fin m, sparseRowTerm (edges e) coords cutoffSq u p) +
    τ * u p

theorem sparseRowTerm_eq {n : ℕ} (e : Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ)
    (u : Fin n × Fin 3 → ℚ) (p : Fin n × Fin 3) :
    sparseRowTerm e coords cutoffSq u p =
      (contactWeight e coords cutoffSq * incidence e coords p) *
        dotQ (incidence e coords) u := by
  have hdot : dotQ (incidence e coords) u =
      stretch e coords (fun j k => u (j,k)) :=
    incidence_dot e coords (fun j k => u (j,k))
  rw [hdot]
  by_cases hs : p.1 = e.source
  · by_cases ht : p.1 = e.target
    · simp [sparseRowTerm, incidence, hs, ht]
    · simp [sparseRowTerm, sparseEdgeTerm, incidence, hs, ht]
      ring
  · by_cases ht : p.1 = e.target
    · simp [sparseRowTerm, sparseEdgeTerm, incidence, hs, ht]
      ring
    · simp [sparseRowTerm, incidence, hs, ht]

/-- The sparse exact-rational computation equals the complete weighted
3N×3N Cartesian Hessian action for EVERY rational input and witness. -/
theorem sparseMatVec_eq {n m : ℕ} (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (u : Fin n × Fin 3 → ℚ) :
    sparseMatVec edges coords cutoffSq τ u =
      matVecQ (proteinHessian edges coords cutoffSq τ) u := by
  funext p
  have hh := matVecQ_hessian
    (fun e => incidence (edges e) coords)
    (fun e => contactWeight (edges e) coords cutoffSq)
    τ u p
  change matVecQ (proteinHessian edges coords cutoffSq τ) u p =
    (∑ e : Fin m,
      (contactWeight (edges e) coords cutoffSq *
        incidence (edges e) coords p) *
      dotQ (incidence (edges e) coords) u) + τ * u p at hh
  rw [hh]
  unfold sparseMatVec
  congr 1
  apply Finset.sum_congr rfl
  intro e _
  exact sparseRowTerm_eq (edges e) coords cutoffSq u p

/-- Fast exact-rational residual certificate endpoints; the mathematical
definition is deliberately independent of any floating-point matrix solver. -/
def sparseBounds {n m : ℕ} (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ)
    (pair : Spring n) (candidate : Fin n × Fin 3 → ℚ) :
    ℚ × ℚ :=
  let z := pairForce pair coords
  let hv := sparseMatVec edges coords cutoffSq τ candidate
  let r := vsub z hv
  let lo := 2 * dotQ z candidate - dotQ candidate hv
  (lo / distSq pair coords,
    (lo + dotQ r r / τ) / distSq pair coords)

/-- Fail-closed sparse score verifier. It shares all structural guards
with the already kernel-checked checker and recomputes the residual via
the proved sparse Hessian action. -/
def checkSparseScaledPair {n m : ℕ} (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ)
    (pair : Spring n) (candidate : Fin n × Fin 3 → ℚ)
    (proposal : ScaledInterval) : Bool :=
  if decide (0 < cutoffSq ∧ 0 < τ ∧
      0 < distSq pair coords ∧
      ∀ e : Fin m, 0 < distSq (edges e) coords) then
    let b := sparseBounds edges coords cutoffSq τ pair candidate
    checkEnclosure proposal b.1 b.2
  else false

/-- Sound sparse checker: a successful result guarantees existence and
uniqueness of the exact rational protein response, and encloses its
normalized pair score by the accepted 2^-96 interval. -/
theorem checkSparseScaledPair_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ)
    (pair : Spring n) (candidate : Fin n × Fin 3 → ℚ)
    (proposal : ScaledInterval)
    (h : checkSparseScaledPair edges coords cutoffSq τ
      pair candidate proposal = true) :
    ∃! y : Fin n × Fin 3 → ℚ,
      matVecQ (proteinHessian edges coords cutoffSq τ) y =
        pairForce pair coords ∧
      lowerQ proposal ≤ pairScore pair coords y ∧
      pairScore pair coords y ≤ upperQ proposal := by
  have hvalid : decide (0 < cutoffSq ∧ 0 < τ ∧
      0 < distSq pair coords ∧
      ∀ e : Fin m, 0 < distSq (edges e) coords) = true := by
    cases hv : decide (0 < cutoffSq ∧ 0 < τ ∧
        0 < distSq pair coords ∧
        ∀ e : Fin m, 0 < distSq (edges e) coords) with
    | false => simp [checkSparseScaledPair, hv] at h
    | true => rfl
  obtain ⟨hc, hτ, hp, he⟩ := of_decide_eq_true hvalid
  have hcheck : checkEnclosure proposal
      (sparseBounds edges coords cutoffSq τ pair candidate).1
      (sparseBounds edges coords cutoffSq τ pair candidate).2 = true := by
    simpa [checkSparseScaledPair, hvalid] using h
  have hin := checkEnclosure_sound proposal
    (sparseBounds edges coords cutoffSq τ pair candidate).1
    (sparseBounds edges coords cutoffSq τ pair candidate).2 hcheck
  have hw : ∀ e : Fin m, 0 ≤ contactWeight (edges e) coords cutoffSq := by
    intro e
    exact contactWeight_nonneg (edges e) coords cutoffSq
  let b : Fin m → (Fin n × Fin 3) → ℚ :=
    fun e => incidence (edges e) coords
  let w : Fin m → ℚ := fun e => contactWeight (edges e) coords cutoffSq
  obtain ⟨y, hy, huniq⟩ :=
    rankOne_solvable b w τ hτ hw (pairForce pair coords)
  have hh : sparseMatVec edges coords cutoffSq τ candidate =
      matVecQ (proteinHessian edges coords cutoffSq τ) candidate :=
    sparseMatVec_eq edges coords cutoffSq τ candidate
  have hresp := rankOne_residual_enclosure
    b w τ (pairForce pair coords) candidate y hτ hw
    (fun i => congrFun hy i)
  have hs : (sparseBounds edges coords cutoffSq τ pair candidate).1 ≤
      pairScore pair coords y ∧
      pairScore pair coords y ≤
        (sparseBounds edges coords cutoffSq τ pair candidate).2 := by
    dsimp [sparseBounds, pairScore]
    rw [hh]
    constructor
    · exact (div_le_div_iff_of_pos_right hp).2 hresp.1
    · exact (div_le_div_iff_of_pos_right hp).2 hresp.2
  refine ⟨y, ⟨hy, le_trans hin.1 hs.1,
    le_trans hs.2 hin.2⟩, ?_⟩
  intro y' hy'
  exact huniq y' hy'.1

end FoldBounds.MotionSparseResponse
