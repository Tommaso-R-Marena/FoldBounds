module

public import FoldBounds.MotionSparseResponse

/-!
# FoldBounds.MotionContactCompleteness

Fail-closed exact-rational certification that a supplied contact enumeration is
exactly the full undirected, cutoff-limited contact graph for frozen coordinates.

Earlier score checkers prove the scores for the *supplied* edges only. Without
this additional gate, omitting a real contact (or duplicating one) would alter
the Hessian while leaving the algebraically checked score valid for the
incorrect graph. This module closes that model-identification trust boundary.

An accepted graph has exactly one correctly oriented edge i<j for every
reference-coordinate pair at/below the cutoff and no other edges. The
checker also rejects coincident residue coordinates, zero/negative cutoff,
duplicate edges, and incorrect orientation. Spring.direction does not matter
in the coordinate-derived Hessian; it is not part of contact identity.

The theorem is about rational in-memory coordinates and edge arrays. It is
NOT a byte-level certificate that an NPZ file was decoded correctly.
-/

@[expose] public section
namespace FoldBounds.MotionContactCompleteness

open FoldBounds.MotionUncertainty
open FoldBounds.MotionVariableDirections
open FoldBounds.MotionSparseResponse
open FoldBounds.MotionProteinHessian
open FoldBounds.MotionProteinResponse
open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionHessian

/-- A canonical contact pair, with direction intentionally unused by
the coordinate-derived matrix and pair-distance functions. -/
def canonicalPair {n : ℕ} (i j : Fin n) : Spring n :=
  ⟨i, j, fun _ => 0⟩

/-- Count the submitted oriented edges for a residue pair. This finite
cardinality avoids nonconstructive existential search in the executable gate. -/
def contactMultiplicity {n m : ℕ}
    (edges : Fin m → Spring n) (i j : Fin n) : ℕ :=
  (Finset.univ.filter (fun e : Fin m =>
    (edges e).source = i ∧ (edges e).target = j)).card

/-- Every pair below the cutoff has multiplicity exactly one, while
every submitted edge is active and oriented source<target. Together
these rules exclude missing, duplicate and spurious contacts. -/
def ExactContactNetwork {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq : ℚ) : Prop :=
  0 < cutoffSq ∧
  (∀ i j : Fin n, i < j →
    0 < distSq (canonicalPair i j) coords) ∧
  (∀ e : Fin m,
    (edges e).source < (edges e).target ∧
      distSq (edges e) coords ≤ cutoffSq) ∧
  (∀ i j : Fin n, i < j →
    (distSq (canonicalPair i j) coords ≤ cutoffSq ↔
      contactMultiplicity edges i j = 1))

/-- Executable complete-graph acceptance; all quantifiers range over
finite, explicit types and arithmetic uses exact rationals. -/
def checkContactNetwork {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq : ℚ) : Bool :=
  decide (ExactContactNetwork edges coords cutoffSq)

theorem checkContactNetwork_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ)
    (h : checkContactNetwork edges coords cutoffSq = true) :
    ExactContactNetwork edges coords cutoffSq := by
  exact of_decide_eq_true h

/-- In particular, no candidate network may exclude an active contact. -/
theorem checked_contact_present {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ)
    (h : checkContactNetwork edges coords cutoffSq = true)
    (i j : Fin n) (hij : i < j)
    (hc : distSq (canonicalPair i j) coords ≤ cutoffSq) :
    ∃ e : Fin m,
      (edges e).source = i ∧ (edges e).target = j := by
  have hcard : contactMultiplicity edges i j = 1 :=
    ((checkContactNetwork_sound edges coords cutoffSq h).2.2.2
      i j hij).mp hc
  have hpos : 0 <
      (Finset.univ.filter (fun e : Fin m =>
        (edges e).source = i ∧ (edges e).target = j)).card := by
    simpa [contactMultiplicity] using (show 0 < contactMultiplicity edges i j by
      rw [hcard]; norm_num)
  obtain ⟨e, he⟩ := Finset.card_pos.mp hpos
  have hm : (edges e).source = i ∧ (edges e).target = j := by
    simpa only [Finset.mem_filter, Finset.mem_univ, true_and] using he
  exact ⟨e, hm⟩

/-- Every supplied contact is canonical and under the cutoff. -/
theorem checked_contact_active {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ)
    (h : checkContactNetwork edges coords cutoffSq = true)
    (e : Fin m) :
    (edges e).source < (edges e).target ∧
      distSq (edges e) coords ≤ cutoffSq := by
  exact (checkContactNetwork_sound edges coords cutoffSq h).2.2.1 e

/-- A full-network protein certificate accepts ONLY when its graph is
complete and the sparse, exact-rational physical score certificate passes. -/
def checkCompleteSparseScore {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate : Fin n × Fin 3 → ℚ)
    (proposal : ScaledInterval) : Bool :=
  checkContactNetwork edges coords cutoffSq &&
    checkSparseScaledPair edges coords cutoffSq τ pair candidate proposal

/-- The result now certifies both the *correct network* and a unique
exact response of the corresponding complete reference protein Hessian. -/
theorem checkCompleteSparseScore_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) (pair : Spring n)
    (candidate : Fin n × Fin 3 → ℚ)
    (proposal : ScaledInterval)
    (h : checkCompleteSparseScore edges coords cutoffSq τ
      pair candidate proposal = true) :
    ExactContactNetwork edges coords cutoffSq ∧
      ∃! y : Fin n × Fin 3 → ℚ,
        matVecQ (proteinHessian edges coords cutoffSq τ) y =
          pairForce pair coords ∧
        lowerQ proposal ≤ pairScore pair coords y ∧
          pairScore pair coords y ≤ upperQ proposal := by
  have hgraph : checkContactNetwork edges coords cutoffSq = true := by
    cases hg : checkContactNetwork edges coords cutoffSq with
    | false => simp [checkCompleteSparseScore, hg] at h
    | true => rfl
  have hscore : checkSparseScaledPair edges coords cutoffSq τ
      pair candidate proposal = true := by
    simpa [checkCompleteSparseScore, hgraph] using h
  exact ⟨checkContactNetwork_sound edges coords cutoffSq hgraph,
    checkSparseScaledPair_sound edges coords cutoffSq τ pair
      candidate proposal hscore⟩

end FoldBounds.MotionContactCompleteness
