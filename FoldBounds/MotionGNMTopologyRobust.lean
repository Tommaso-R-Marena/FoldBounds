module

public import FoldBounds.MotionGNMExact

/-!
# FoldBounds.MotionGNMTopologyRobust

The electrical-network Rayleigh principle for scalar GNM networks over a
**fixed candidate supergraph** with uncertain nonnegative spring weights.
Every edge may switch between 0 and 1 (or take any intermediate nonnegative
weight). For fixed residue-pair force z, adding contacts cannot increase
maximum quadratic compliance zᵀ L⁻¹ z, conditional on endpoint/actual
maximizers. This supports two-extreme-network enclosures over arbitrarily
many uncertain contacts without sampling every state.

Our Python geometry classifier derives the permitted weight box from
reference coordinate boxes; this file is the mathematical deduction assuming
that proven weight inclusion. The coordinate parser and 2,614-protein data
replay are separately tested outside Lean; no biological accuracy is claimed.
-/

@[expose] public section
namespace FoldBounds.MotionGNMTopologyRobust

open FoldBounds.MotionHessian

variable {I : Type*} [Fintype I] [DecidableEq I]

/-- Sum-of-squares scalar spring energy with arbitrary nonnegative weights. -/
def scalarNetworkEnergy {m : ℕ} (b : Fin m → I → ℚ)
    (w : Fin m → ℚ) (u : I → ℚ) : ℚ :=
  ∑ e : Fin m, w e * (dotQ (b e) u) ^ 2

/-- Contact insertion cannot decrease quadratic stiffness. -/
theorem scalarNetworkEnergy_mono {m : ℕ}
    (b : Fin m → I → ℚ) (lo hi : Fin m → ℚ)
    (u : I → ℚ) (hw : ∀ e, lo e ≤ hi e) :
    scalarNetworkEnergy b lo u ≤ scalarNetworkEnergy b hi u := by
  unfold scalarNetworkEnergy
  apply Finset.sum_le_sum
  intro e _
  exact mul_le_mul_of_nonneg_right (hw e) (sq_nonneg _)

/-- Variational response of a scalar contact network. -/
def scalarResponseDual {m : ℕ}
    (b : Fin m → I → ℚ) (w : Fin m → ℚ)
    (z u : I → ℚ) : ℚ :=
  2 * dotQ z u - scalarNetworkEnergy b w u

/-- An added-contact network has no larger optimal scalar response.
No independence of optional contacts is needed. -/
theorem scalarResponse_antitone {m : ℕ}
    (b : Fin m → I → ℚ) (lo hi : Fin m → ℚ)
    (z : I → ℚ) (qlo qhi : ℚ)
    (hweights : ∀ e, lo e ≤ hi e)
    (hmax : ∀ u, scalarResponseDual b lo z u ≤ qlo)
    (hattain : ∃ u, scalarResponseDual b hi z u = qhi) :
    qhi ≤ qlo := by
  obtain ⟨u, hu⟩ := hattain
  have he := scalarNetworkEnergy_mono b lo hi u hweights
  have hm := hmax u
  dsimp [scalarResponseDual] at hu hm
  linarith

/-- Uniform two-extreme-network enclosure: for EVERY graph between the
certain-on graph and the possible-on graph, the actual effective resistance
lies between the two endpoint optimum responses. -/
theorem scalarResponse_extremes {m : ℕ}
    (b : Fin m → I → ℚ)
    (wlo wa whi : Fin m → ℚ) (z : I → ℚ)
    (qlo qa qhi : ℚ)
    (hweight : ∀ e, wlo e ≤ wa e ∧ wa e ≤ whi e)
    (hloMax : ∀ u, scalarResponseDual b wlo z u ≤ qlo)
    (haMax : ∀ u, scalarResponseDual b wa z u ≤ qa)
    (haAttain : ∃ u, scalarResponseDual b wa z u = qa)
    (hhiAttain : ∃ u, scalarResponseDual b whi z u = qhi) :
    qhi ≤ qa ∧ qa ≤ qlo := by
  constructor
  · exact scalarResponse_antitone b wa whi z qa qhi
      (fun e => (hweight e).2) haMax hhiAttain
  · exact scalarResponse_antitone b wlo wa z qlo qa
      (fun e => (hweight e).1) hloMax haAttain

end FoldBounds.MotionGNMTopologyRobust