module

public import FoldBounds.MotionGNMExact
public import FoldBounds.MotionExactSolvability

/-!
# Source-only exact electrical-flow certificates for grounded GNM

An arbitrary untrusted `P : vertex → contact → ℚ` is checked by the finite
rational identity `P * incidence = I` (grounded vertices). This identity,
and **not** an external rooted-tree path lemma, supplies the complete
Poincaré/coercivity premise needed by the exact GNM residual theorem.

By finite Cauchy–Schwarz, for every grounded displacement vector u,
    ||u||² ≤ (∑_v∑_e P[v,e]²) * ∑_e (incidence[e] · u)².

Any supplied flow is allowed. Numerical optimization, dyadic rounding and
path-based residual repairs are UNTRUSTED proposals. All input rationals are
checked by `checkFlowCertificate`. A proven solution exists uniquely, with no
positive ridge term and no unproved coercivity/solvability assumption.

A single certificate on *certain-on* edges controls every supergraph whose
additional contact spring coefficients are nonnegative; the separate Rayleigh
monotonicity theorem supplies the graph-inclusion implication.

The theorem covers in-memory rational data. It does NOT claim that the Lean
kernel decodes original NPZ/ZIP bytes or independently audits deposited NMR
structures; a separate integer byte verifier performs those checks.
-/

@[expose] public section
namespace FoldBounds.MotionGNMFlowCoercivity

open FoldBounds.MotionHessian
open FoldBounds.MotionResidual
open FoldBounds.MotionGNMExact
open FoldBounds.MotionExactSolvability

/-- All candidate grounded vertex basis vectors are reconstructed exactly
from untrusted oriented edge flows. -/
def checkFlowIdentity {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ) : Bool :=
  decide (∀ i j : Fin n,
    (∑ e : Fin m, flow i e * incidence e j) =
      if i = j then (1 : ℚ) else 0)

/-- Exact rational coercivity constant derived from the proposed flow. -/
def flowBound {n m : ℕ} (flow : Fin n → Fin m → ℚ) : ℚ :=
  ∑ i : Fin n, ∑ e : Fin m, (flow i e) ^ 2

/-- A complete executable source-only flow gate. -/
def checkFlowCertificate {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ) : Bool :=
  decide (0 < flowBound flow) && checkFlowIdentity incidence flow

theorem checkFlowIdentity_sound {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ)
    (h : checkFlowIdentity incidence flow = true) :
    ∀ i j : Fin n,
      (∑ e : Fin m, flow i e * incidence e j) =
        if i = j then (1 : ℚ) else 0 := by
  exact of_decide_eq_true h

theorem checkFlowCertificate_sound {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ)
    (h : checkFlowCertificate incidence flow = true) :
    0 < flowBound flow ∧
      (∀ i j : Fin n,
        (∑ e : Fin m, flow i e * incidence e j) =
          if i = j then (1 : ℚ) else 0) := by
  have hh : decide (0 < flowBound flow) = true ∧
      checkFlowIdentity incidence flow = true := by
    simpa only [checkFlowCertificate, Bool.and_eq_true] using h
  exact ⟨of_decide_eq_true hh.1,
    checkFlowIdentity_sound incidence flow hh.2⟩

/-- Exact flow conservation reconstructs every grounded coordinate. -/
theorem flow_reconstruct {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ)
    (h : ∀ i j : Fin n,
      (∑ e : Fin m, flow i e * incidence e j) =
        if i = j then (1 : ℚ) else 0)
    (u : Fin n → ℚ) (i : Fin n) :
    u i = ∑ e : Fin m, flow i e * dotQ (incidence e) u := by
  calc
    u i = ∑ j : Fin n,
        (if i = j then (1 : ℚ) else 0) * u j := by simp
    _ = ∑ j : Fin n,
        (∑ e : Fin m, flow i e * incidence e j) * u j := by
          apply Finset.sum_congr rfl
          intro j _
          rw [h i j]
    _ = ∑ j : Fin n, ∑ e : Fin m,
        (flow i e * incidence e j) * u j := by
          apply Finset.sum_congr rfl
          intro j _
          rw [Finset.sum_mul]
    _ = ∑ e : Fin m, ∑ j : Fin n,
        (flow i e * incidence e j) * u j := by rw [Finset.sum_comm]
    _ = ∑ e : Fin m, flow i e * dotQ (incidence e) u := by
          apply Finset.sum_congr rfl
          intro e _
          unfold dotQ
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro j _
          ring

/-- The model-relative Cauchy–Schwarz coercivity bound is now discharged
by concrete finite equality checks, with no assumed graph connectivity. -/
theorem flow_coercivity {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ)
    (h : ∀ i j : Fin n,
      (∑ e : Fin m, flow i e * incidence e j) =
        if i = j then (1 : ℚ) else 0)
    (u : Fin n → ℚ) :
    dotQ u u ≤ flowBound flow *
      dotQ u (matVecQ (hessianQ incidence (fun _ => 1) 0) u) := by
  let energy : ℚ := ∑ e : Fin m, (dotQ (incidence e) u) ^ 2
  have hcs (i : Fin n) :
      (u i) ^ 2 ≤ (∑ e : Fin m, (flow i e) ^ 2) * energy := by
    rw [flow_reconstruct incidence flow h u i]
    exact Finset.sum_mul_sq_le_sq_mul_sq Finset.univ (flow i)
      (fun e => dotQ (incidence e) u)
  have hmatrix :
      dotQ u (matVecQ (hessianQ incidence (fun _ => 1) 0) u) =
        energy := by
    have he := dot_matVec_hessian incidence (fun _ => (1 : ℚ)) 0 u
    simpa [springEnergyQ, energy] using he
  have hnorm : dotQ u u = ∑ i : Fin n, (u i) ^ 2 := by
    unfold dotQ
    apply Finset.sum_congr rfl
    intro i _
    ring
  calc
    dotQ u u = ∑ i : Fin n, (u i) ^ 2 := hnorm
    _ ≤ ∑ i : Fin n,
        (∑ e : Fin m, (flow i e) ^ 2) * energy := by
          apply Finset.sum_le_sum
          intro i _
          exact hcs i
    _ = flowBound flow * energy := by
          simp only [flowBound, Finset.sum_mul]
    _ = flowBound flow *
        dotQ u (matVecQ (hessianQ incidence (fun _ => 1) 0) u) := by
          rw [hmatrix]

/-- Correctly constructed source flow closes the coercivity and existence
premises for a *unique* exact rational GNM displacement response. -/
theorem flow_response_exists_unique {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ)
    (z approx : Fin n → ℚ)
    (h : checkFlowCertificate incidence flow = true) :
    ∃! y : Fin n → ℚ,
      matVecQ (hessianQ incidence (fun _ => 1) 0) y = z ∧
      (2 * dotQ z approx -
          dotQ approx (matVecQ (hessianQ incidence (fun _ => 1) 0) approx)) ≤
        dotQ z y ∧
      dotQ z y ≤
        (2 * dotQ z approx -
          dotQ approx (matVecQ (hessianQ incidence (fun _ => 1) 0) approx)) +
        flowBound flow *
          dotQ (vsub z (matVecQ (hessianQ incidence (fun _ => 1) 0) approx))
            (vsub z (matVecQ (hessianQ incidence (fun _ => 1) 0) approx)) := by
  obtain ⟨hB,hidentity⟩ :=
    checkFlowCertificate_sound incidence flow h
  let H := hessianQ incidence (fun _ => (1 : ℚ)) 0
  have hcoercive : ∀ u : Fin n → ℚ,
      dotQ u u ≤ flowBound flow * dotQ u (matVecQ H u) := by
    intro u
    exact flow_coercivity incidence flow hidentity u
  have hpos : 0 < (1 / flowBound flow : ℚ) := one_div_pos.mpr hB
  have hsmall : ∀ u : Fin n → ℚ,
      (1 / flowBound flow : ℚ) * dotQ u u ≤ dotQ u (matVecQ H u) := by
    intro u
    have hdiv : dotQ u u / flowBound flow ≤ dotQ u (matVecQ H u) :=
      (div_le_iff₀ hB).2 (by simpa [mul_comm] using hcoercive u)
    simpa [one_div, div_eq_mul_inv, mul_comm] using hdiv
  obtain ⟨y,hy,hu⟩ :=
    coercive_matrix_solvable H (1 / flowBound flow) hpos hsmall z
  have hbounds := gnm_residual_enclosure H (flowBound flow) z approx y
    hB (hessianQ_symmetric incidence (fun _ => 1) 0)
    hcoercive (fun i => congrFun hy i)
  refine ⟨y, ⟨hy, hbounds.1, hbounds.2⟩, ?_⟩
  intro y' hy'
  exact hu y' hy'.1

/-- Readable executable form: the proposed flow and approximate response are
both untrusted, and a false exact conservation check returns none. -/
def checkFlowResponse {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ)
    (z approx : Fin n → ℚ) : Option (ℚ × ℚ) :=
  if checkFlowCertificate incidence flow then
    let H := hessianQ incidence (fun _ => (1 : ℚ)) 0
    let lower := 2 * dotQ z approx - dotQ approx (matVecQ H approx)
    let residual := vsub z (matVecQ H approx)
    some (lower, lower + flowBound flow * dotQ residual residual)
  else none

/-- Fully checked interval includes a uniquely existent rational response.
No untrusted inverse and no assumed solution or Poincaré inequality. -/
theorem checkFlowResponse_sound {n m : ℕ}
    (incidence : Fin m → Fin n → ℚ)
    (flow : Fin n → Fin m → ℚ)
    (z approx : Fin n → ℚ) (bounds : ℚ × ℚ)
    (h : checkFlowResponse incidence flow z approx = some bounds) :
    ∃! y : Fin n → ℚ,
      matVecQ (hessianQ incidence (fun _ => 1) 0) y = z ∧
      bounds.1 ≤ dotQ z y ∧ dotQ z y ≤ bounds.2 := by
  have hgate : checkFlowCertificate incidence flow = true := by
    cases hc : checkFlowCertificate incidence flow with
    | false => simp [checkFlowResponse, hc] at h
    | true => rfl
  have hb :
      (2 * dotQ z approx -
        dotQ approx (matVecQ (hessianQ incidence (fun _ => 1) 0) approx),
       (2 * dotQ z approx -
        dotQ approx (matVecQ (hessianQ incidence (fun _ => 1) 0) approx)) +
        flowBound flow *
          dotQ (vsub z (matVecQ (hessianQ incidence (fun _ => 1) 0) approx))
            (vsub z (matVecQ (hessianQ incidence (fun _ => 1) 0) approx))) =
      bounds := by
    simpa [checkFlowResponse, hgate] using h
  cases hb
  exact flow_response_exists_unique incidence flow z approx hgate

end FoldBounds.MotionGNMFlowCoercivity