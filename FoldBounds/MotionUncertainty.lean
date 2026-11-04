module

public import Mathlib

/-!
# FoldBounds motion-certificate and uncertain-contact mathematics

`residual_gap_bound` is the coercivity-based residual enclosure.
An arbitrary candidate contact supergraph supports optional on/off contacts
by assigning weight zero. Theorems hold for ALL allowable contact assignments,
not only sample network draws. Generated witnesses remain untrusted.
-/

@[expose] public section

namespace FoldBounds.MotionUncertainty

/-- Algebraic core of a matrix-residual certificate. A quadratic completion
gives gap E, coercivity supplies τ*N ≤ E, and the final premise follows from
the sum of squares of the residual minus τ times the displacement. -/
theorem residual_gap_bound (τ E N R : ℝ) (hτ : 0 < τ)
    (hN : 0 ≤ N) (hcoercive : τ * N ≤ E)
    (hsquare : 0 ≤ R - 2 * τ * E + τ ^ 2 * N) :
    0 ≤ E ∧ E ≤ R / τ := by
  have hEN : 0 ≤ E := le_trans (mul_nonneg hτ.le hN) hcoercive
  have hmul : τ * (τ * N) ≤ τ * E :=
    mul_le_mul_of_nonneg_left hcoercive hτ.le
  have hbound : τ * E ≤ R := by
    nlinarith
  constructor
  · exact hEN
  · exact (le_div_iff₀ hτ).2 (by nlinarith)

/-- Normalized response bound. For an exact model, L = 2zᵀx - xᵀAx,
E = zᵀA⁻¹z - L and R = ||z-Ax||², with normalization D > 0. -/
theorem normalized_residual_enclosure (τ D q L E N R : ℝ)
    (hτ : 0 < τ) (hD : 0 < D)
    (hcomplete : q * D - L = E)
    (hN : 0 ≤ N) (hcoercive : τ * N ≤ E)
    (hsquare : 0 ≤ R - 2 * τ * E + τ ^ 2 * N) :
    L / D ≤ q ∧ q ≤ (L + R / τ) / D := by
  obtain ⟨hE, hER⟩ := residual_gap_bound τ E N R hτ hN hcoercive hsquare
  constructor
  · apply (div_le_iff₀ hD).2
    nlinarith
  · apply (le_div_iff₀ hD).2
    nlinarith

/-- Sum of residual squares is nonnegative without numeric approximation. -/
theorem residual_square_nonneg {n : ℕ} (r d : Fin n → ℝ) (τ : ℝ) :
    0 ≤ ∑ i : Fin n, (r i - τ * d i) ^ 2 := by
  apply Finset.sum_nonneg
  intro i _
  positivity

/-- One rank-one elastic-network spring. Its weight is stored separately. -/
structure Spring (n : ℕ) where
  source : Fin n
  target : Fin n
  direction : Fin 3 → ℝ

/-- Extension of a spring under a displacement field. -/
def springStretch {n : ℕ} (e : Spring n)
    (u : Fin n → Fin 3 → ℝ) : ℝ :=
  ∑ k : Fin 3, e.direction k * (u e.source k - u e.target k)

/-- Quadratic form of an arbitrary finite contact graph. -/
def springEnergy {n m : ℕ} (edges : Fin m → Spring n)
    (w : Fin m → ℝ) (u : Fin n → Fin 3 → ℝ) : ℝ :=
  ∑ i : Fin m, w i * (springStretch (edges i) u) ^ 2

/-- Nonnegative contact stiffness implies a positive semidefinite energy. -/
theorem springEnergy_nonneg {n m : ℕ} (edges : Fin m → Spring n)
    (w : Fin m → ℝ) (u : Fin n → Fin 3 → ℝ)
    (hw : ∀ i, 0 ≤ w i) : 0 ≤ springEnergy edges w u := by
  unfold springEnergy
  apply Finset.sum_nonneg
  intro i _
  exact mul_nonneg (hw i) (sq_nonneg _)

/-- Loewner ordering on any fixed candidate supergraph. -/
theorem springEnergy_mono {n m : ℕ} (edges : Fin m → Spring n)
    (a b : Fin m → ℝ) (u : Fin n → Fin 3 → ℝ)
    (hab : ∀ i, a i ≤ b i) :
    springEnergy edges a u ≤ springEnergy edges b u := by
  unfold springEnergy
  apply Finset.sum_le_sum
  intro i _
  exact mul_le_mul_of_nonneg_right (hab i) (sq_nonneg _)

/-- Squared norm of a 3D displacement field. -/
def displacementSq {n : ℕ} (u : Fin n → Fin 3 → ℝ) : ℝ :=
  ∑ i : Fin n, ∑ k : Fin 3, (u i k) ^ 2

theorem displacementSq_nonneg {n : ℕ} (u : Fin n → Fin 3 → ℝ) :
    0 ≤ displacementSq u := by
  unfold displacementSq
  apply Finset.sum_nonneg
  intro i _
  apply Finset.sum_nonneg
  intro k _
  positivity

/-- Ridge stabilization. -/
def ridgeEnergy {n m : ℕ} (edges : Fin m → Spring n)
    (w : Fin m → ℝ) (τ : ℝ) (u : Fin n → Fin 3 → ℝ) : ℝ :=
  springEnergy edges w u + τ * displacementSq u

/-- Ridge coercivity is UNIFORM over all allowed nonnegative contact graphs. -/
theorem ridgeEnergy_coercive {n m : ℕ} (edges : Fin m → Spring n)
    (w : Fin m → ℝ) (τ : ℝ) (u : Fin n → Fin 3 → ℝ)
    (hw : ∀ i, 0 ≤ w i) :
    τ * displacementSq u ≤ ridgeEnergy edges w τ u := by
  dsimp [ridgeEnergy]
  have hs := springEnergy_nonneg edges w u hw
  linarith

/-- Uncertain network weights (missing edges have weight zero). -/
structure WeightEnvelope (m : ℕ) where
  lower : Fin m → ℝ
  upper : Fin m → ℝ
  lower_nonneg : ∀ i, 0 ≤ lower i
  lower_le_upper : ∀ i, lower i ≤ upper i

/-- Every network within the weight envelope obeys the same PSD sandwich. -/
theorem uncertain_contact_energy {n m : ℕ}
    (edges : Fin m → Spring n) (env : WeightEnvelope m)
    (actual : Fin m → ℝ) (u : Fin n → Fin 3 → ℝ)
    (hactual : ∀ i, env.lower i ≤ actual i ∧ actual i ≤ env.upper i) :
    springEnergy edges env.lower u ≤ springEnergy edges actual u ∧
      springEnergy edges actual u ≤ springEnergy edges env.upper u := by
  constructor
  · exact springEnergy_mono edges env.lower actual u (fun i => (hactual i).1)
  · exact springEnergy_mono edges actual env.upper u (fun i => (hactual i).2)

/-- The PSD sandwich is preserved by uniform ridge stabilization. -/
theorem uncertain_ridge_energy {n m : ℕ}
    (edges : Fin m → Spring n) (env : WeightEnvelope m)
    (actual : Fin m → ℝ) (τ : ℝ) (u : Fin n → Fin 3 → ℝ)
    (hactual : ∀ i, env.lower i ≤ actual i ∧ actual i ≤ env.upper i) :
    ridgeEnergy edges env.lower τ u ≤ ridgeEnergy edges actual τ u ∧
      ridgeEnergy edges actual τ u ≤ ridgeEnergy edges env.upper τ u := by
  obtain ⟨h1, h2⟩ := uncertain_contact_energy edges env actual u hactual
  constructor
  · dsimp [ridgeEnergy]
    linarith
  · dsimp [ridgeEnergy]
    linarith

/-- The quadratic response functional. -/
def responseFunctional {n : ℕ} (force : Fin n → Fin 3 → ℝ)
    (energy : (Fin n → Fin 3 → ℝ) → ℝ)
    (u : Fin n → Fin 3 → ℝ) : ℝ :=
  2 * (∑ i : Fin n, ∑ k : Fin 3, force i k * u i k) - energy u

/-- The optimum response is antitone in quadratic energy: larger spring
stiffness cannot increase the optimized compliance response. -/
theorem uncertain_response_antitone {n : ℕ}
    (force : Fin n → Fin 3 → ℝ)
    (eLower eActual eUpper : (Fin n → Fin 3 → ℝ) → ℝ)
    (qLower qActual qUpper : ℝ)
    (henergy : ∀ u, eLower u ≤ eActual u ∧ eActual u ≤ eUpper u)
    (hoptLower : ∀ u, responseFunctional force eLower u ≤ qLower)
    (hoptActual : ∀ u, responseFunctional force eActual u ≤ qActual)
    (hattainUpper : ∃ u, responseFunctional force eUpper u = qUpper)
    (hattainActual : ∃ u, responseFunctional force eActual u = qActual) :
    qUpper ≤ qActual ∧ qActual ≤ qLower := by
  obtain ⟨uUpper, huUpper⟩ := hattainUpper
  obtain ⟨uActual, huActual⟩ := hattainActual
  constructor
  · have he := (henergy uUpper).2
    have hf : responseFunctional force eUpper uUpper ≤
        responseFunctional force eActual uUpper := by
      dsimp [responseFunctional]
      linarith
    linarith [hoptActual uUpper]
  · have he := (henergy uActual).1
    have hf : responseFunctional force eActual uActual ≤
        responseFunctional force eLower uActual := by
      dsimp [responseFunctional]
      linarith
    linarith [hoptLower uActual]

/-- Combine verified endpoint responses into a bound valid for all interior
contact graphs without enumerating them. -/
theorem response_interval_from_endpoints (qLower qActual qUpper L U : ℝ)
    (horder : qUpper ≤ qActual ∧ qActual ≤ qLower)
    (hLower : L ≤ qUpper) (hUpper : qLower ≤ U) :
    L ≤ qActual ∧ qActual ≤ U := by
  constructor <;> linarith

end FoldBounds.MotionUncertainty
