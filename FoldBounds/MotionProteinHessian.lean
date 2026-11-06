module

public import FoldBounds.MotionHessian

/-!
# Physical Hessian bridge for the inverse-fourth protein spring model

This module constructs the *explicit 3N by 3N matrix* (indexed by
residue x Cartesian axis) from the coordinate-dependent incidence
vectors. It proves that its exact quadratic form equals the existing
inverse-fourth-distance elastic energy plus ridge regularization.

The coefficient 49 / (squared distance)^2 matches MotionInverseFourth.
No matrix inverse, floating-point eigensolver, or untrusted oracle occurs
in the proof.
-/

@[expose] public section

namespace FoldBounds.MotionProteinHessian

open FoldBounds.MotionHessian
open FoldBounds.MotionUncertainty
open FoldBounds.MotionInverseFourth
open FoldBounds.MotionVariableDirections

def flat {n : ℕ} (u : Fin n → Fin 3 → ℚ) :
    (Fin n × Fin 3) → ℚ := fun p => u p.1 p.2

def coordDelta {n : ℕ} (e : Spring n)
    (x : Fin n → Fin 3 → ℚ) (a : Fin 3) : ℚ :=
  x e.source a - x e.target a

/-- For edge e and coordinate pair (i,a), the signed Cartesian gradient
coefficient of its extension. Self-edges cancel automatically. -/
def incidence {n : ℕ} (e : Spring n)
    (x : Fin n → Fin 3 → ℚ) (p : Fin n × Fin 3) : ℚ :=
  (if p.1 = e.source then coordDelta e x p.2 else 0) -
  (if p.1 = e.target then coordDelta e x p.2 else 0)

/-- Sum of one-hot contributions over residue indices. -/
private theorem sum_indicator_mul {n : ℕ}
    (k : Fin n) (c : ℚ) (v : Fin n → ℚ) :
    (∑ i : Fin n, (if i = k then c else 0) * v i) = c * v k := by
  calc
    (∑ i : Fin n, (if i = k then c else 0) * v i) =
        ∑ i : Fin n, if i = k then c * v i else 0 := by
          apply Finset.sum_congr rfl
          intro i _
          split_ifs <;> ring
    _ = c * v k := by simp

/-- This is the key Cartesian incidence identity: multiplying the
flattened displacement by the matrix incidence vector gives exactly
the contact extension from the original physical source code. -/
theorem incidence_dot {n : ℕ} (e : Spring n)
    (x u : Fin n → Fin 3 → ℚ) :
    dotQ (incidence e x) (flat u) = stretch e x u := by
  have hcomponent (a : Fin 3) :
      (∑ i : Fin n, incidence e x (i, a) * u i a) =
        (u e.source a - u e.target a) * coordDelta e x a := by
    change (∑ i : Fin n,
        ((if i = e.source then coordDelta e x a else 0) -
        (if i = e.target then coordDelta e x a else 0)) * u i a) =
        (u e.source a - u e.target a) * coordDelta e x a
    calc
      _ = (∑ i : Fin n, (if i = e.source then coordDelta e x a else 0) * u i a) -
            (∑ i : Fin n, (if i = e.target then coordDelta e x a else 0) * u i a) := by
              simp only [sub_mul, Finset.sum_sub_distrib]
      _ = coordDelta e x a * u e.source a -
            coordDelta e x a * u e.target a := by
              rw [sum_indicator_mul, sum_indicator_mul]
      _ = (u e.source a - u e.target a) * coordDelta e x a := by ring
  calc
    dotQ (incidence e x) (flat u) =
        ∑ p : Fin n × Fin 3, incidence e x p * u p.1 p.2 := rfl
    _ = ∑ i : Fin n, ∑ a : Fin 3, incidence e x (i, a) * u i a := by
          exact Fintype.sum_prod_type' (fun i a =>
            incidence e x (i, a) * u i a)
    _ = ∑ a : Fin 3, ∑ i : Fin n, incidence e x (i, a) * u i a := by
          rw [Finset.sum_comm]
    _ = ∑ a : Fin 3, (u e.source a - u e.target a) * coordDelta e x a := by
          apply Finset.sum_congr rfl
          intro a _
          exact hcomponent a
    _ = stretch e x u := rfl

def contactWeight {n : ℕ} (e : Spring n)
    (x : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) : ℚ :=
  if distSq e x ≤ cutoffSq then 49 / (distSq e x) ^ 2 else 0

theorem contactWeight_nonneg {n : ℕ} (e : Spring n)
    (x : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) :
    0 ≤ contactWeight e x cutoffSq := by
  unfold contactWeight
  split_ifs
  · positivity
  · exact le_refl _

theorem weight_stretch_energy {n : ℕ} (e : Spring n)
    (x u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) :
    contactWeight e x cutoffSq * (stretch e x u) ^ 2 =
      weightedEnergyTerm e x u cutoffSq := by
  by_cases h : distSq e x ≤ cutoffSq
  · simp only [contactWeight, weightedEnergyTerm, if_pos h]
    unfold stretch
    ring
  · simp [contactWeight, weightedEnergyTerm, h]

/-- Full 3N by 3N canonical protein Hessian with exact contact weights. -/
def proteinHessian {n m : ℕ} (edges : Fin m → Spring n)
    (x : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ) :
    Matrix (Fin n × Fin 3) (Fin n × Fin 3) ℚ :=
  hessianQ (fun e => incidence (edges e) x)
    (fun e => contactWeight (edges e) x cutoffSq) τ

/-- Flattened displacement squared norm equals the original 3D norm. -/
theorem flat_dot {n : ℕ} (u : Fin n → Fin 3 → ℚ) :
    dotQ (flat u) (flat u) = displacementSqQ u := by
  change (∑ p : Fin n × Fin 3, u p.1 p.2 * u p.1 p.2) =
    ∑ i : Fin n, ∑ a : Fin 3, (u i a) ^ 2
  calc
    (∑ p : Fin n × Fin 3, u p.1 p.2 * u p.1 p.2) =
        ∑ i : Fin n, ∑ a : Fin 3, u i a * u i a := by
          exact Fintype.sum_prod_type' (fun i a => u i a * u i a)
    _ = ∑ i : Fin n, ∑ a : Fin 3, (u i a) ^ 2 := by
          apply Finset.sum_congr rfl
          intro i _
          apply Finset.sum_congr rfl
          intro a _
          ring

/-- Main bridge: the actual physical energy is exactly the quadratic
form of the explicit Cartesian Hessian, not merely an analogy. -/
theorem proteinHessian_energy {n m : ℕ}
    (edges : Fin m → Spring n)
    (x u : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ) :
    dotQ (flat u) (matVecQ (proteinHessian edges x cutoffSq τ) (flat u)) =
      weightedRidgeEnergy edges x u cutoffSq τ := by
  change dotQ (flat u)
      (matVecQ (hessianQ (fun e => incidence (edges e) x)
        (fun e => contactWeight (edges e) x cutoffSq) τ) (flat u)) =
      weightedRidgeEnergy edges x u cutoffSq τ
  rw [dot_matVec_hessian]
  unfold springEnergyQ weightedRidgeEnergy weightedNetworkEnergy
  have hterms :
      (∑ e : Fin m, contactWeight (edges e) x cutoffSq *
        (dotQ (incidence (edges e) x) (flat u)) ^ 2) =
      (∑ e : Fin m, weightedEnergyTerm (edges e) x u cutoffSq) := by
    apply Finset.sum_congr rfl
    intro e _
    rw [incidence_dot]
    exact weight_stretch_energy (edges e) x u cutoffSq
  rw [hterms, flat_dot]

theorem proteinHessian_symm {n m : ℕ}
    (edges : Fin m → Spring n)
    (x : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ)
    (i j : Fin n × Fin 3) :
    proteinHessian edges x cutoffSq τ i j =
      proteinHessian edges x cutoffSq τ j i :=
  hessianQ_symmetric _ _ τ i j

/-- This is the coercivity premise needed for the matrix residual checker. -/
theorem proteinHessian_coercive {n m : ℕ}
    (edges : Fin m → Spring n)
    (x u : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ) :
    τ * dotQ (flat u) (flat u) ≤
      dotQ (flat u) (matVecQ (proteinHessian edges x cutoffSq τ) (flat u)) := by
  rw [proteinHessian_energy, flat_dot]
  exact weightedRidgeEnergy_coercive edges x u cutoffSq τ

end FoldBounds.MotionProteinHessian
