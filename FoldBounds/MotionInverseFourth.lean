module

public import FoldBounds.MotionVariableDirections

/-!
# Certified inverse-fourth-distance elastic contacts

This file closes the rational inverse-square-distance *squared* denominator
gap for the previously studied weighted ANM contact term. With
d² = ∑ₐ (xᵢₐ - xⱼₐ)² and s = ∑ₐ (xᵢₐ - xⱼₐ)(uᵢₐ - uⱼₐ),
the per-contact quadratic energy is 49*s²/(d²)² when d² ≤ cutoff²,
and zero otherwise.

All rational interval witnesses are computed from coordinate boxes in Lean.
The executable safety gate insists that every candidate edge has a strictly
positive squared-distance lower bound; a possibly colliding input must FAIL.
All results quantify over *every* rational coordinate assignment in its box.

The theorem certifies the energy functional for any fixed displacement field,
not yet the floating-point eigenvectors or inverse-Hessian response.
-/

@[expose] public section

namespace FoldBounds.MotionInverseFourth

open FoldBounds.MotionCoordinates
open FoldBounds.MotionVariableDirections
open FoldBounds.MotionUncertainty

/-- An interval-square lower endpoint is always nonnegative. -/
theorem squared_lower_nonneg (b : Box) :
    0 ≤ (squared b).lower := by
  unfold squared
  split_ifs
  · exact le_refl _
  · exact le_min (sq_nonneg _) (sq_nonneg _)

/-- Exact rational enclosure of a ratio of squares. -/
theorem squared_ratio_enclosure
    (nLo n nHi dLo d dHi : ℚ)
    (hnLo : 0 ≤ nLo) (hn1 : nLo ≤ n) (hn2 : n ≤ nHi)
    (hdLo : 0 < dLo) (hd1 : dLo ≤ d) (hd2 : d ≤ dHi) :
    nLo / dHi ^ 2 ≤ n / d ^ 2 ∧
      n / d ^ 2 ≤ nHi / dLo ^ 2 := by
  have hdpos : 0 < d := lt_of_lt_of_le hdLo hd1
  have hdhipos : 0 < dHi := lt_of_lt_of_le hdpos hd2
  have hn : 0 ≤ n := le_trans hnLo hn1
  have hnHi : 0 ≤ nHi := le_trans hn hn2
  have hdSqLo : dLo ^ 2 ≤ d ^ 2 := by
    have hp : 0 ≤ (d - dLo) * (d + dLo) :=
      mul_nonneg (by linarith) (by linarith)
    nlinarith
  have hdSqHi : d ^ 2 ≤ dHi ^ 2 := by
    have hp : 0 ≤ (dHi - d) * (dHi + d) :=
      mul_nonneg (by linarith) (by linarith)
    nlinarith
  constructor
  · apply (div_le_div_iff₀ (pow_pos hdhipos 2) (pow_pos hdpos 2)).2
    have h1 : nLo * d ^ 2 ≤ n * d ^ 2 :=
      mul_le_mul_of_nonneg_right hn1 (sq_nonneg _)
    have h2 : n * d ^ 2 ≤ n * dHi ^ 2 :=
      mul_le_mul_of_nonneg_left hdSqHi hn
    exact le_trans h1 h2
  · apply (div_le_div_iff₀ (pow_pos hdpos 2) (pow_pos hdLo 2)).2
    have h1 : n * dLo ^ 2 ≤ nHi * dLo ^ 2 :=
      mul_le_mul_of_nonneg_right hn2 (sq_nonneg _)
    have h2 : nHi * dLo ^ 2 ≤ nHi * d ^ 2 :=
      mul_le_mul_of_nonneg_left hdSqLo hnHi
    exact le_trans h1 h2

/-- Physical weighted-ANM spring energy: inverse fourth power in distance. -/
def weightedEnergyTerm {n : ℕ} (e : Spring n)
    (x u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) : ℚ :=
  if distSq e x ≤ cutoffSq then
    49 * ((stretch e x u) ^ 2 / (distSq e x) ^ 2)
  else 0

/-- Lower enclosure; only a contact that is definitely ON contributes. -/
def weightedLowerTerm {n : ℕ} (e : Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) : ℚ :=
  if (distInterval e boxes).upper ≤ cutoffSq then
    49 * ((squared (stretchInterval e boxes u)).lower /
      (distInterval e boxes).upper ^ 2)
  else 0

/-- Upper enclosure; possibly ON contacts use the positive minimum distance. -/
def weightedUpperTerm {n : ℕ} (e : Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) : ℚ :=
  if cutoffSq < (distInterval e boxes).lower then 0 else
    49 * ((squared (stretchInterval e boxes u)).upper /
      (distInterval e boxes).lower ^ 2)

/-- Exact executable safety gate, rejecting inconsistent boxes, nonpositive
cutoffs and any candidate edge whose distance may vanish. -/
def safeNetwork {n m : ℕ} (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box) (cutoffSq : ℚ) : Bool :=
  decide (0 < cutoffSq ∧
    (∀ j a, (boxes j a).lower ≤ (boxes j a).upper) ∧
    (∀ i, 0 < (distInterval (edges i) boxes).lower))

theorem safeNetwork_valid {n m : ℕ} (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box) (cutoffSq : ℚ)
    (h : safeNetwork edges boxes cutoffSq = true) :
    0 < cutoffSq ∧
    (∀ j a, (boxes j a).lower ≤ (boxes j a).upper) ∧
    (∀ i, 0 < (distInterval (edges i) boxes).lower) := by
  exact of_decide_eq_true h

/-- Soundness of the full variable-direction, inverse-fourth-distance
contact *for any rational coordinate configuration inside the boxes*. -/
theorem weightedEnergyTerm_sound {n : ℕ} (e : Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (x u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ)
    (hdmin : 0 < (distInterval e boxes).lower)
    (hx : ∀ j a, contains (boxes j a) (x j a)) :
    weightedLowerTerm e boxes u cutoffSq ≤
      weightedEnergyTerm e x u cutoffSq ∧
    weightedEnergyTerm e x u cutoffSq ≤
      weightedUpperTerm e boxes u cutoffSq := by
  have hd : contains (distInterval e boxes) (distSq e x) :=
    squaredDistance_sound (boxes e.source) (boxes e.target)
      (x e.source) (x e.target)
      (fun a => hx e.source a) (fun a => hx e.target a)
  have hs : contains (stretchInterval e boxes u) (stretch e x u) :=
    stretchBox_sound (boxes e.source) (boxes e.target)
      (fun a => u e.source a - u e.target a)
      (x e.source) (x e.target)
      (fun a => hx e.source a) (fun a => hx e.target a)
  have hsq := squared_sound (stretchInterval e boxes u) (stretch e x u) hs
  have hRatio := squared_ratio_enclosure
    (squared (stretchInterval e boxes u)).lower
    ((stretch e x u) ^ 2)
    (squared (stretchInterval e boxes u)).upper
    (distInterval e boxes).lower
    (distSq e x)
    (distInterval e boxes).upper
    (squared_lower_nonneg _) hsq.1 hsq.2 hdmin hd.1 hd.2
  have h0 : (0 : ℚ) ≤ 49 := by norm_num
  constructor
  · by_cases hAlways : (distInterval e boxes).upper ≤ cutoffSq
    · have hOn : distSq e x ≤ cutoffSq := le_trans hd.2 hAlways
      simpa [weightedLowerTerm, weightedEnergyTerm, hAlways, hOn] using
        (mul_le_mul_of_nonneg_left hRatio.1 h0)
    · simp only [weightedLowerTerm, if_neg hAlways]
      by_cases hOn : distSq e x ≤ cutoffSq
      · simp only [weightedEnergyTerm, if_pos hOn]
        exact mul_nonneg h0 (div_nonneg (sq_nonneg _) (sq_nonneg _))
      · simp [weightedEnergyTerm, hOn]
  · by_cases hOff : cutoffSq < (distInterval e boxes).lower
    · have hNotOn : ¬ distSq e x ≤ cutoffSq := by linarith [hd.1]
      simp [weightedUpperTerm, weightedEnergyTerm, hOff, hNotOn]
    · simp only [weightedUpperTerm, if_neg hOff]
      by_cases hOn : distSq e x ≤ cutoffSq
      · simpa [weightedEnergyTerm, hOn] using
          (mul_le_mul_of_nonneg_left hRatio.2 h0)
      · simp only [weightedEnergyTerm, if_neg hOn]
        have hUpperN : 0 ≤ (squared (stretchInterval e boxes u)).upper :=
          le_trans (squared_lower_nonneg _) (le_trans hsq.1 hsq.2)
        exact mul_nonneg h0 (div_nonneg hUpperN (sq_nonneg _))

/-- Complete finite-network enclosure for the weighted elastic ANM energy,
uniform over every rational realization and every fixed displacement u.
All division safety assumptions are checked by the executable gate. -/
theorem safeNetwork_energy_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box) (cutoffSq : ℚ)
    (u : Fin n → Fin 3 → ℚ)
    (h : safeNetwork edges boxes cutoffSq = true) :
    ∀ (x : Fin n → Fin 3 → ℚ),
      (∀ j a, contains (boxes j a) (x j a)) →
      (∑ i : Fin m, weightedLowerTerm (edges i) boxes u cutoffSq) ≤
      (∑ i : Fin m, weightedEnergyTerm (edges i) x u cutoffSq) ∧
      (∑ i : Fin m, weightedEnergyTerm (edges i) x u cutoffSq) ≤
      (∑ i : Fin m, weightedUpperTerm (edges i) boxes u cutoffSq) := by
  intro x hx
  have hmin := (safeNetwork_valid edges boxes cutoffSq h).2.2
  constructor
  · apply Finset.sum_le_sum
    intro i _
    exact (weightedEnergyTerm_sound (edges i) boxes x u cutoffSq
      (hmin i) hx).1
  · apply Finset.sum_le_sum
    intro i _
    exact (weightedEnergyTerm_sound (edges i) boxes x u cutoffSq
      (hmin i) hx).2


/-- Exact weighted network quadratic form for a fixed configuration. -/
def weightedNetworkEnergy {n m : ℕ} (edges : Fin m → Spring n)
    (x u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) : ℚ :=
  ∑ i : Fin m, weightedEnergyTerm (edges i) x u cutoffSq

/-- Every individual inverse-fourth spring contributes nonnegative energy. -/
theorem weightedEnergyTerm_nonneg {n : ℕ} (e : Spring n)
    (x u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) :
    0 ≤ weightedEnergyTerm e x u cutoffSq := by
  unfold weightedEnergyTerm
  split_ifs
  · positivity
  · exact le_refl _

/-- The physical weighted contact network is positive semidefinite in its
rational displacement field, including arbitrary allowed contact topologies. -/
theorem weightedNetworkEnergy_nonneg {n m : ℕ}
    (edges : Fin m → Spring n) (x u : Fin n → Fin 3 → ℚ)
    (cutoffSq : ℚ) :
    0 ≤ weightedNetworkEnergy edges x u cutoffSq := by
  unfold weightedNetworkEnergy
  apply Finset.sum_nonneg
  intro i _
  exact weightedEnergyTerm_nonneg (edges i) x u cutoffSq

/-- Rational 3D displacement-field squared norm. -/
def displacementSqQ {n : ℕ} (u : Fin n → Fin 3 → ℚ) : ℚ :=
  ∑ i : Fin n, ∑ a : Fin 3, (u i a) ^ 2

theorem displacementSqQ_nonneg {n : ℕ}
    (u : Fin n → Fin 3 → ℚ) : 0 ≤ displacementSqQ u := by
  unfold displacementSqQ
  apply Finset.sum_nonneg
  intro i _
  apply Finset.sum_nonneg
  intro a _
  exact sq_nonneg _

/-- Ridge regularization supplies coercivity for any nonnegative stiffness
network, including coordinate-dependent directions and switching contacts. -/
def weightedRidgeEnergy {n m : ℕ} (edges : Fin m → Spring n)
    (x u : Fin n → Fin 3 → ℚ) (cutoffSq τ : ℚ) : ℚ :=
  weightedNetworkEnergy edges x u cutoffSq + τ * displacementSqQ u

theorem weightedRidgeEnergy_coercive {n m : ℕ}
    (edges : Fin m → Spring n) (x u : Fin n → Fin 3 → ℚ)
    (cutoffSq τ : ℚ) :
    τ * displacementSqQ u ≤
      weightedRidgeEnergy edges x u cutoffSq τ := by
  dsimp [weightedRidgeEnergy]
  have h := weightedNetworkEnergy_nonneg edges x u cutoffSq
  linarith

/-- The active-network interval theorem stated directly in terms of the
physically weighted network energy. -/
theorem weightedNetworkEnergy_interval {n m : ℕ}
    (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box) (cutoffSq : ℚ)
    (u : Fin n → Fin 3 → ℚ)
    (h : safeNetwork edges boxes cutoffSq = true) :
    ∀ (x : Fin n → Fin 3 → ℚ),
      (∀ j a, contains (boxes j a) (x j a)) →
      (∑ i : Fin m, weightedLowerTerm (edges i) boxes u cutoffSq) ≤
        weightedNetworkEnergy edges x u cutoffSq ∧
      weightedNetworkEnergy edges x u cutoffSq ≤
        (∑ i : Fin m, weightedUpperTerm (edges i) boxes u cutoffSq) := by
  intro x hx
  exact safeNetwork_energy_sound edges boxes cutoffSq u h x hx


/-- Executable, fail-closed rational network-energy certificate checker.
Untrusted box descriptions are accepted only if all edges are collision-safe;
the successful result contains exact lower and upper energy endpoints. -/
def checkNetworkEnergy {n m : ℕ}
    (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ) :
    Option (ℚ × ℚ) :=
  if safeNetwork edges boxes cutoffSq then
    some (∑ i : Fin m, weightedLowerTerm (edges i) boxes u cutoffSq,
          ∑ i : Fin m, weightedUpperTerm (edges i) boxes u cutoffSq)
  else none

/-- **End-to-end executable soundness:** every successful return from the
actual rational checker encloses the weighted physical network energy for
all rational coordinate configurations inside the supplied boxes. The proof
does not assume that untrusted input boxes have passed a separate gate. -/
theorem checkNetworkEnergy_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (u : Fin n → Fin 3 → ℚ) (cutoffSq : ℚ)
    (bounds : ℚ × ℚ)
    (h : checkNetworkEnergy edges boxes u cutoffSq = some bounds) :
    ∀ (x : Fin n → Fin 3 → ℚ),
      (∀ j a, contains (boxes j a) (x j a)) →
      bounds.1 ≤ weightedNetworkEnergy edges x u cutoffSq ∧
      weightedNetworkEnergy edges x u cutoffSq ≤ bounds.2 := by
  have hs : safeNetwork edges boxes cutoffSq = true := by
    cases hval : safeNetwork edges boxes cutoffSq with
    | false => simp [checkNetworkEnergy, hval] at h
    | true => rfl
  have heq :
      (∑ i : Fin m, weightedLowerTerm (edges i) boxes u cutoffSq,
       ∑ i : Fin m, weightedUpperTerm (edges i) boxes u cutoffSq) = bounds := by
    simpa [checkNetworkEnergy, hs] using h
  cases heq
  exact weightedNetworkEnergy_interval edges boxes cutoffSq u hs

end FoldBounds.MotionInverseFourth
