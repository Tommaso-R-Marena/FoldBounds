module

public import FoldBounds.MotionCoordinates

/-!
# Variable-direction elastic-contact enclosure

This file certifies the coordinate-dependent SPRING DIRECTIONS, not just
the on/off contact topology. For a fixed nonnegative rational stiffness w,
a coordinate-dependent rank-one contact contributes

  w * (sum_a (x_source,a - x_target,a) *
               (u_source,a - u_target,a))^2

when its squared distance is within cutoff; otherwise it contributes zero.
All endpoints are computed from rational coordinate intervals and all
operations are mathematically sound for every rational coordinate realization.
The physically used 49/d^4 *stretch^2 needs an additional division theorem.
-/

@[expose] public section

namespace FoldBounds.MotionVariableDirections

open FoldBounds.MotionCoordinates
open FoldBounds.MotionUncertainty

def addBox (a b : Box) : Box :=
  ⟨a.lower + b.lower, a.upper + b.upper⟩

theorem addBox_sound (a b : Box) (x y : ℚ)
    (hx : contains a x) (hy : contains b y) :
    contains (addBox a b) (x + y) := by
  dsimp [contains, addBox] at *
  constructor <;> linarith [hx.1, hx.2, hy.1, hy.2]

def scaleBox (a : Box) (c : ℚ) : Box :=
  if 0 ≤ c then ⟨c * a.lower, c * a.upper⟩
    else ⟨c * a.upper, c * a.lower⟩

theorem scaleBox_sound (a : Box) (c x : ℚ)
    (hx : contains a x) :
    contains (scaleBox a c) (c * x) := by
  dsimp [contains] at *
  by_cases hc : 0 ≤ c
  · simp only [scaleBox, if_pos hc]
    exact ⟨mul_le_mul_of_nonneg_left hx.1 hc,
           mul_le_mul_of_nonneg_left hx.2 hc⟩
  · have hcn : c ≤ 0 := le_of_lt (lt_of_not_ge hc)
    simp only [scaleBox, if_neg hc]
    exact ⟨mul_le_mul_of_nonpos_left hx.2 hcn,
           mul_le_mul_of_nonpos_left hx.1 hcn⟩

/-- Box enclosing the coordinate-dependent spring extension for fixed u. -/
def stretchBox (a b : Fin 3 → Box) (du : Fin 3 → ℚ) : Box :=
  { lower := ∑ k : Fin 3,
      (scaleBox (difference (a k) (b k)) (du k)).lower
    upper := ∑ k : Fin 3,
      (scaleBox (difference (a k) (b k)) (du k)).upper }

theorem stretchBox_sound (a b : Fin 3 → Box) (du x y : Fin 3 → ℚ)
    (hx : ∀ k, contains (a k) (x k))
    (hy : ∀ k, contains (b k) (y k)) :
    contains (stretchBox a b du)
      (∑ k : Fin 3, du k * (x k - y k)) := by
  constructor
  · change (∑ k : Fin 3,
      (scaleBox (difference (a k) (b k)) (du k)).lower) ≤
      ∑ k : Fin 3, du k * (x k - y k)
    apply Finset.sum_le_sum
    intro k _
    exact (scaleBox_sound _ _ _
      (difference_sound _ _ _ _ (hx k) (hy k))).1
  · change (∑ k : Fin 3, du k * (x k - y k)) ≤
      ∑ k : Fin 3,
        (scaleBox (difference (a k) (b k)) (du k)).upper
    apply Finset.sum_le_sum
    intro k _
    exact (scaleBox_sound _ _ _
      (difference_sound _ _ _ _ (hx k) (hy k))).2

/-- Raw rank-one stretch in exact rational coordinates. -/
def stretch {n : ℕ} (e : Spring n)
    (x u : Fin n → Fin 3 → ℚ) : ℚ :=
  ∑ a : Fin 3, (u e.source a - u e.target a) *
    (x e.source a - x e.target a)

def distSq {n : ℕ} (e : Spring n)
    (x : Fin n → Fin 3 → ℚ) : ℚ :=
  ∑ a : Fin 3, (x e.source a - x e.target a) ^ 2

def stretchInterval {n : ℕ} (e : Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (u : Fin n → Fin 3 → ℚ) : Box :=
  stretchBox (boxes e.source) (boxes e.target)
    (fun a => u e.source a - u e.target a)

def distInterval {n : ℕ} (e : Spring n)
    (boxes : Fin n → Fin 3 → Box) : Box :=
  squaredDistance (boxes e.source) (boxes e.target)

/-- Energy of a contact that may switch on/off. -/
def energyTerm {n : ℕ} (e : Spring n)
    (x u : Fin n → Fin 3 → ℚ)
    (cutoffSq w : ℚ) : ℚ :=
  if distSq e x ≤ cutoffSq then w * (stretch e x u) ^ 2 else 0

def lowerTerm {n : ℕ} (e : Spring n)
    (boxes : Fin n → Fin 3 → Box) (u : Fin n → Fin 3 → ℚ)
    (cutoffSq w : ℚ) : ℚ :=
  if (distInterval e boxes).upper ≤ cutoffSq then
    w * (squared (stretchInterval e boxes u)).lower else 0

def upperTerm {n : ℕ} (e : Spring n)
    (boxes : Fin n → Fin 3 → Box) (u : Fin n → Fin 3 → ℚ)
    (cutoffSq w : ℚ) : ℚ :=
  if cutoffSq < (distInterval e boxes).lower then 0 else
    w * (squared (stretchInterval e boxes u)).upper

/-- Fixed nonnegative stiffness with completely varying spring direction
and contact membership. -/
theorem energyTerm_sound {n : ℕ} (e : Spring n)
    (boxes : Fin n → Fin 3 → Box) (x u : Fin n → Fin 3 → ℚ)
    (cutoffSq w : ℚ) (hw : 0 ≤ w)
    (hx : ∀ j a, contains (boxes j a) (x j a)) :
    lowerTerm e boxes u cutoffSq w ≤ energyTerm e x u cutoffSq w ∧
      energyTerm e x u cutoffSq w ≤ upperTerm e boxes u cutoffSq w := by
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
  constructor
  · by_cases hAlways : (distInterval e boxes).upper ≤ cutoffSq
    · have hOn : distSq e x ≤ cutoffSq := le_trans hd.2 hAlways
      simpa [lowerTerm, energyTerm, hAlways, hOn] using
        (mul_le_mul_of_nonneg_left hsq.1 hw)
    · simp only [lowerTerm, if_neg hAlways]
      by_cases hOn : distSq e x ≤ cutoffSq
      · simp only [energyTerm, if_pos hOn]
        exact mul_nonneg hw (sq_nonneg _)
      · simp [energyTerm, hOn]
  · by_cases hOff : cutoffSq < (distInterval e boxes).lower
    · have hNotOn : ¬ distSq e x ≤ cutoffSq := by linarith [hd.1]
      simp [upperTerm, energyTerm, hOff, hNotOn]
    · simp only [upperTerm, if_neg hOff]
      by_cases hOn : distSq e x ≤ cutoffSq
      · simpa [energyTerm, hOn] using
          (mul_le_mul_of_nonneg_left hsq.2 hw)
      · simp only [energyTerm, if_neg hOn]
        exact mul_nonneg hw (le_trans (sq_nonneg _) hsq.2)

/-- Finite network soundness with position-dependent directions, even when
the set of active edges changes. No untrusted matrix inverse is involved. -/
theorem variable_direction_energy_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (x u : Fin n → Fin 3 → ℚ)
    (cutoffSq : ℚ) (w : Fin m → ℚ)
    (hw : ∀ i, 0 ≤ w i)
    (hx : ∀ j a, contains (boxes j a) (x j a)) :
    (∑ i : Fin m, lowerTerm (edges i) boxes u cutoffSq (w i)) ≤
      (∑ i : Fin m, energyTerm (edges i) x u cutoffSq (w i)) ∧
    (∑ i : Fin m, energyTerm (edges i) x u cutoffSq (w i)) ≤
      (∑ i : Fin m, upperTerm (edges i) boxes u cutoffSq (w i)) := by
  constructor
  · apply Finset.sum_le_sum
    intro i _
    exact (energyTerm_sound (edges i) boxes x u cutoffSq (w i)
      (hw i) hx).1
  · apply Finset.sum_le_sum
    intro i _
    exact (energyTerm_sound (edges i) boxes x u cutoffSq (w i)
      (hw i) hx).2

end FoldBounds.MotionVariableDirections
