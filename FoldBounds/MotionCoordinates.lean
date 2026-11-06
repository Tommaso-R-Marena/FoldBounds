module

public import FoldBounds.MotionContactChecker

/-!
# FoldBounds.MotionCoordinates

Exact rational interval arithmetic for squared C-alpha distances, connected
directly to the finite executable contact-network gate. This proves that
axis-aligned coordinate boxes imply the supplied distance enclosures.

The result concerns **contact inclusion/exclusion and fixed spring directions**.
It does not yet bound the changes in spring *directions* caused by coordinate
perturbations; that is a separate model obligation.
-/

@[expose] public section

namespace FoldBounds.MotionCoordinates

open FoldBounds.MotionContactChecker
open FoldBounds.MotionUncertainty

/-- Small, transparent rational interval, with no trusted validity invariant. -/
structure Box where
  lower : ℚ
  upper : ℚ
  deriving DecidableEq, Repr

def contains (b : Box) (x : ℚ) : Prop :=
  b.lower ≤ x ∧ x ≤ b.upper

def difference (a b : Box) : Box :=
  ⟨a.lower - b.upper, a.upper - b.lower⟩

theorem difference_sound (a b : Box) (x y : ℚ)
    (hx : contains a x) (hy : contains b y) :
    contains (difference a b) (x - y) := by
  dsimp [contains, difference] at *
  constructor <;> linarith [hx.1, hx.2, hy.1, hy.2]

/-- Exact squared interval. If a component can cross zero its lower
squared bound is zero; otherwise it is the smaller endpoint square. -/
def squared (a : Box) : Box :=
  ⟨if a.lower ≤ 0 ∧ 0 ≤ a.upper then 0
      else min (a.lower ^ 2) (a.upper ^ 2),
   max (a.lower ^ 2) (a.upper ^ 2)⟩

theorem squared_sound (a : Box) (x : ℚ)
    (hx : contains a x) : contains (squared a) (x ^ 2) := by
  dsimp [contains, squared]
  constructor
  · by_cases hCross : a.lower ≤ 0 ∧ 0 ≤ a.upper
    · simp only [if_pos hCross]
      positivity
    · simp only [if_neg hCross]
      by_cases hLeft : 0 < a.lower
      · have hprod : 0 ≤ (x - a.lower) * (x + a.lower) :=
          mul_nonneg (by linarith [hx.1]) (by linarith [hx.1])
        have hsq : a.lower ^ 2 ≤ x ^ 2 := by nlinarith
        exact le_trans (min_le_left _ _) hsq
      · have hRight : a.upper < 0 := by
          by_contra hno
          have hle : a.lower ≤ 0 := by linarith
          have hge : 0 ≤ a.upper := by linarith
          exact hCross ⟨hle, hge⟩
        have hprod : 0 ≤ (a.upper - x) * (-a.upper - x) :=
          mul_nonneg (by linarith [hx.2]) (by linarith [hx.2, hRight])
        have hsq : a.upper ^ 2 ≤ x ^ 2 := by nlinarith
        exact le_trans (min_le_right _ _) hsq
  · by_cases hPos : 0 ≤ x
    · have hprod : 0 ≤ (a.upper - x) * (a.upper + x) :=
        mul_nonneg (by linarith [hx.2]) (by linarith [hx.2, hPos])
      have hsq : x ^ 2 ≤ a.upper ^ 2 := by nlinarith
      exact le_trans hsq (le_max_right _ _)
    · have hprod : 0 ≤ (x - a.lower) * (-x - a.lower) :=
        mul_nonneg (by linarith [hx.1]) (by linarith [hx.1, hPos])
      have hsq : x ^ 2 ≤ a.lower ^ 2 := by nlinarith
      exact le_trans hsq (le_max_left _ _)

/-- Pythagorean squared-distance interval from two 3D coordinate boxes. -/
def squaredDistance (a b : Fin 3 → Box) : Box :=
  { lower := ∑ k : Fin 3, (squared (difference (a k) (b k))).lower
    upper := ∑ k : Fin 3, (squared (difference (a k) (b k))).upper }

theorem squaredDistance_sound
    (a b : Fin 3 → Box) (x y : Fin 3 → ℚ)
    (hx : ∀ k, contains (a k) (x k))
    (hy : ∀ k, contains (b k) (y k)) :
    contains (squaredDistance a b)
      (∑ k : Fin 3, (x k - y k) ^ 2) := by
  constructor
  · change (∑ k : Fin 3, (squared (difference (a k) (b k))).lower) ≤
      ∑ k : Fin 3, (x k - y k) ^ 2
    apply Finset.sum_le_sum
    intro k _
    exact (squared_sound _ _ (difference_sound _ _ _ _ (hx k) (hy k))).1
  · change (∑ k : Fin 3, (x k - y k) ^ 2) ≤
      (∑ k : Fin 3, (squared (difference (a k) (b k))).upper)
    apply Finset.sum_le_sum
    intro k _
    exact (squared_sound _ _ (difference_sound _ _ _ _ (hx k) (hy k))).2

/-- Convert coordinate boxes for all residues into rational contact intervals.
The edge identifiers and stiffness intervals are untrusted but explicit. -/
def contactBoxes {n m : ℕ} (edges : Fin m → Spring n)
    (coords : Fin n → Fin 3 → Box) (cutoffSq : ℚ)
    (kLower kUpper : Fin m → ℚ) : RawContacts m :=
  { cutoffSq := cutoffSq
    d2Lower := fun i =>
      (squaredDistance (coords (edges i).source) (coords (edges i).target)).lower
    d2Upper := fun i =>
      (squaredDistance (coords (edges i).source) (coords (edges i).target)).upper
    kLower := kLower
    kUpper := kUpper }

/-- Squared distances computed from any coordinates in the boxes are
admissible for the contact gate, provided the stiffness interval holds. -/
theorem contactBoxes_admissible {n m : ℕ}
    (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (x : Fin n → Fin 3 → ℚ)
    (cutoffSq : ℚ) (kLower kUpper k : Fin m → ℚ)
    (hx : ∀ j a, contains (boxes j a) (x j a))
    (hk : ∀ i, kLower i ≤ k i ∧ k i ≤ kUpper i) :
    admissible (contactBoxes edges boxes cutoffSq kLower kUpper)
      (fun i => ∑ a : Fin 3, (x (edges i).source a -
        x (edges i).target a) ^ 2) k := by
  intro i
  have hs := squaredDistance_sound
    (boxes (edges i).source) (boxes (edges i).target)
    (x (edges i).source) (x (edges i).target)
    (fun a => hx (edges i).source a)
    (fun a => hx (edges i).target a)
  exact ⟨hs.1, hs.2, (hk i).1, (hk i).2⟩

/-- End-to-end exact-rational coordinate-box contact soundness: if Lean's
executable finite interval gate accepts, then EVERY rational coordinate
assignment in the source boxes produces spring energy between the two
certified endpoint models, even as contacts cross the cutoff. -/
theorem coordinate_contact_energy_sound {n m : ℕ}
    (edges : Fin m → Spring n)
    (boxes : Fin n → Fin 3 → Box)
    (cutoffSq : ℚ) (kLower kUpper : Fin m → ℚ)
    (u : Fin n → Fin 3 → ℝ)
    (hacc : accepts (contactBoxes edges boxes cutoffSq kLower kUpper) = true) :
    ∀ (x : Fin n → Fin 3 → ℚ), (∀ j a, contains (boxes j a) (x j a)) →
      ∀ (k : Fin m → ℚ), (∀ i, kLower i ≤ k i ∧ k i ≤ kUpper i) →
      springEnergy edges (fun i =>
        ((lowerWeight (contactBoxes edges boxes cutoffSq kLower kUpper) i) : ℝ)) u ≤
      springEnergy edges (fun i =>
        ((actualWeight (contactBoxes edges boxes cutoffSq kLower kUpper)
          (fun j => ∑ a : Fin 3,
            (x (edges j).source a - x (edges j).target a) ^ 2) k i) : ℝ)) u ∧
      springEnergy edges (fun i =>
        ((actualWeight (contactBoxes edges boxes cutoffSq kLower kUpper)
          (fun j => ∑ a : Fin 3,
            (x (edges j).source a - x (edges j).target a) ^ 2) k i) : ℝ)) u ≤
      springEnergy edges (fun i =>
        ((upperWeight (contactBoxes edges boxes cutoffSq kLower kUpper) i) : ℝ)) u := by
  intro x hx k hk
  exact all_admissible_energy_sound edges
    (contactBoxes edges boxes cutoffSq kLower kUpper) u hacc
    (fun i => ∑ a : Fin 3,
      (x (edges i).source a - x (edges i).target a) ^ 2) k
    (contactBoxes_admissible edges boxes x cutoffSq kLower kUpper k hx hk)

end FoldBounds.MotionCoordinates
