module

public import FoldBounds.MotionScaledCertificates

/-!
# FoldBounds.MotionBulkChain

Kernel-checked ranking arithmetic for arbitrary lists of 2^96-scaled
score enclosures, with no external solver or Python arithmetic in the
acceptance logic.

A chain accepts only if EVERY neighboring interval is strictly
separated. Passing a chain certificate does NOT by itself establish
that the proposed intervals enclose the actual physical-model scores:
the per-pair MotionScaledCertificates checker must separately accept
each interval. This deliberately separates ranking arithmetic from
the more expensive 3N x 3N score verification.
-/

@[expose] public section

namespace FoldBounds.MotionBulkChain

open FoldBounds.MotionScaledCertificates

/-- Exact-integer comparison of two adjacent proposed score enclosures. -/
def separated (a b : ScaledInterval) : Bool :=
  decide (a.upper < b.lower)

/-- All adjacent comparisons, tested using signed integer arithmetic.
The trailing interval overlaps the first interval of the next chunk
when exported in the multi-chunk replay format. -/
def checkAdjacent : List ScaledInterval → Bool
  | [] => true
  | [_] => true
  | a :: b :: xs =>
      separated a b && checkAdjacent (b :: xs)
termination_by
  xs => xs.length + 1

/-- Logical meaning of an all-adjacent integer separation certificate. -/
inductive StrictAdjacent : List ScaledInterval → Prop
  | nil : StrictAdjacent []
  | singleton (a : ScaledInterval) : StrictAdjacent [a]
  | cons {a b : ScaledInterval} {xs : List ScaledInterval}
      (h : a.upper < b.lower) (htail : StrictAdjacent (b :: xs)) :
      StrictAdjacent (a :: b :: xs)

/-- Every accepted adjacent comparison really is strict. -/
theorem separated_sound (a b : ScaledInterval)
    (h : separated a b = true) : a.upper < b.lower :=
  of_decide_eq_true h

/-- The recursive executable checker implies the full chain of strict
adjacent integer inequalities; it does not trust the proposed endpoint
ordering or any precomputed sorting witness. -/
theorem checkAdjacent_sound :
    ∀ xs : List ScaledInterval,
      checkAdjacent xs = true → StrictAdjacent xs
  | [], _ => StrictAdjacent.nil
  | [a], _ => StrictAdjacent.singleton a
  | a :: b :: xs, h => by
      have hfirst : separated a b = true := by
        cases hd : separated a b with
        | false => simp [checkAdjacent, hd] at h
        | true => rfl
      have hrest : checkAdjacent (b :: xs) = true := by
        simpa [checkAdjacent, hfirst] using h
      exact StrictAdjacent.cons (separated_sound a b hfirst)
        (checkAdjacent_sound (b :: xs) hrest)
termination_by
  xs => xs.length

/-- Strict separation in an integer unit implies strict separation in
the exact rational 2^-96 units used by the protein-response theorems. -/
theorem separated_scaled_sound (a b : ScaledInterval)
    (h : separated a b = true) :
    upperQ a < lowerQ b := by
  have hi := separated_sound a b h
  have hq : (a.upper : ℚ) < (b.lower : ℚ) := by exact_mod_cast hi
  unfold upperQ lowerQ
  apply (div_lt_div_iff₀ binaryScale_pos binaryScale_pos).2
  exact mul_lt_mul_of_pos_right hq binaryScale_pos

/-- A successful pair comparison plus the two per-pair accepted score
bounds entails an actual ordering of the true rational model scores. -/
theorem scores_ordered_of_certificates
    (a b : ScaledInterval)
    (qa qb : ℚ)
    (ha : lowerQ a ≤ qa ∧ qa ≤ upperQ a)
    (hb : lowerQ b ≤ qb ∧ qb ≤ upperQ b)
    (h : separated a b = true) : qa < qb :=
  lt_of_lt_of_le
    (lt_of_le_of_lt ha.2 (separated_scaled_sound a b h))
    hb.1

end FoldBounds.MotionBulkChain
