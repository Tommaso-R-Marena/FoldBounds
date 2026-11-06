import Mathlib

/-!
# Generic Lean kernel consequence of certified finite-family mass enclosures

This module does not formalize the cubic lattice, local tile enumeration,
or Python code. It proves the abstract positive-rational inference from
event/complement masses to population guarantees, including a uniform
infinite-tail consequence of a checked shrinking-odds recurrence.
-/

namespace FoldBounds.CubicFamilyMass

/-- If event and complement masses have separate positive bounds,
the coupled lower probability endpoint is sound. -/
theorem probability_lower_of_masses {E N L U : ℚ}
    (hL : 0 < L) (hLE : L ≤ E) (hN : 0 ≤ N) (hNU : N ≤ U) :
    L / (L + U) ≤ E / (E + N) := by
  have hU : 0 ≤ U := le_trans hN hNU
  have hd1 : 0 < L + U := by linarith
  have hd2 : 0 < E + N := by linarith
  apply (div_le_div_iff₀ hd1 hd2).2
  have h1 : L * N ≤ L * U := mul_le_mul_of_nonneg_left hNU hL.le
  have h2 : L * U ≤ E * U := mul_le_mul_of_nonneg_right hLE hU
  nlinarith

/-- A strict lower-event versus upper-complement gap certifies a majority. -/
theorem majority_of_mass_gap {E N L U : ℚ}
    (hN : 0 ≤ N) (hLE : L ≤ E) (hNU : N ≤ U) (hgap : U < L) :
    (1 / 2 : ℚ) < E / (E + N) := by
  have hE : 0 < E := by linarith
  have hden : 0 < E + N := by linarith
  apply (lt_div_iff₀ hden).2
  nlinarith

/-- An exact-rational odds majorant with geometric decay of successive ratios. -/
def tailOdds (F0 r d : ℚ) : ℕ → ℚ
  | 0 => F0
  | k + 1 => tailOdds F0 r d k * r * d ^ k

theorem power_le_one (d : ℚ) (hd : 0 ≤ d) (hu : d ≤ 1) (k : ℕ) :
    d ^ k ≤ 1 := by
  induction k with
  | zero => norm_num
  | succ k ih =>
    rw [pow_succ]
    have hp : 0 ≤ d ^ k := pow_nonneg hd k
    nlinarith [mul_nonneg (sub_nonneg.mpr hu) hp]

theorem tailOdds_nonneg (F0 r d : ℚ) (hF : 0 ≤ F0) (hr : 0 ≤ r)
    (hd : 0 ≤ d) (k : ℕ) : 0 ≤ tailOdds F0 r d k := by
  induction k with
  | zero => simpa [tailOdds] using hF
  | succ k ih =>
    change 0 ≤ tailOdds F0 r d k * r * d ^ k
    exact mul_nonneg (mul_nonneg ih hr) (pow_nonneg hd k)

theorem tailOdds_le_initial (F0 r d : ℚ) (hF : 0 ≤ F0)
    (hr : 0 ≤ r) (hr1 : r ≤ 1)
    (hd : 0 ≤ d) (hd1 : d ≤ 1) (k : ℕ) :
    tailOdds F0 r d k ≤ F0 := by
  induction k with
  | zero => simp [tailOdds]
  | succ k ih =>
    have ht : 0 ≤ tailOdds F0 r d k :=
      tailOdds_nonneg F0 r d hF hr hd k
    have hpr : d ^ k ≤ 1 := power_le_one d hd hd1 k
    have htR : 0 ≤ tailOdds F0 r d k * r := mul_nonneg ht hr
    have h1 : tailOdds F0 r d k * r ≤ tailOdds F0 r d k := by
      nlinarith [mul_nonneg ht (sub_nonneg.mpr hr1)]
    have h2 : tailOdds F0 r d k * r * d ^ k ≤ tailOdds F0 r d k * r := by
      nlinarith [mul_nonneg htR (sub_nonneg.mpr hpr)]
    calc
      tailOdds F0 r d (k + 1) =
          tailOdds F0 r d k * r * d ^ k := rfl
      _ ≤ tailOdds F0 r d k * r := h2
      _ ≤ tailOdds F0 r d k := h1
      _ ≤ F0 := ih

/-- The event's exact population is bounded by the complement/event odds cap. -/
theorem population_from_odds {E N F : ℚ}
    (hE : 0 < E) (hN : 0 ≤ N) (hF : 0 ≤ F)
    (hNE : N ≤ E * F) :
    (1 : ℚ) / (1 + F) ≤ E / (E + N) := by
  have hden1 : (0 : ℚ) < 1 + F := by linarith
  have hden2 : (0 : ℚ) < E + N := by linarith
  apply (div_le_div_iff₀ hden1 hden2).2
  nlinarith

/-- Uniform infinite-index consequence, conditional on a valid odds premise.

The concrete 3D SAW counting and H/P contact inequalities are still
verified by independent Python calculations, not by this Lean theorem.
-/
theorem uniform_population_from_tail (F0 r d : ℚ) (hF : 0 ≤ F0)
    (hr : 0 ≤ r) (hr1 : r ≤ 1)
    (hd : 0 ≤ d) (hd1 : d ≤ 1)
    (k : ℕ) {E N : ℚ} (hE : 0 < E) (hN : 0 ≤ N)
    (hOdds : N ≤ E * tailOdds F0 r d k) :
    (1 : ℚ) / (1 + F0) ≤ E / (E + N) := by
  have htail : tailOdds F0 r d k ≤ F0 :=
    tailOdds_le_initial F0 r d hF hr hr1 hd hd1 k
  have hEF : N ≤ E * F0 := by
    calc
      N ≤ E * tailOdds F0 r d k := hOdds
      _ ≤ E * F0 := mul_le_mul_of_nonneg_left htail hE.le
  exact population_from_odds hE hN hF hEF

#print axioms FoldBounds.CubicFamilyMass.probability_lower_of_masses
#print axioms FoldBounds.CubicFamilyMass.majority_of_mass_gap
#print axioms FoldBounds.CubicFamilyMass.uniform_population_from_tail

end FoldBounds.CubicFamilyMass
