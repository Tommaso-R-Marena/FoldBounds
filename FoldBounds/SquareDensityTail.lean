import Mathlib

/-!
# A kernel-checked algebraic core for infinite-tail population certificates

This is the *generic consequence* of a nonnegative odds-majorant recurrence.
It does NOT assert that the Python SAW enumeration, geometric construction,
Boltzmann weight, or decoder satisfies that recurrence. Those remain outside
Lean's trust boundary until explicitly formalized and checked.

If T(0)=F0>=0 and T(k+1)=T(k)*r*d^k with r,d in [0,1], then every
T(k)<=F0. Whenever a model supplies event mass E>0 and complement mass
N>=0 with N<=E*T(k), event population is at least 1/(1+F0),
simultaneously for all k, without sampling finite ranges of k.
-/

namespace FoldBounds.SquareDensityTail

/-- A nonnegative geometric-ratio tail envelope. -/
def tailOdds (F0 r d : ℚ) : ℕ → ℚ
  | 0 => F0
  | k + 1 => tailOdds F0 r d k * r * d ^ k

/-- Powers of an exact rational in [0,1] remain at most one. -/
theorem power_le_one (d : ℚ) (hd : 0 ≤ d) (hu : d ≤ 1) (k : ℕ) :
    d ^ k ≤ 1 := by
  induction k with
  | zero => norm_num
  | succ k ih =>
    rw [pow_succ]
    have hp : 0 ≤ d ^ k := pow_nonneg hd k
    nlinarith [mul_nonneg (sub_nonneg.mpr hu) hp]

/-- Every geometric-ratio envelope term is nonnegative. -/
theorem tailOdds_nonneg (F0 r d : ℚ) (hF : 0 ≤ F0) (hr : 0 ≤ r)
    (hd : 0 ≤ d) (k : ℕ) : 0 ≤ tailOdds F0 r d k := by
  induction k with
  | zero => simpa [tailOdds] using hF
  | succ k ih =>
    change 0 ≤ tailOdds F0 r d k * r * d ^ k
    exact mul_nonneg (mul_nonneg ih hr) (pow_nonneg hd k)

/-- An exact rational recurrence is *uniformly* bounded by its initial term. -/
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

/-- If N/E <= F, the event probability is at least 1/(1+F). -/
theorem odds_probability_sound {E N F : ℚ} (hE : 0 < E)
    (hN : 0 ≤ N) (hF : 0 ≤ F) (h : N ≤ E * F) :
    (1 : ℚ) / (1 + F) ≤ E / (E + N) := by
  have hden1 : (0 : ℚ) < 1 + F := by linarith
  have hden2 : (0 : ℚ) < E + N := by linarith
  apply (div_le_div_iff₀ hden1 hden2).2
  nlinarith

/-- **Infinite-family soundness core:** one finite rational odds gate
implies a uniform population guarantee at every nonnegative index k,
provided the caller proves the model-specific N<=E*T(k) premise. -/
theorem uniform_population_bound (F0 r d : ℚ) (hF : 0 ≤ F0)
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
  exact odds_probability_sound hE hN hF hEF

-- The theorem has no admits, custom axioms, native_decide, or unsafe code.
#print axioms FoldBounds.SquareDensityTail.uniform_population_bound

end FoldBounds.SquareDensityTail
