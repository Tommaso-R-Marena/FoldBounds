import Mathlib

/-!
# Sound population reweighting under adversarial energy errors

For any finite physical model whose event and complement masses are independently
certified to lie in positive rational intervals, this algebraic lemma bounds
the reweighted population after arbitrary per-conformation multiplicative
Boltzmann-weight errors. It **does not** certify the input mass enclosures,
calibrate a physical force field, or prove that Python generated values match
the formal mathematical inputs. Its strict positivity hypotheses exclude the
zero-lower-mass cases, handled conservatively in Python.
-/

namespace FoldBounds.RobustReweighting

/-- If event and complement masses have positive certified lower bounds and
    their multiplicative correction factors lie in a common positive interval,
    then the reweighted event population lies in the advertised enclosure. -/
theorem robust_population_bounds
    (eventLo eventHi compLo compHi event comp facLo facHi facE facC : ℚ)
    (hEventLo : 0 < eventLo) (hCompLo : 0 < compLo)
    (heLo : eventLo ≤ event) (heHi : event ≤ eventHi)
    (hcLo : compLo ≤ comp) (hcHi : comp ≤ compHi)
    (hfLo : 0 < facLo)
    (hfELo : facLo ≤ facE) (hfEHi : facE ≤ facHi)
    (hfCLo : facLo ≤ facC) (hfCHi : facC ≤ facHi) :
    facLo * eventLo / (facLo * eventLo + facHi * compHi) ≤
      facE * event / (facE * event + facC * comp) ∧
    facE * event / (facE * event + facC * comp) ≤
      facHi * eventHi / (facHi * eventHi + facLo * compLo) := by
  have hEvent : 0 < event := lt_of_lt_of_le hEventLo heLo
  have hComp : 0 < comp := lt_of_lt_of_le hCompLo hcLo
  have hEventHi : 0 < eventHi := lt_of_lt_of_le hEvent heHi
  have hCompHi : 0 < compHi := lt_of_lt_of_le hComp hcHi
  have hfE : 0 < facE := lt_of_lt_of_le hfLo hfELo
  have hfC : 0 < facC := lt_of_lt_of_le hfLo hfCLo
  have hfHi : 0 < facHi := lt_of_lt_of_le hfE hfEHi
  have hdenLo : 0 < facLo * eventLo + facHi * compHi := by positivity
  have hden : 0 < facE * event + facC * comp := by positivity
  have hdenHi : 0 < facHi * eventHi + facLo * compLo := by positivity
  have hEventDown : facLo * eventLo ≤ facE * event := by
    calc
      facLo * eventLo ≤ facE * eventLo :=
        mul_le_mul_of_nonneg_right hfELo (le_of_lt hEventLo)
      _ ≤ facE * event :=
        mul_le_mul_of_nonneg_left heLo (le_of_lt hfE)
  have hCompUp : facC * comp ≤ facHi * compHi := by
    calc
      facC * comp ≤ facHi * comp :=
        mul_le_mul_of_nonneg_right hfCHi (le_of_lt hComp)
      _ ≤ facHi * compHi :=
        mul_le_mul_of_nonneg_left hcHi (le_of_lt hfHi)
  have hProdLo : (facLo * eventLo) * (facC * comp) ≤
        (facE * event) * (facHi * compHi) := by
    calc
      (facLo * eventLo) * (facC * comp) ≤
          (facE * event) * (facC * comp) :=
        mul_le_mul_of_nonneg_right hEventDown (le_of_lt (mul_pos hfC hComp))
      _ ≤ (facE * event) * (facHi * compHi) :=
        mul_le_mul_of_nonneg_left hCompUp (le_of_lt (mul_pos hfE hEvent))
  have hEventUp : facE * event ≤ facHi * eventHi := by
    calc
      facE * event ≤ facHi * event :=
        mul_le_mul_of_nonneg_right hfEHi (le_of_lt hEvent)
      _ ≤ facHi * eventHi :=
        mul_le_mul_of_nonneg_left heHi (le_of_lt hfHi)
  have hCompDown : facLo * compLo ≤ facC * comp := by
    calc
      facLo * compLo ≤ facC * compLo :=
        mul_le_mul_of_nonneg_right hfCLo (le_of_lt hCompLo)
      _ ≤ facC * comp :=
        mul_le_mul_of_nonneg_left hcLo (le_of_lt hfC)
  have hProdHi : (facE * event) * (facLo * compLo) ≤
        (facHi * eventHi) * (facC * comp) := by
    calc
      (facE * event) * (facLo * compLo) ≤
          (facHi * eventHi) * (facLo * compLo) :=
        mul_le_mul_of_nonneg_right hEventUp (le_of_lt (mul_pos hfLo hCompLo))
      _ ≤ (facHi * eventHi) * (facC * comp) :=
        mul_le_mul_of_nonneg_left hCompDown (le_of_lt (mul_pos hfHi hEventHi))
  constructor
  · apply (div_le_div_iff₀ hdenLo hden).2
    nlinarith [hProdLo]
  · apply (div_le_div_iff₀ hden hdenHi).2
    nlinarith [hProdHi]

end FoldBounds.RobustReweighting