import Mathlib

/-!
# Order-statistic margin robustness for mutation-invariant distance codes

If at most B reporter channels have signed geometric margin <= delta,
then perturbing every pair-distance readout by strictly less than delta
can create no more than B wrong (or ambiguous) votes from geometry.
Allow a further T arbitrary wrong channels plus E erasures. Correct
majority is guaranteed if 2*(B+T)+E<K.

Python separately calculates the exact order statistic on ALL named
NMR conformers. These theorems are the sound *abstract inference*
from its bound, NOT end-to-end verified PDB decoding or the Lipschitz
Euclidean coordinate-distance calculations.
-/

namespace FoldBounds.MutationSafeMargins

/-- A strictly positive signed classification margin stays positive
under an additive perturbation of smaller absolute magnitude. -/
theorem signed_margin_stays_correct
    {margin perturbed delta : ℚ}
    (hPositive : delta < margin)
    (hError : |perturbed| ≤ delta) :
    0 < margin + perturbed := by
  have hLower : -delta ≤ perturbed := (abs_le.mp hError).1
  linarith

/-- A strict class threshold gap guarantees a distance decision under
symmetric signed measurement error, without probabilistic assumptions. -/
theorem stable_distance_threshold
    {trueDistance threshold observed margin noise : ℚ}
    (hTrue : trueDistance + margin ≤ threshold)
    (hNoise : |observed - trueDistance| ≤ noise)
    (hGap : noise < margin) :
    observed < threshold := by
  have hUp := (abs_le.mp hNoise).2
  linarith

/-- A robust code that has at most B geometrically ambiguous channels
and T further adversarially incorrect channels survives E missing
reporters whenever the **combined fault budget** is strict. -/
theorem quantile_error_erasure_majority
    {K E T B good wrong : ℕ}
    (hTotal : good + wrong + E = K)
    (hWrong : wrong ≤ B + T)
    (hBudget : 2*(B+T) + E < K) :
    wrong < good := by
  omega

/-- If the total wrong count includes existing, geometry-induced and
adversarial sources, they all consume the SAME error budget. -/
theorem three_source_error_budget
    {K E Bexisting Bgeometry T good wrong : ℕ}
    (hTotal : good + wrong + E = K)
    (hWrong : wrong ≤ Bexisting + Bgeometry + T)
    (hBudget : 2*(Bexisting+Bgeometry+T) + E < K) :
    wrong < good := by
  omega

/-- Strict geometric margin robustness is monotone in the declared
measurement error budget. -/
theorem smaller_perturbation_inherits_certificate
    {delta smaller bound : ℚ}
    (h : delta < bound)
    (hSmaller : smaller ≤ delta) :
    smaller < bound := by
  linarith

/-- The combined budget is tight: K=2*(B+T)+E admits an exact
tie with B+T wrong votes, B+T correct votes and E erasures. -/
theorem quantile_limit_is_sharp (B T E : ℕ) :
    ∃ good wrong : ℕ,
      good + wrong + E = 2*(B+T)+E ∧
      wrong = B+T ∧ good = wrong := by
  refine ⟨B+T,B+T,?_,rfl,rfl⟩
  omega

#print axioms FoldBounds.MutationSafeMargins.signed_margin_stays_correct
#print axioms FoldBounds.MutationSafeMargins.stable_distance_threshold
#print axioms FoldBounds.MutationSafeMargins.quantile_error_erasure_majority
#print axioms FoldBounds.MutationSafeMargins.three_source_error_budget
#print axioms FoldBounds.MutationSafeMargins.smaller_perturbation_inherits_certificate
#print axioms FoldBounds.MutationSafeMargins.quantile_limit_is_sharp

end FoldBounds.MutationSafeMargins
