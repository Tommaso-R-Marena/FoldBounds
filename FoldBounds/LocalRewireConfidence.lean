import FoldBounds

/-!
# Conditional exact-rational structural-confidence ceiling

This proves the *generic mass-to-probability implication* used by the
Python 3D local-rewiring confidence checker:

If the proposed target's true mass T is <=hi*W and the mass O of
all other conformations is >=lo*W*(S-1), where W,lo,hi>0 and S>=1,
then the Gibbs probability of that exact target is <=

    hi / (hi + lo*(S-1)).

The checker separately asserts that S was constructed using mutually
nonadjacent cube rewiring patches. **Lean does not yet prove the
3D geometry, exact local contact-weight calculations, or that S
matches the Python-reconstructed family.** Those remain separate
assurance obligations.
-/

namespace FoldBounds.LocalRewireConfidence

/-- Nonnegative constructive alternative mass limits the target's
equilibrium probability, even with arbitrary microscopic corrections. -/
theorem confidence_ceiling
    {T O W S lo hi : ℚ}
    (hT : 0 < T) (hW : 0 < W)
    (hlo : 0 < lo) (hhi : 0 < hi)
    (hS : 1 ≤ S)
    (hTarget : T ≤ hi * W)
    (hOther : lo * W * (S - 1) ≤ O) :
    T / (T + O) ≤ hi / (hi + lo * (S - 1)) := by
  have hSm : (0 : ℚ) ≤ S - 1 := by linarith
  have hB : (0 : ℚ) ≤ lo * W * (S - 1) :=
    mul_nonneg (mul_nonneg hlo.le hW.le) hSm
  have hOp : 0 ≤ O := le_trans hB hOther
  have hA : 0 < hi * W := mul_pos hhi hW
  have hDen1 : 0 < T + O := by linarith
  have hDen2 : 0 < hi * W + lo * W * (S - 1) := by linarith
  have hScaled :
      T / (T + O) ≤
        (hi * W) / (hi * W + lo * W * (S - 1)) :=
    FoldBounds.ratio_le_ratio hT.le hB hTarget hOther hDen1 hDen2
  have hD : 0 < hi + lo * (S - 1) := by
    have : 0 ≤ lo * (S - 1) := mul_nonneg hlo.le hSm
    linarith
  have hcancel :
      (hi * W) / (hi * W + lo * W * (S - 1)) =
        hi / (hi + lo * (S - 1)) := by
    have hden : hi * W + lo * W * (S - 1) =
        W * (hi + lo * (S - 1)) := by ring
    rw [hden]
    field_simp [ne_of_gt hW, ne_of_gt hD]
  calc
    T / (T + O) ≤ (hi * W) / (hi * W + lo * W * (S - 1)) := hScaled
    _ = hi / (hi + lo * (S - 1)) := hcancel

/-- Equal-energy alternatives are the special case S=D with no
microscopic correction uncertainty. The cap 1/D is independent of
temperature and the baseline target weight. -/
theorem isoenergetic_ceiling
    {T O W D : ℚ}
    (hT : 0 < T) (hW : 0 < W)
    (hD : 1 ≤ D) (hTarget : T ≤ W)
    (hOther : W * (D - 1) ≤ O) :
    T / (T + O) ≤ 1 / D := by
  have h := confidence_ceiling hT hW (show (0:ℚ)<1 by norm_num)
      (show (0:ℚ)<1 by norm_num) hD (by simpa using hTarget)
      (by simpa using hOther)
  simpa [show (1 : ℚ) + (D - 1) = D by ring] using h

/-- A model-independent strict no-majority test. -/
theorem excludes_majority
    {T O W S lo hi : ℚ}
    (hT : 0 < T) (hW : 0 < W)
    (hlo : 0 < lo) (hhi : 0 < hi)
    (hS : 1 ≤ S) (hTarget : T ≤ hi * W)
    (hOther : lo * W * (S - 1) ≤ O)
    (hgap : hi < lo * (S - 1)) :
    T / (T + O) < 1 / 2 := by
  have hcap := confidence_ceiling hT hW hlo hhi hS hTarget hOther
  have hSm : (0 : ℚ) ≤ S - 1 := by linarith
  have hden : 0 < hi + lo * (S - 1) := by
    have : 0 ≤ lo * (S - 1) := mul_nonneg hlo.le hSm
    linarith
  have hhalf : hi / (hi + lo * (S - 1)) < (1 / 2 : ℚ) := by
    apply (div_lt_div_iff₀ hden (by norm_num : (0 : ℚ) < 2)).2
    nlinarith
  exact lt_of_le_of_lt hcap hhalf

#print axioms FoldBounds.LocalRewireConfidence.confidence_ceiling
#print axioms FoldBounds.LocalRewireConfidence.isoenergetic_ceiling
#print axioms FoldBounds.LocalRewireConfidence.excludes_majority

end FoldBounds.LocalRewireConfidence
