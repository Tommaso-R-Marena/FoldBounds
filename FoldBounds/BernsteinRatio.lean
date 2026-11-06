import Mathlib

/-!
# A kernel-checked abstract convex-ratio bound

This is the finite weighted-average part of the Bernstein ensemble certificate.
No claim is made here that untrusted Python-generated coefficient grids
correspond to a protein Hamiltonian, or that Bernstein tensor coefficients
were computed from their polynomial inputs correctly. Those are separate
refinement obligations.
-/

namespace FoldBounds.BernsteinRatio

/-- If each coefficient term satisfies a local numerator/denominator ratio
    enclosure and each mixing weight is nonnegative, then the ratio of total
    weighted masses lies in that same enclosure whenever the total
    denominator is positive. Used by Bernstein convex-hull population bounds. -/
theorem weighted_ratio_enclosure {ι : Type*} (s : Finset ι)
    (event complement weight : ι → ℚ) (lo hi : ℚ)
    (hweight : ∀ i ∈ s, 0 ≤ weight i)
    (hden : 0 < ∑ i ∈ s, (event i + complement i) * weight i)
    (hlo : ∀ i ∈ s, lo * (event i + complement i) ≤ event i)
    (hhi : ∀ i ∈ s, event i ≤ hi * (event i + complement i)) :
    lo ≤ (∑ i ∈ s, event i * weight i) /
      (∑ i ∈ s, (event i + complement i) * weight i) ∧
    (∑ i ∈ s, event i * weight i) /
      (∑ i ∈ s, (event i + complement i) * weight i) ≤ hi := by
  have hLower :
      lo * (∑ i ∈ s, (event i + complement i) * weight i) ≤
      (∑ i ∈ s, event i * weight i) := by
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum
    intro i hi_mem
    have h := mul_le_mul_of_nonneg_right (hlo i hi_mem) (hweight i hi_mem)
    nlinarith [h]
  have hUpper :
      (∑ i ∈ s, event i * weight i) ≤
      hi * (∑ i ∈ s, (event i + complement i) * weight i) := by
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum
    intro i hi_mem
    have h := mul_le_mul_of_nonneg_right (hhi i hi_mem) (hweight i hi_mem)
    nlinarith [h]
  constructor
  · exact (le_div_iff₀ hden).2 hLower
  · exact (div_le_iff₀ hden).2 hUpper

end FoldBounds.BernsteinRatio