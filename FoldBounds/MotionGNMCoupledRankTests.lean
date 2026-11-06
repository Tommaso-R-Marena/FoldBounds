module
public import FoldBounds.MotionGNMCoupledRank

@[expose] public section
namespace FoldBounds.MotionGNMCoupledRankTests

open FoldBounds.MotionGNMCoupledRank

-- Each of K and M-K can be the identity (all cross terms zero).
-- The coupled sign test certifies any endpoint-midpoint sum s>2.
example : (0:ℚ) < (3/2:ℚ) := by
  apply coupled_strict_rank_sound
    (1:ℚ) 0 1 (1:ℚ) 0 1
    (1:ℚ) 1 0 3 (3/2:ℚ)
  · norm_num
  · intro t
    nlinarith [sq_nonneg t]
  · intro t
    nlinarith [sq_nonneg t]
  · norm_num
  · norm_num
  · norm_num
  · norm_num
  · norm_num

#print axioms FoldBounds.MotionGNMCoupledRank.quadratic_discriminant_bound
#print axioms FoldBounds.MotionGNMCoupledRank.psd_pair_contraction
#print axioms FoldBounds.MotionGNMCoupledRank.positive_of_coupled_square
#print axioms FoldBounds.MotionGNMCoupledRank.coupled_strict_rank_sound

end FoldBounds.MotionGNMCoupledRankTests