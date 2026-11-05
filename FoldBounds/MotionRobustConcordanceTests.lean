import FoldBounds.MotionRobustConcordance

open FoldBounds.MotionRobustConcordance

example : alignedSign (1 : ℚ) 2 = 1 := by
  norm_num [alignedSign]

example : alignedSign (2 : ℚ) 1 = -1 := by
  norm_num [alignedSign]

example : alignedSign (1 : ℚ) 1 = 0 := by
  norm_num [alignedSign]

#print axioms FoldBounds.MotionRobustConcordance.alignedSign_bounds
#print axioms FoldBounds.MotionRobustConcordance.alignedSign_interval
#print axioms FoldBounds.MotionRobustConcordance.sum_alignedSign_interval
