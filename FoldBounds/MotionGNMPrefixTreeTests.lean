import FoldBounds.MotionGNMPrefixTree

/-! Reported rational witness data; this tests tree coverage and positivity,
not the separate 124-residue numeric producer or the NMR interpretation. -/

namespace FoldBounds.MotionGNMPrefixTreeTests
open FoldBounds.MotionGNMPrefixTree

private def reportedRationalTree : SplitTree ℚ :=
(SplitTree.branch
 (SplitTree.branch
 (SplitTree.branch
 (SplitTree.leaf (497/3415191 : ℚ))
 (SplitTree.branch
 (SplitTree.leaf (12329/3415191 : ℚ))
 (SplitTree.leaf (15779/3415191 : ℚ))))
 (SplitTree.branch
 (SplitTree.branch
 (SplitTree.leaf (23573/3415191 : ℚ))
 (SplitTree.leaf (27499/3415191 : ℚ)))
 (SplitTree.branch
 (SplitTree.leaf (11495/3415191 : ℚ))
 (SplitTree.leaf (14959/3415191 : ℚ)))))
 (SplitTree.branch
 (SplitTree.branch
 (SplitTree.leaf (1043/1138397 : ℚ))
 (SplitTree.branch
 (SplitTree.leaf (14801/3415191 : ℚ))
 (SplitTree.leaf (18359/3415191 : ℚ))))
 (SplitTree.branch
 (SplitTree.leaf (2369/3415191 : ℚ))
 (SplitTree.branch
 (SplitTree.leaf (4689/1138397 : ℚ))
 (SplitTree.leaf (263/50973 : ℚ))))))

private theorem all_reported_leaf_counts_positive :
    AllLeaves (fun q : ℚ => 0 < q) reportedRationalTree := by
  norm_num [reportedRationalTree, AllLeaves]

example (world : ℕ → Bool) : 0 < pickLeaf reportedRationalTree world 0 := by
  exact pickLeaf_sound (fun q : ℚ => 0 < q) reportedRationalTree world 0
    all_reported_leaf_counts_positive

#print axioms FoldBounds.MotionGNMPrefixTree.pickLeaf_sound
#print axioms FoldBounds.MotionGNMPrefixTree.leafCount_eq_splitCount_add_one
#print axioms FoldBounds.MotionGNMPrefixTree.lowerVote_mono
#print axioms FoldBounds.MotionGNMPrefixTree.leafLower_mono
#print axioms FoldBounds.MotionGNMPrefixTree.intersect_valid
#print axioms FoldBounds.MotionGNMPrefixTree.prefixTree_concordance_sound

end FoldBounds.MotionGNMPrefixTreeTests