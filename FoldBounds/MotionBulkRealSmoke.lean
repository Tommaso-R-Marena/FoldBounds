module

public import FoldBounds.MotionBulkChain

/-!
# Real frozen CspB / crambin interval-order smoke test

These are the five lowest-lower-endpoint frozen score intervals from each
protein, extracted from the SHA-256-verified 2026-11-04 audit archive.
They cover 4 adjacent comparisons per protein, under exact signed integer
arithmetic, without any numerical sorting oracle inside the Lean proof.

The archived endpoint-to-physical-score inclusion remains a separate,
unverified real-data obligation: these are ranking *arithmetic* tests.
-/

@[expose] public section

namespace FoldBounds.MotionBulkRealSmoke

open FoldBounds.MotionScaledCertificates
open FoldBounds.MotionBulkChain

def cspb : List ScaledInterval := [
  ⟨14712621271019661326076997341, 14712621271019677981865279856⟩,
  ⟨14807566128725553564387109138, 14807566128725571402186899464⟩,
  ⟨16260503532813406288034728204, 16260503532813418558151786880⟩,
  ⟨16508778583327083471082377410, 16508778583327097233462611725⟩,
  ⟨16703103069675000519680543264, 16703103069675014217150447418⟩
]

def crambin : List ScaledInterval := [
  ⟨15756618156690966767311523900, 15756618156690974349917324978⟩,
  ⟨16111180574249126797278368313, 16111180574249134258379367842⟩,
  ⟨16544332648389772331303999404, 16544332648389778334476827205⟩,
  ⟨17184025700892225792723436925, 17184025700892229266952058874⟩,
  ⟨17450141670767598317215934846, 17450141670767608600857274835⟩
]

theorem cspb_frozen_rank_arithmetic : checkAdjacent cspb = true := by native_decide
theorem crambin_frozen_rank_arithmetic :
    checkAdjacent crambin = true := by native_decide

theorem cspb_frozen_strict : StrictAdjacent cspb :=
  checkAdjacent_sound cspb cspb_frozen_rank_arithmetic

theorem crambin_frozen_strict : StrictAdjacent crambin :=
  checkAdjacent_sound crambin crambin_frozen_rank_arithmetic

end FoldBounds.MotionBulkRealSmoke
