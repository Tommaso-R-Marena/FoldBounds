module

public import FoldBounds.MotionProteinResponse

/-!
# Frozen real-protein coordinate samples (exact IEEE-754 rationals)

These four C-alpha coordinates were exported from the frozen reference NPZs in
archive SHA-256
666a5d734fbd61be0f182c954474841bb552b71b3a1e9156d6c21b8c2f095f95.

Each sample is the FIRST eligible reference pair, zero-based array
positions (1,5), corresponding to PDB residue IDs (2,6) for these archives.
The exact rational input and squared-distance arithmetic are checked by Lean.

IMPORTANT: PDB/NPZ-to-Lean source provenance is NOT kernel-checked here.
The exporter/manifest remain separate trust boundary. This does not certify
full-protein motion scores for the original datasets.
-/

@[expose] public section

namespace FoldBounds.MotionRealPairSamples

open FoldBounds.MotionVariableDirections
open FoldBounds.MotionUncertainty

def pair : Spring 2 :=
  ⟨0, 1, ![1, 0, 0]⟩

/-- CspB crystal 1CSP, reference array positions 1 and 5. -/
def cspb : Fin 2 → Fin 3 → ℚ :=
  ![![1108607 / 131072, 370475 / 16384, 10617881 / 1048576],
    ![9162457 / 524288, 8750891 / 524288, 2173567 / 131072]]

/-- Crambin crystal 1CRN, reference array positions 1 and 5. -/
def crambin : Fin 2 → Fin 3 → ℚ :=
  ![![14529069 / 1048576, 6013059 / 524288, 3180331 / 524288],
    ![2273575 / 262144, 5572133 / 1048576, 16021193 / 1048576]]

theorem cspb_frozen_sq_distance :
    distSq pair cspb = (173805739498313 : ℚ) / 1099511627776 := by
  norm_num [distSq, Fin.sum_univ_three, pair, cspb, crambin]

theorem crambin_frozen_sq_distance :
    distSq pair crambin = (164516495665547 : ℚ) / 1099511627776 := by
  norm_num [distSq, Fin.sum_univ_three, pair, cspb, crambin]

theorem cspb_reference_pair_eligible :
    (25 : ℚ) ≤ distSq pair cspb ∧ distSq pair cspb ≤ 256 := by
  rw [cspb_frozen_sq_distance]
  constructor <;> norm_num

theorem crambin_reference_pair_eligible :
    (25 : ℚ) ≤ distSq pair crambin ∧ distSq pair crambin ≤ 256 := by
  rw [crambin_frozen_sq_distance]
  constructor <;> norm_num

end FoldBounds.MotionRealPairSamples
