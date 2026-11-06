import FoldBounds.FrontierNextFlat

/-!
# Kernel-checked tests of the version-2 frontier checker

Every equality is decided by the Lean kernel (`decide +kernel`; no
`native_decide`, no extra axioms).  The tests cover accepted certificates
with exact expected endpoints, the monotone effect of memory, the strict
improvement from disjoint family + witness counting, backward
compatibility with the version-1 leaves, and a catalogue of rejected
malformed or dishonest certificates.

`HHHHHH` (`steps = 4`) has 71 anchored conformations; at `q = 2` its exact values
are `Z = 117`, `Z_{k≥2} = 32`, `Z_{k<2} = 85`.
-/

open FoldBounds FoldBounds.FrontierNext FoldBounds.LatticeModel

namespace FoldBounds.FrontierNextTests

def hhhh : List Bool := [true, true, true, true]
def h6 : List Bool := List.replicate 6 true
def mix6 : List Bool := [true, false, true, true, false, true]

/-! ### Accepted certificates -/

/-- Exhaustive local evaluation at the root is exact (`HHHH`, `q = 2`). -/
theorem accept_enum_hhhh :
    check hhhh 2 1 .enum = some ⟨11, 11, 4, 4, 7, 7, 4/11, 4/11⟩ := by decide +kernel

/-- Exhaustive local evaluation at the root is exact (`HHHHHH`, `q = 2`). -/
theorem accept_enum_h6 :
    check h6 2 2 .enum = some ⟨117, 117, 32, 32, 85, 85, 32/117, 32/117⟩ := by decide +kernel

/-- Mixed sequence, `q = 3/2`. -/
theorem accept_enum_mix6 :
    check mix6 (3/2) 1 .enum = some ⟨175/2, 175/2, 81/2, 81/2, 47, 47, 81/175, 81/175⟩ := by
  decide +kernel

/-- `q = 1` counts conformations: the 6-residue chain has 71 anchored
conformations (no reflection quotient), 8 of them with at least two contacts. -/
theorem accept_enum_q1 :
    check h6 1 2 .enum = some ⟨71, 71, 8, 8, 63, 63, 8/71, 8/71⟩ := by decide +kernel

/-- Legacy leaves keep exactly their version-1 meaning. -/
theorem legacy_agrees_with_v1 :
    check h6 2 2 (.legacy true []) =
      FrontierCheck.check h6 2 2 (.bound true []) := by decide +kernel

theorem accept_root_legacy_transfer :
    check h6 2 2 (.legacy true []) = some ⟨16, 202, 0, 104, 0, 170, 0, 1⟩ := by decide +kernel

/-- One-node memory-0 certificate: geometry and degree caps already improve
the legacy transfer bound (`135 < 202`); the partially directed family
gives `Z ≥ 41` (versus `16` for the diagonal family). -/
theorem accept_root_mem0 :
    check h6 2 2 (.mem 0 []) = some ⟨41, 135, 0, 48, 0, 119, 0, 1⟩ := by decide +kernel

/-- Memory 2 is strictly tighter, and here exact for the total mass. -/
theorem accept_root_mem2 :
    check h6 2 2 (.mem 2 []) = some ⟨41, 117, 0, 32, 0, 85, 0, 1⟩ := by decide +kernel

/-- **Disjoint family + witness counting.**  The witness `2113` is not a
partially directed continuation, so its exact weight `4` is *added* to the
family mass `41` (kernel-checked disjointness): `Z ≥ 45 > max 41 4`. -/
theorem accept_root_mem3_witness :
    check h6 2 2 (.mem 3 [[2, 1, 1, 3]]) = some ⟨45, 117, 4, 32, 0, 85, 4/89, 1⟩ := by
  decide +kernel

/-- A mixed trie with every new leaf kind. -/
theorem accept_mixed_trie :
    check h6 2 2 (.node (.mem 2 []) .dead .enum (.legacy false [])) =
      some ⟨57, 509, 12, 452, 28, 489, 4/167, 113/120⟩ := by decide +kernel

/-- The same certificate in flat form. -/
theorem accept_flat :
    checkFlat2 h6 2 2 [([0], .mem 2 []), ([2], .enum), ([3], .legacy false [])] =
      some ⟨57, 509, 12, 452, 28, 489, 4/167, 113/120⟩ := by decide +kernel

/-! ### Rejected certificates -/

/-- `q < 1` is rejected. -/
theorem reject_small_q : check h6 (1/2) 2 (.mem 3 []) = none := by decide +kernel

/-- Too-short sequences are rejected. -/
theorem reject_short_seq : check [true] 2 0 .enum = none := by decide +kernel

/-- Duplicate witnesses are rejected (no double counting). -/
theorem reject_duplicate_witness :
    check h6 2 2 (.mem 3 [[2, 1, 1, 3], [2, 1, 1, 3]]) = none := by decide +kernel

/-- A colliding witness is rejected. -/
theorem reject_colliding_witness : check h6 2 2 (.mem 3 [[2, 1, 3, 3]]) = none := by
  decide +kernel

/-- A witness of the wrong length is rejected. -/
theorem reject_short_witness : check h6 2 2 (.mem 3 [[2, 1, 3]]) = none := by decide +kernel

/-- A nonterminal node marked exact is rejected. -/
theorem reject_false_exact : check h6 2 2 .exact = none := by decide +kernel

/-- A dead mark on a prefix with legal extensions is rejected. -/
theorem reject_false_dead : check h6 2 2 .dead = none := by decide +kernel

/-- A memory leaf on a colliding (reversal) prefix is rejected. -/
theorem reject_mem_on_collision :
    check h6 2 2 (.node (.mem 0 []) (.mem 0 []) .enum .enum) = none := by decide +kernel

/-- An exact leaf on a nonterminal child is rejected. -/
theorem reject_exact_child : check h6 2 2 (.node .enum .dead .enum .exact) = none := by
  decide +kernel

/-- A flat certificate omitting a legal branch is rejected. -/
theorem reject_flat_missing_branch :
    checkFlat2 h6 2 2 [([0], .mem 2 []), ([2], .enum)] = none := by decide +kernel

/-- A flat certificate with a duplicated entry is rejected. -/
theorem reject_flat_duplicate :
    checkFlat2 h6 2 2 [([0], .mem 2 []), ([0], .mem 2 []), ([2], .enum), ([3], .enum)] = none := by
  decide +kernel

/-! ### From acceptance to a verified exact probability -/

/-- The equilibrium probability of at least two H–H contacts for `HHHHHH`
at `q = 2` is exactly `32/117`, obtained from an accepted one-line
certificate and the flagship theorem `FrontierNext.check_sound`. -/
theorem h6_probability_exact :
    Executable.massQ (weightQ 4 2 (fun _ => true)) (contactEvent 4 (fun _ => true) 2) /
      Executable.massQ (weightQ 4 2 (fun _ => true)) Finset.univ = 32 / 117 := by
  have hseq : List.ofFn (fun _ : Fin (4 + 2) => true) = h6 := by decide
  have h := check_sound 4 (fun _ => true) 2 2 .enum ⟨117, 117, 32, 32, 85, 85, 32/117, 32/117⟩
    (by rw [hseq]; exact accept_enum_h6)
  exact le_antisymm h.2.2.2.2 h.2.2.2.1

end FoldBounds.FrontierNextTests
