import FoldBounds.FrontierFlat

/-!
# Concrete kernel-checked tests of the frontier checker

All equalities below are decided by the Lean *kernel* (`decide +kernel`,
no `native_decide`, no extra axioms).  They cover accepted certificates
with exact expected endpoints and a catalogue of rejected malformed or
dishonest certificates.  The last theorem turns an accepted certificate
into a kernel-verified *exact* equilibrium probability via `check_sound`.

Sequence `HHHH` (`steps = 2`) has 9 conformations; exactly two (the two
U-shapes) have one H–H contact, so for `q = 2`: `Z = 7 + 2·2 = 11`,
`Z_{k≥1} = 4`, `Z_{k<1} = 7`, `P(k ≥ 1) = 4/11`.
-/

open FoldBounds FoldBounds.FrontierCheck FoldBounds.LatticeModel

namespace FoldBounds.FrontierTests

def hhhh : List Bool := [true, true, true, true]

/-- Complete trie for `HHHH`: directions `0:+x 1:-x 2:+y 3:-y`; reversal
children (`1` at the root, the opposite of the previous step below) are
`dead` because they collide. -/
def fullTree : FTree :=
  .node (.node .exact .dead .exact .exact) .dead
    (.node .exact .exact .exact .dead) (.node .exact .exact .dead .exact)

/-! ### Accepted certificates -/

theorem accept_full_exact :
    check hhhh 2 1 fullTree = some ⟨11, 11, 4, 4, 7, 7, 4/11, 4/11⟩ := by decide +kernel

/-- One-node certificate: basic bound `3^2 · 2^1 = 18`. -/
theorem accept_root_basic :
    check hhhh 2 1 (.bound false []) = some ⟨0, 18, 0, 18, 0, 18, 0, 1⟩ := by decide +kernel

/-- One-node certificate with the verified transfer relaxation: here the
relaxation is already exact for the total mass, and the certified diagonal
family (`2^2` directed completions, weight `≥ q^0`) gives `Z ≥ 4`. -/
theorem accept_root_transfer :
    check hhhh 2 1 (.bound true []) = some ⟨4, 11, 0, 4, 0, 7, 0, 1⟩ := by decide +kernel

/-- Threshold `0`: every completion is in the event, so the family lower
bound also lower-bounds the event mass. -/
theorem accept_root_transfer_thr0 :
    check hhhh 2 0 (.bound true []) = some ⟨4, 11, 4, 11, 0, 0, 1, 1⟩ := by decide +kernel

/-- Threshold `3` exceeds the contact cap: the event is certified empty and the
family lower bound goes to the complement; `P = 0` exactly. -/
theorem accept_root_transfer_thr3 :
    check hhhh 2 3 (.bound true []) = some ⟨4, 11, 0, 0, 4, 11, 0, 0⟩ := by decide +kernel

/-- Witnesses raise lower bounds and tighten upper bounds:
`18 - (2 - 1) - (2 - 2) = 17`. -/
theorem accept_root_basic_witnesses :
    check hhhh 2 1 (.bound false [[0, 0], [2, 1]]) =
      some ⟨3, 17, 2, 16, 1, 15, 2/17, 16/17⟩ := by decide +kernel

theorem accept_root_transfer_witnesses :
    check hhhh 2 1 (.bound true [[0, 0], [2, 1]]) =
      some ⟨4, 11, 2, 4, 1, 7, 2/9, 4/5⟩ := by decide +kernel

theorem accept_rational_q :
    check hhhh (3/2) 1 (.bound true [[0, 0], [2, 1]]) =
      some ⟨4, 10, 3/2, 3, 1, 7, 3/17, 3/4⟩ := by decide +kernel

/-- Mixed frontier given in the flat (Python-style) format. -/
theorem accept_flat_mixed :
    checkFlat hhhh 2 1 [([0, 0], .exact), ([0, 2], .exact), ([0, 3], .exact),
      ([2], .bound false []), ([3], .bound true [])] =
      some ⟨5, 13, 0, 8, 3, 11, 0, 8/11⟩ := by decide +kernel

/-! ### Rejected certificates -/

/-- `exact` at a nonterminal prefix. -/
theorem reject_exact_nonterminal : check hhhh 2 1 .exact = none := by decide +kernel

/-- `dead` on a prefix that has valid extensions. -/
theorem reject_dead_with_extensions : check hhhh 2 1 .dead = none := by decide +kernel

/-- An omitted nonempty branch (a legal child marked `dead`). -/
theorem reject_omitted_branch :
    check hhhh 2 1 (.node (.node .exact .dead .exact .dead) .dead
      (.node .exact .exact .exact .dead) (.node .exact .exact .dead .exact)) = none := by
  decide +kernel

/-- A colliding (immediate reversal) child claimed as a real conformation. -/
theorem reject_invalid_geometry :
    check hhhh 2 1 (.node (.node .exact .exact .exact .exact) .dead
      (.node .exact .exact .exact .dead) (.node .exact .exact .dead .exact)) = none := by
  decide +kernel

/-- Internal node at a terminal prefix. -/
theorem reject_node_too_deep :
    check hhhh 2 1 (.node (.node (.node .dead .dead .dead .dead) .dead .exact .exact) .dead
      (.node .exact .exact .exact .dead) (.node .exact .exact .dead .exact)) = none := by
  decide +kernel

theorem reject_duplicate_witness :
    check hhhh 2 1 (.bound false [[0, 0], [0, 0]]) = none := by decide +kernel

theorem reject_self_intersecting_witness :
    check hhhh 2 1 (.bound false [[1, 0]]) = none := by decide +kernel

theorem reject_incomplete_witness :
    check hhhh 2 1 (.bound false [[0]]) = none := by decide +kernel

theorem reject_q_below_one : check hhhh (1/2) 1 (.bound false []) = none := by decide +kernel

theorem reject_short_sequence : check [true] 2 1 (.bound false []) = none := by decide +kernel

/-- Flat format: duplicate frontier node. -/
theorem reject_flat_duplicate :
    checkFlat hhhh 2 1 [([0, 0], .exact), ([0, 2], .exact), ([0, 3], .exact),
      ([2], .bound false []), ([3], .bound true []), ([3], .dead)] = none := by decide +kernel

/-- Flat format: hidden descendant of an already covered node. -/
theorem reject_flat_overlap :
    checkFlat hhhh 2 1 [([0, 0], .exact), ([0, 2], .exact), ([0, 3], .exact),
      ([2], .bound false []), ([3], .bound true []), ([3, 0], .exact)] = none := by
  decide +kernel

/-- Flat format: omitted nonempty branch. -/
theorem reject_flat_missing :
    checkFlat hhhh 2 1 [([0, 0], .exact), ([0, 2], .exact),
      ([2], .bound false []), ([3], .bound true [])] = none := by decide +kernel

/-! ### From an accepted certificate to a verified exact probability -/

def allH4 : Fin 4 → Bool := fun _ => true

theorem ofFn_allH4 : List.ofFn allH4 = hhhh := by decide

/-- Kernel-verified exact equilibrium probability `P(k ≥ 1) = 4/11` for
`HHHH` at `q = 2`, obtained from the accepted certificate via the
flagship soundness theorem. -/
theorem hhhh_probability_exact :
    Executable.massQ (weightQ 2 2 allH4) (contactEvent 2 allH4 1) /
      Executable.massQ (weightQ 2 2 allH4) Finset.univ = 4 / 11 := by
  have h : check (List.ofFn allH4) 2 1 fullTree = some ⟨11, 11, 4, 4, 7, 7, 4/11, 4/11⟩ := by
    rw [ofFn_allH4]; exact accept_full_exact
  have hs := (check_sound 2 allH4 2 1 fullTree _ h).2.2.2
  exact le_antisymm hs.2 hs.1

end FoldBounds.FrontierTests
