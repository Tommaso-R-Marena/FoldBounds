module

public import FoldBounds.FrontierNext

/-!
# Flat version-2 certificates

Version-2 analogue of `FrontierFlat`: an untrusted converter from a flat
list of `(direction word, leaf)` entries to an `NTree`.  Duplicate
entries, ancestor/descendant overlaps and omitted legal branches are
rejected exactly as in version 1.  `checkFlat_sound` follows from
`FrontierNext.check_sound` whatever the converter returns, because the
trie is re-validated by the verified checker.
-/
@[expose] public section

namespace FoldBounds
namespace FrontierNext

open LatticeModel Frontier FrontierCheck

/-- Kind of a flat version-2 frontier entry (any non-`node` tree). -/
abbrev FlatCert2 := List (List (Fin 4) × NTree)

/-- Entries whose relative word starts with direction `d`, with `d` removed. -/
def childEntries2 (d : Fin 4) (es : FlatCert2) : FlatCert2 :=
  es.filterMap (fun e => match e.1 with
    | d' :: rest => if d' = d then some (rest, e.2) else none
    | [] => none)

/-- Untrusted converter from flat entries (relative to prefix `p`) to a trie. -/
def buildTrie2 : ℕ → List (Fin 4) → FlatCert2 → Option NTree
  | 0, _, _ => none
  | fuel + 1, p, es =>
      match es with
      | [] => if (walkFrom p).Nodup then none else some .dead
      | [([], k)] => some k
      | _ =>
          if es.any (fun e => e.1.isEmpty) then none
          else
            match buildTrie2 fuel (p ++ [0]) (childEntries2 0 es),
                buildTrie2 fuel (p ++ [1]) (childEntries2 1 es),
                buildTrie2 fuel (p ++ [2]) (childEntries2 2 es),
                buildTrie2 fuel (p ++ [3]) (childEntries2 3 es) with
            | some a, some b, some c, some d => some (.node a b c d)
            | _, _, _, _ => none

/-- Checker for flat version-2 certificates. -/
def checkFlat2 (S : List Bool) (q : ℚ) (thr : ℕ) (es : FlatCert2) : Option FrontierResult :=
  match buildTrie2 (S.length + 1) [] es with
  | some t => check S q thr t
  | none => none

/-- **Soundness of the flat version-2 checker.** -/
theorem checkFlat2_sound (steps : ℕ) (seq : Fin (steps + 2) → Bool) (q : ℚ) (thr : ℕ)
    (es : FlatCert2) (r : FrontierResult) (h : checkFlat2 (List.ofFn seq) q thr es = some r) :
    Executable.massQ (weightQ steps q seq) Finset.univ ∈ Set.Icc r.zLo r.zHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) ∈
      Set.Icc r.eLo r.eHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ ∈
      Set.Icc r.nLo r.nHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) /
      Executable.massQ (weightQ steps q seq) Finset.univ ∈ Set.Icc r.pLo r.pHi := by
  unfold checkFlat2 at h
  split at h
  · exact check_sound steps seq q thr _ r h
  · exact absurd h (by simp)

end FrontierNext
end FoldBounds
