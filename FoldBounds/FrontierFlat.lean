module

public import FoldBounds.FrontierPartition

/-!
# Flat (Python-style) frontier certificates

The Python prototype represents a frontier as a flat list of nodes
`(path, kind)`.  `buildTrie` converts such a list (with paths written as
direction words) into an `FTree`, rejecting

* duplicate entries and ancestor/descendant overlaps (an entry with an
  empty relative suffix must be the *only* entry in its subtree),
* omitted branches: a missing child is filled with `dead` only when its
  prefix already collides; a missing *legal* child is rejected, exactly
  like Python's "certificate omits a nonempty subtree",
* entries deeper than the chain.

The converter is *not* trusted: `checkFlat_sound` follows from
`check_sound` whatever `buildTrie` returns, because the resulting trie is
re-validated by `check`.
-/
@[expose] public section

namespace FoldBounds
namespace FrontierCheck

open LatticeModel Frontier

/-- Kind of a flat frontier entry. -/
inductive NodeKind where
  | exact : NodeKind
  | bound (transfer : Bool) (witnesses : List (List (Fin 4))) : NodeKind
  | dead : NodeKind
  deriving Repr, DecidableEq

def NodeKind.toTree : NodeKind → FTree
  | .exact => .exact
  | .bound tr ws => .bound tr ws
  | .dead => .dead

/-- Flat certificate: list of (direction word after the fixed bond, kind). -/
abbrev FlatCert := List (List (Fin 4) × NodeKind)

/-- Entries whose relative word starts with direction `d`, with `d` removed. -/
def childEntries (d : Fin 4) (es : FlatCert) : FlatCert :=
  es.filterMap (fun e => match e.1 with
    | d' :: rest => if d' = d then some (rest, e.2) else none
    | [] => none)

/-- Untrusted converter from flat entries (relative to prefix `p`) to a trie. -/
def buildTrie : ℕ → List (Fin 4) → FlatCert → Option FTree
  | 0, _, _ => none
  | fuel + 1, p, es =>
      match es with
      | [] => if (walkFrom p).Nodup then none else some .dead
      | [([], k)] => some k.toTree
      | _ =>
          if es.any (fun e => e.1.isEmpty) then none
          else
            match buildTrie fuel (p ++ [0]) (childEntries 0 es),
                buildTrie fuel (p ++ [1]) (childEntries 1 es),
                buildTrie fuel (p ++ [2]) (childEntries 2 es),
                buildTrie fuel (p ++ [3]) (childEntries 3 es) with
            | some a, some b, some c, some d => some (.node a b c d)
            | _, _, _, _ => none

/-- Checker for flat certificates: build the trie, then run `check`. -/
def checkFlat (S : List Bool) (q : ℚ) (thr : ℕ) (es : FlatCert) : Option FrontierResult :=
  match buildTrie (S.length + 1) [] es with
  | some t => check S q thr t
  | none => none

/-- **Soundness of the flat checker** (inherits `check_sound`; the
converter needs no proof). -/
theorem checkFlat_sound (steps : ℕ) (seq : Fin (steps + 2) → Bool) (q : ℚ) (thr : ℕ)
    (es : FlatCert) (r : FrontierResult) (h : checkFlat (List.ofFn seq) q thr es = some r) :
    Executable.massQ (weightQ steps q seq) Finset.univ ∈ Set.Icc r.zLo r.zHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) ∈
      Set.Icc r.eLo r.eHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ ∈
      Set.Icc r.nLo r.nHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) /
      Executable.massQ (weightQ steps q seq) Finset.univ ∈ Set.Icc r.pLo r.pHi := by
  unfold checkFlat at h
  split at h
  · exact check_sound steps seq q thr _ r h
  · exact absurd h (by simp)

end FrontierCheck
end FoldBounds
