module

public import FoldBounds.MotionGNMAdaptiveCases

/-!
# Heterogeneous-depth, proof-carrying contact partitions

For any finite Boolean decision tree, *different leaves may fix different
numbers of contact switches*. The tree covers every full contact realization
without requiring all leaves to live at the same depth. A leaf provides pair
score intervals; if the intervals soundly enclose the true model scores,
existing signed-concordance inequalities give a rigorous global lower bound.

This theorem is independent of the protein coordinates and numeric certificate
bytes. The verified Python replay supplies finite ground-graph calculations.
-/

@[expose] public section
namespace FoldBounds.MotionGNMPrefixTree

open FoldBounds.MotionRobustConcordance

universe u

inductive SplitTree (α : Type u) where
  | leaf : α → SplitTree α
  | branch : SplitTree α → SplitTree α → SplitTree α

/-- Select the unique terminal leaf by reading subsequent contact bits. -/
def pickLeaf {α : Type u} : SplitTree α → (ℕ → Bool) → ℕ → α
  | .leaf a, _, _ => a
  | .branch left right, world, depth =>
      if world depth then pickLeaf right world (depth + 1)
      else pickLeaf left world (depth + 1)

/-- A property checked at every leaf of a possibly unbalanced tree. -/
def AllLeaves {α : Type u} (P : α → Prop) : SplitTree α → Prop
  | .leaf a => P a
  | .branch left right => AllLeaves P left ∧ AllLeaves P right

/-- No full Boolean contact assignment is omitted by an unbalanced tree. -/
theorem pickLeaf_sound {α : Type u} (P : α → Prop) :
    ∀ (tree : SplitTree α) (world : ℕ → Bool) (depth : ℕ),
      AllLeaves P tree → P (pickLeaf tree world depth)
  | .leaf a, _, _, h => h
  | .branch left right, world, depth, h => by
      change AllLeaves P left ∧ AllLeaves P right at h
      cases hb : world depth with
      | true =>
          simpa [pickLeaf, hb] using
            pickLeaf_sound P right world (depth + 1) h.2
      | false =>
          simpa [pickLeaf, hb] using
            pickLeaf_sound P left world (depth + 1) h.1

/-- Number of terminal certificate leaves (possibly at different depths). -/
def leafCount {α : Type u} : SplitTree α → ℕ
  | .leaf _ => 1
  | .branch left right => leafCount left + leafCount right

/-- Number of internal split nodes. -/
def splitCount {α : Type u} : SplitTree α → ℕ
  | .leaf _ => 0
  | .branch left right => 1 + splitCount left + splitCount right

/-- Each split replaces one leaf by two; hence a full binary prefix tree
has exactly L-1 internal nodes, independently of its varying depths. -/
theorem leafCount_eq_splitCount_add_one {α : Type u} (tree : SplitTree α) :
    leafCount tree = splitCount tree + 1 := by
  induction tree with
  | leaf _ => simp [leafCount, splitCount]
  | branch left right ihLeft ihRight =>
      simp only [leafCount, splitCount]
      omega

/-- Sound score envelopes for two oriented members of each comparison. -/
structure SignedEnvelopes (n : ℕ) where
  aLo : Fin n → ℚ
  aHi : Fin n → ℚ
  bLo : Fin n → ℚ
  bHi : Fin n → ℚ

/-- A pessimistic signed vote: +1 only for forced correct strict ordering. -/
def lowerVote (aHi bLo : ℚ) : ℚ :=
  if aHi < bLo then 1 else -1

/-- Sum the pessimistic votes of all comparisons in one leaf. -/
def leafLower {n : ℕ} (e : SignedEnvelopes n) : ℚ :=
  ∑ i : Fin n, lowerVote (e.aHi i) (e.bLo i)

/-- A refinement can only improve a pessimistic signed vote when it lowers
one upper score endpoint and raises the other lower score endpoint. -/
theorem lowerVote_mono (aHi oldLo newHi newLo : ℚ)
    (hHi : newHi ≤ aHi) (hLo : oldLo ≤ newLo) :
    lowerVote aHi oldLo ≤ lowerVote newHi newLo := by
  by_cases h : aHi < oldLo
  · have h' : newHi < newLo :=
      lt_of_lt_of_le (lt_of_le_of_lt hHi h) hLo
    simp [lowerVote, h, h']
  · by_cases h' : newHi < newLo
    · simp [lowerVote, h, h']
    · simp [lowerVote, h, h']

/-- Any interval refinement, including intersecting old and new approximate
certificates, cannot reduce the certified signed-concordance lower bound. -/
theorem leafLower_mono {n : ℕ} (old new : SignedEnvelopes n)
    (hHi : ∀ i, new.aHi i ≤ old.aHi i)
    (hLo : ∀ i, old.bLo i ≤ new.bLo i) :
    leafLower old ≤ leafLower new := by
  unfold leafLower
  apply Finset.sum_le_sum
  intro i _
  exact lowerVote_mono (old.aHi i) (old.bLo i)
    (new.aHi i) (new.bLo i) (hHi i) (hLo i)

/-- Intersecting two valid rational score certificates is valid and no wider. -/
theorem intersect_valid (lo hi lo' hi' trueValue : ℚ)
    (h : lo ≤ trueValue ∧ trueValue ≤ hi)
    (h' : lo' ≤ trueValue ∧ trueValue ≤ hi') :
    max lo lo' ≤ trueValue ∧ trueValue ≤ min hi hi' := by
  exact ⟨max_le h.1 h'.1, le_min h.2 h'.2⟩

/-- General soundness of a variable-depth protein rank concordance tree.
Worlds include arbitrary *jointly correlated* Boolean contact choices.
For each world, source-verified graph bounds must enclose the same world's
true two scores for every oriented comparison. All leaves need only certify
one common pessimistic lower count. The tree may stop at any depth. -/
theorem prefixTree_concordance_sound {n : ℕ}
    (tree : SplitTree (SignedEnvelopes n))
    (scoreA scoreB : (ℕ → Bool) → Fin n → ℚ)
    (guaranteed : ℚ)
    (hA : ∀ world i,
      (pickLeaf tree world 0).aLo i ≤ scoreA world i ∧
      scoreA world i ≤ (pickLeaf tree world 0).aHi i)
    (hB : ∀ world i,
      (pickLeaf tree world 0).bLo i ≤ scoreB world i ∧
      scoreB world i ≤ (pickLeaf tree world 0).bHi i)
    (hLeaf : AllLeaves (fun e => guaranteed ≤ leafLower e) tree) :
    ∀ world : ℕ → Bool,
      guaranteed ≤ ∑ i : Fin n,
        alignedSign (scoreA world i) (scoreB world i) := by
  intro world
  let e := pickLeaf tree world 0
  have hc : guaranteed ≤ leafLower e :=
    pickLeaf_sound (fun r => guaranteed ≤ leafLower r) tree world 0 hLeaf
  have hpair := sum_alignedSign_interval
    e.aLo e.aHi e.bLo e.bHi (scoreA world) (scoreB world)
    (fun i => hA world i) (fun i => hB world i)
  have heq : leafLower e =
      ∑ i : Fin n, if e.aHi i < e.bLo i then (1 : ℚ) else -1 := rfl
  rw [heq] at hc
  exact le_trans hc hpair.1

end FoldBounds.MotionGNMPrefixTree