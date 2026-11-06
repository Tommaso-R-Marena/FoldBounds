module

public import FoldBounds.Frontier

/-!
# A verified weighted transfer relaxation for subtree partition functions

For a fixed prefix `p`, every genuine completion is in particular a
*relaxed walk*: a non-reversing continuation that never steps onto a
site of the fixed prefix (but may revisit its own new sites).  Contacts
between a new residue and a fixed prefix residue depend only on the
current site and index, so a forward dynamic program over states
`(site, last direction, cross score, next index)` computes the exact
cross-score-weighted sum over all relaxed walks *without enumerating
them*: duplicate states are merged (`mergeBy`: sort, then combine equal
neighbours) so the table size is
polynomial in the number of remaining steps.

Main results:

* `wsum_dp` — the merged forward DP equals the (exponential) backward
  path-sum specification, for an arbitrary finite-branching transition
  system (generic transfer-matrix soundness);
* `relaxed_sound` — genuine completions are dominated by the relaxed
  path-sum, for any nonnegative function of the realised cross score;
* `sum_ext_le_dp` — the corollary used by the frontier checker.
-/
@[expose] public section

namespace FoldBounds
namespace Transfer

open LatticeModel Frontier Finset

/-! ### Generic merged forward dynamic programming -/

section Generic

variable {σ : Type*} [DecidableEq σ]

/-- Weighted sum of a function over a weighted state list. -/
def wsum (F : σ → ℚ) (l : List (σ × ℚ)) : ℚ := (l.map (fun x => x.2 * F x.1)).sum

/-- Insert a weighted state, adding weights of equal states. -/
def insertAdd (x : σ × ℚ) : List (σ × ℚ) → List (σ × ℚ)
  | [] => [x]
  | y :: ys => if y.1 = x.1 then (y.1, y.2 + x.2) :: ys else y :: insertAdd x ys

/-- Merge duplicate states (weights add). -/
def merge (l : List (σ × ℚ)) : List (σ × ℚ) := l.foldr insertAdd []

omit [DecidableEq σ] in
theorem wsum_nil (F : σ → ℚ) : wsum F [] = 0 := rfl

omit [DecidableEq σ] in
theorem wsum_cons (F : σ → ℚ) (x : σ × ℚ) (l : List (σ × ℚ)) :
    wsum F (x :: l) = x.2 * F x.1 + wsum F l := by
  simp [wsum]

omit [DecidableEq σ] in
theorem wsum_append (F : σ → ℚ) (l l' : List (σ × ℚ)) :
    wsum F (l ++ l') = wsum F l + wsum F l' := by
  simp [wsum]

theorem wsum_insertAdd (F : σ → ℚ) (x : σ × ℚ) (l : List (σ × ℚ)) :
    wsum F (insertAdd x l) = x.2 * F x.1 + wsum F l := by
  induction l with
  | nil => simp [insertAdd, wsum]
  | cons y ys ih =>
      by_cases h : y.1 = x.1
      · simp only [insertAdd, h, if_true, wsum_cons]
        ring
      · simp only [insertAdd, h, if_false, wsum_cons, ih]
        ring

theorem wsum_merge (F : σ → ℚ) (l : List (σ × ℚ)) : wsum F (merge l) = wsum F l := by
  induction l with
  | nil => rfl
  | cons x l ih =>
      simp only [merge, List.foldr_cons] at ih ⊢
      rw [wsum_insertAdd, ih, wsum_cons]

/-- Combine adjacent entries with equal states (accumulator form). -/
def combineAux (acc : σ × ℚ) : List (σ × ℚ) → List (σ × ℚ)
  | [] => [acc]
  | y :: rest =>
      if acc.1 = y.1 then combineAux (acc.1, acc.2 + y.2) rest else acc :: combineAux y rest

def combineAdj : List (σ × ℚ) → List (σ × ℚ)
  | [] => []
  | x :: xs => combineAux x xs

section Sorting

variable {α : Type*}

/-- Fuel-bounded merge (structural, kernel-reducible). -/
def mergeF (le : α → α → Bool) : ℕ → List α → List α → List α
  | 0, xs, ys => xs ++ ys
  | _ + 1, [], ys => ys
  | _ + 1, x :: xs, [] => x :: xs
  | n + 1, x :: xs, y :: ys =>
      if le x y then x :: mergeF le n xs (y :: ys) else y :: mergeF le n (x :: xs) ys

/-- Split a list into alternate elements. -/
def splitAlt : List α → List α × List α
  | [] => ([], [])
  | [x] => ([x], [])
  | x :: y :: rest => (x :: (splitAlt rest).1, y :: (splitAlt rest).2)

/-- Fuel-bounded merge sort (structural, kernel-reducible). -/
def msortF (le : α → α → Bool) : ℕ → List α → List α
  | 0, l => l
  | _ + 1, [] => []
  | _ + 1, [x] => [x]
  | n + 1, x :: y :: rest =>
      let p := splitAlt (x :: y :: rest)
      mergeF le (rest.length + 2) (msortF le n p.1) (msortF le n p.2)

theorem mergeF_perm (le : α → α → Bool) :
    ∀ (n : ℕ) (xs ys : List α), (mergeF le n xs ys).Perm (xs ++ ys)
  | 0, xs, ys => List.Perm.refl _
  | _ + 1, [], ys => List.Perm.refl _
  | _ + 1, x :: xs, [] => by simp [mergeF]
  | n + 1, x :: xs, y :: ys => by
      unfold mergeF
      split_ifs
      · exact (mergeF_perm le n xs (y :: ys)).cons x
      · refine ((mergeF_perm le n (x :: xs) ys).cons y).trans ?_
        exact List.perm_middle.symm

theorem splitAlt_perm : ∀ l : List α, ((splitAlt l).1 ++ (splitAlt l).2).Perm l
  | [] => List.Perm.refl _
  | [x] => List.Perm.refl _
  | x :: y :: rest => by
      simp only [splitAlt, List.cons_append]
      refine List.Perm.cons x ?_
      refine List.perm_middle.trans ?_
      exact (splitAlt_perm rest).cons y

theorem msortF_perm (le : α → α → Bool) :
    ∀ (n : ℕ) (l : List α), (msortF le n l).Perm l
  | 0, l => List.Perm.refl _
  | _ + 1, [] => List.Perm.refl _
  | _ + 1, [x] => List.Perm.refl _
  | n + 1, x :: y :: rest => by
      simp only [msortF]
      refine (mergeF_perm le _ _ _).trans ?_
      refine ((msortF_perm le n _).append (msortF_perm le n _)).trans ?_
      exact splitAlt_perm _

end Sorting

/-- Sort by an arbitrary comparator, then combine adjacent equal states.
Soundness never depends on `le`; only the amount of merging does. -/
def mergeBy (le : σ × ℚ → σ × ℚ → Bool) (l : List (σ × ℚ)) : List (σ × ℚ) :=
  combineAdj (msortF le l.length l)

theorem wsum_combineAux (F : σ → ℚ) (acc : σ × ℚ) (l : List (σ × ℚ)) :
    wsum F (combineAux acc l) = acc.2 * F acc.1 + wsum F l := by
  induction l generalizing acc with
  | nil => simp [combineAux, wsum]
  | cons y rest ih =>
      by_cases h : acc.1 = y.1
      · simp only [combineAux, h, if_true, ih, wsum_cons]
        ring
      · simp only [combineAux, h, if_false, wsum_cons, ih]

theorem wsum_mergeBy (F : σ → ℚ) (le : σ × ℚ → σ × ℚ → Bool) (l : List (σ × ℚ)) :
    wsum F (mergeBy le l) = wsum F l := by
  have hperm : wsum F (msortF le l.length l) = wsum F l := by
    unfold wsum
    exact ((msortF_perm le l.length l).map (fun x => x.2 * F x.1)).sum_eq
  rw [mergeBy, ← hperm]
  cases h : msortF le l.length l with
  | nil => rfl
  | cons x xs => simp only [combineAdj, wsum_combineAux, wsum_cons]

/-- One forward layer: every state is replaced by its successors. -/
def expand (next : σ → List σ) (l : List (σ × ℚ)) : List (σ × ℚ) :=
  l.flatMap (fun x => (next x.1).map (fun s => (s, x.2)))

/-- Backward path-sum specification (exponential if evaluated naively). -/
def spec (next : σ → List σ) (G : σ → ℚ) : ℕ → σ → ℚ
  | 0, s => G s
  | r + 1, s => ((next s).map (spec next G r)).sum

omit [DecidableEq σ] in
theorem wsum_expand (next : σ → List σ) (G : σ → ℚ) (r : ℕ) (l : List (σ × ℚ)) :
    wsum (spec next G r) (expand next l) = wsum (spec next G (r + 1)) l := by
  induction l with
  | nil => rfl
  | cons x l ih =>
      simp only [expand, List.flatMap_cons, wsum_append] at ih ⊢
      rw [ih, wsum_cons, spec]
      congr 1
      simp only [wsum, List.map_map, Function.comp_def]
      rw [List.sum_map_mul_left]

/-- The merged forward DP, run for `r` layers (states merged after sorting
with the comparator `le`). -/
def dp (next : σ → List σ) (le : σ × ℚ → σ × ℚ → Bool) : ℕ → List (σ × ℚ) → List (σ × ℚ)
  | 0, l => l
  | r + 1, l => dp next le r (mergeBy le (expand next l))

/-- **Transfer-matrix soundness.** The merged forward DP computes exactly
the backward path-sum specification, for any comparator. -/
theorem wsum_dp (next : σ → List σ) (le : σ × ℚ → σ × ℚ → Bool) (G : σ → ℚ) (r : ℕ)
    (l : List (σ × ℚ)) :
    wsum G (dp next le r l) = wsum (spec next G r) l := by
  induction r generalizing l with
  | zero => rfl
  | succ r ih =>
      rw [dp, ih, wsum_mergeBy, wsum_expand]

omit [DecidableEq σ] in
theorem spec_nonneg (next : σ → List σ) (G : σ → ℚ) (hG : ∀ s, 0 ≤ G s) :
    ∀ r s, 0 ≤ spec next G r s := by
  intro r
  induction r with
  | zero => exact hG
  | succ r ih =>
      intro s
      rw [spec]
      exact List.sum_nonneg (by
        intro x hx
        obtain ⟨y, _, rfl⟩ := List.mem_map.1 hx
        exact ih y)

end Generic

/-! ### The lattice relaxation -/

/-- State: current site, last direction, accumulated cross score, index
of the next residue to place. -/
abbrev RState := Site × Fin 4 × ℕ × ℕ

/-- All four directions, in constructor order. -/
def allDirs : List (Fin 4) := [0, 1, 2, 3]

theorem sum_allDirs {M : Type*} [AddCommMonoid M] (f : Fin 4 → M) :
    (allDirs.map f).sum = ∑ d : Fin 4, f d := by
  simp [allDirs, Fin.sum_univ_four, add_assoc]

/-- Cross contacts gained by placing residue `j` at site `s`, against the
fixed prefix walk `Wp` with `m` placed sites. -/
def crossAt (S : List Bool) (Wp : List Site) (m : ℕ) (j : ℕ) (s : Site) : ℕ :=
  ∑ i ∈ range S.length,
    if candidatePair S i j && decide (i < m) && manhattan1 (siteAt Wp i) s then 1 else 0

/-- Admissible relaxed step: no immediate reversal and no fixed prefix site. -/
def allowed (Wp : List Site) (st : RState) (e : Fin 4) : Bool :=
  decide (e ≠ opp st.2.1) && !(Wp.contains (advance st.1 e))

/-- Relaxed successor states. -/
def relaxedNext (S : List Bool) (Wp : List Site) (m : ℕ) (st : RState) : List RState :=
  (allDirs.filter (allowed Wp st)).map (fun e =>
    (advance st.1 e, e, st.2.2.1 + crossAt S Wp m st.2.2.2 (advance st.1 e), st.2.2.2 + 1))

/-- Relaxed state reached by a (genuine) prefix `p'` extending `p`. -/
def stateOf (S : List Bool) (p p' : List (Fin 4)) : RState :=
  (endFrom first p', lastDir p', crossContacts S (p.length + 2) (walkFrom p'), p'.length + 2)

theorem sum_map_filter_eq {α : Type*} (l : List α) (P : α → Bool) (f : α → ℚ) :
    ((l.filter P).map f).sum = (l.map (fun x => if P x then f x else 0)).sum := by
  induction l with
  | nil => rfl
  | cons x l ih =>
      by_cases h : P x
      · simp [h, ih]
      · simp [h, ih]

/-- Appending one site adds exactly the cross contacts of the new residue. -/
theorem crossContacts_snoc (S : List Bool) (p : List (Fin 4)) (W : List Site) (s : Site)
    (hpre : walkFrom p <+: W) (hlt : W.length < S.length) :
    crossContacts S (p.length + 2) (W ++ [s]) =
      crossContacts S (p.length + 2) W +
        crossAt S (walkFrom p) (p.length + 2) W.length s := by
  have hm : (walkFrom p).length = p.length + 2 := walkFrom_length p
  have hle : p.length + 2 ≤ W.length := by have := hpre.length_le; omega
  rw [crossContacts, crossContacts, pairSum_split S (Q := fun i j => candidatePair S i j &&
      decide (i < p.length + 2) && decide (p.length + 2 ≤ j) && decide (j < W.length) &&
      manhattan1 (siteAt W i) (siteAt W j))
    (R := fun i j => candidatePair S i j && decide (i < p.length + 2) && decide (j = W.length) &&
      manhattan1 (siteAt (walkFrom p) i) s)]
  · congr 1
    unfold pairSum crossAt
    refine Finset.sum_congr rfl (fun i _ => ?_)
    have : ∀ j, (if (candidatePair S i j && decide (i < p.length + 2) && decide (j = W.length) &&
        manhattan1 (siteAt (walkFrom p) i) s) = true then 1 else 0) =
        if j = W.length then (if (candidatePair S i W.length && decide (i < p.length + 2) &&
          manhattan1 (siteAt (walkFrom p) i) s) = true then 1 else 0) else 0 := by
      intro j
      by_cases hj : j = W.length
      · subst hj; simp
      · simp [hj]
    simp only [this]
    rw [Finset.sum_ite_eq' (range S.length) W.length]
    simp [hlt]
  · intro i j _ _
    by_cases hc : candidatePair S i j = true
    · have hij := candidatePair_lt hc
      by_cases hi : i < p.length + 2
      · by_cases hjW : j = W.length
        · subst hjW
          have h1 : siteAt (W ++ [s]) W.length = s := by simp [siteAt]
          have h2 : siteAt (W ++ [s]) i = siteAt (walkFrom p) i :=
            (siteAt_prefix (hpre.trans (List.prefix_append W [s])) (by omega)).symm
          rw [h1, h2]
          cases hb : manhattan1 (siteAt (walkFrom p) i) s <;> simp [hc, hi]
          omega
        · by_cases hjlt : j < W.length
          · rw [siteAt_prefix (List.prefix_append W [s]) (show i < W.length by omega),
              siteAt_prefix (List.prefix_append W [s]) hjlt]
            cases hb : manhattan1 (siteAt W i) (siteAt W j) <;> simp [hc, hi, hjW, hjlt] <;> omega
          · have : ¬ j < W.length + 1 := by omega
            simp [hjW, hjlt, this]
      · simp [hi]
    · simp [hc]

theorem crossContacts_self (S : List Bool) (p : List (Fin 4)) :
    crossContacts S (p.length + 2) (walkFrom p) = 0 := by
  have hm : (walkFrom p).length = p.length + 2 := walkFrom_length p
  unfold crossContacts pairSum
  refine Finset.sum_eq_zero (fun i _ => Finset.sum_eq_zero (fun j _ => ?_))
  rw [hm]
  have : ¬ ((p.length + 2 ≤ j) ∧ j < p.length + 2) := by omega
  simp only [Bool.and_eq_true, decide_eq_true_eq, ite_eq_right_iff, one_ne_zero, imp_false]
  tauto

/-- **Relaxation soundness.** For every genuine prefix `p'` extending `p`,
the sum over genuine completions of any nonnegative function of the
realised cross score is bounded by the relaxed path-sum from the state of
`p'`. -/
theorem relaxed_sound (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2)
    (p : List (Fin 4)) (g : ℕ → ℚ) (hg : ∀ k, 0 ≤ g k) :
    ∀ r (p' : List (Fin 4)), p <+: p' → (walkFrom p').Nodup → p'.length ≤ steps →
      steps - p'.length = r →
      ∑ l ∈ ext steps p', g (crossContacts S (p.length + 2) (walkFrom l)) ≤
        spec (relaxedNext S (walkFrom p) (p.length + 2)) (fun st => g st.2.2.1) r
          (stateOf S p p') := by
  intro r
  induction r with
  | zero =>
      intro p' _ hn hle hr
      rw [ext_eq_singleton (by omega) hn, Finset.sum_singleton]
      rfl
  | succ r ih =>
      intro p' hpp' hn hle hr
      have hlt : p'.length < steps := by omega
      rw [sum_ext_split hlt, spec, relaxedNext, List.map_map, sum_map_filter_eq, sum_allDirs]
      refine Finset.sum_le_sum (fun d _ => ?_)
      by_cases hnd : (walkFrom (p' ++ [d])).Nodup
      · have hrev : d ≠ opp (lastDir p') := by
          rintro rfl; exact not_nodup_reverse p' hnd
        have hnew : advance (endFrom first p') d ∉ walkFrom p' := by
          intro hm
          rw [walkFrom_append_single, List.nodup_append] at hnd
          exact hnd.2.2 _ hm _ (List.mem_singleton_self _) rfl
        have hpre : walkFrom p <+: walkFrom p' := walkFrom_prefix hpp'
        have hnotp : advance (endFrom first p') d ∉ walkFrom p :=
          fun hm => hnew (hpre.subset hm)
        have hall : allowed (walkFrom p) (stateOf S p p') d = true := by
          have hc : (walkFrom p).contains (advance (endFrom first p') d) = false := by
            simpa using hnotp
          simp only [allowed, stateOf, Bool.and_eq_true, decide_eq_true_eq, Bool.not_eq_true', hc]
          exact ⟨hrev, trivial⟩
        rw [if_pos hall]
        have hstate : stateOf S p (p' ++ [d]) =
            (advance (stateOf S p p').1 d, d,
              (stateOf S p p').2.2.1 + crossAt S (walkFrom p) (p.length + 2)
                (stateOf S p p').2.2.2 (advance (stateOf S p p').1 d),
              (stateOf S p p').2.2.2 + 1) := by
          simp only [stateOf, endFrom_append, endFrom, lastDir_append_single,
            List.length_append, List.length_singleton]
          refine Prod.ext rfl (Prod.ext rfl (Prod.ext ?_ (by ring)))
          simp only
          rw [walkFrom_append_single, crossContacts_snoc S p _ _ hpre
            (by rw [walkFrom_length]; omega), walkFrom_length]
        have := ih (p' ++ [d]) (hpp'.trans (List.prefix_append _ _)) hnd
          (by simp; omega) (by simp; omega)
        rw [hstate] at this
        exact this
      · rw [ext_eq_empty_of_not_nodup hnd, Finset.sum_empty]
        split_ifs
        · exact spec_nonneg _ _ (fun st => hg _) _ _
        · exact le_refl _

/-- Initial relaxed table for prefix `p`. -/
def initTable (S : List Bool) (p : List (Fin 4)) : List (RState × ℚ) :=
  [(stateOf S p p, 1)]

/-- Lexicographic comparator on relaxed states (efficiency only). -/
def rstateLe (a b : RState × ℚ) : Bool :=
  let x := a.1
  let y := b.1
  decide (x.1.1 < y.1.1) || (decide (x.1.1 = y.1.1) && (decide (x.1.2 < y.1.2) ||
    (decide (x.1.2 = y.1.2) && (decide (x.2.1.val < y.2.1.val) ||
      (decide (x.2.1 = y.2.1) && (decide (x.2.2.1 < y.2.2.1) ||
        (decide (x.2.2.1 = y.2.2.1) && decide (x.2.2.2 ≤ y.2.2.2))))))))

/-- The final relaxed table: merged forward DP over the remaining steps. -/
def relaxedTable (steps : ℕ) (S : List Bool) (p : List (Fin 4)) : List (RState × ℚ) :=
  dp (relaxedNext S (walkFrom p) (p.length + 2)) rstateLe (steps - p.length) (initTable S p)

/-- **Corollary used by the checker.** The DP table bounds the sum over all
genuine completions of any nonnegative function of the cross score. -/
theorem sum_ext_le_dp (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2)
    (p : List (Fin 4)) (hn : (walkFrom p).Nodup) (hle : p.length ≤ steps)
    (g : ℕ → ℚ) (hg : ∀ k, 0 ≤ g k) :
    ∑ l ∈ ext steps p, g (crossContacts S (p.length + 2) (walkFrom l)) ≤
      wsum (fun st => g st.2.2.1) (relaxedTable steps S p) := by
  rw [relaxedTable, wsum_dp, initTable, wsum_cons, wsum_nil, one_mul, add_zero]
  exact relaxed_sound steps S hS p g hg _ p (List.prefix_refl p) hn hle rfl

end Transfer
end FoldBounds
