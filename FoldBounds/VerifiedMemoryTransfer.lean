module

public import FoldBounds.FrontierCheck
public import FoldBounds.VerifiedGeometry

/-!
# A verified finite-memory transfer relaxation

The transfer relaxation of `FoldBounds.Transfer` forbids only the sites of
the fixed prefix (and immediate reversal).  Here a relaxed walk also
remembers the last `k` sites it placed itself (with their chain indices)
and may not step onto them.  Every genuine self-avoiding completion
avoids *all* earlier sites, so it is in particular a `k`-memory relaxed
walk (`memCover`).

The retained history is also used in both logical directions for the
unplaced–unplaced contacts of the new residue `j` at site `s`:

* **observed contacts (lower bound, complement direction).**  A retained
  residue `i` (`i ≥ m`, candidate pair) whose remembered site is
  adjacent to `s` is a genuine contact: `loAt`.
* **excluded contacts (upper bound, event direction).**  A retained
  residue whose remembered site is *not* adjacent to `s` is genuinely
  not a contact.  Hence the contacts of `j` with earlier unplaced
  residues are at most `loAt + outAt` (`outAt` counts the candidate
  partners that are not retained), and at most `slots j` by the
  degree lemma: `upAt = min (slots j) (loAt + outAt)`.

Main results:

* `memState_inv` — invariants of the deterministic state of a genuine
  partial conformation (endpoint, direction, exact cross score, genuine
  history, `ffLo ≤ ff ≤ ffHi`);
* `memCover` — the sum of any nonnegative state function over genuine
  completions is at most its sum over the memory-`k` DP table (computed by
  the merged forward DP, exact by `Transfer.wsum_dp`);
* `memContacts_bounds` — `known + cross + ffLo ≤ k(c) ≤ known + cross + ffHi`
  and `k(c) ≤ geomCap N (prefix ++ endpoint :: history)` for every
  genuine completion `c` with final state `st`.
-/
@[expose] public section

namespace FoldBounds
namespace MemTransfer

open LatticeModel Frontier Transfer Geometry Finset

/-! ### States and transitions -/

/-- A finite-memory relaxed state. -/
structure MState where
  /-- Current endpoint. -/
  cur : Site
  /-- Last direction. -/
  dir : Fin 4
  /-- Retained `(site, index)` of earlier residues, newest first (at most `k`). -/
  hist : List (Site × ℕ)
  /-- Exact accumulated placed–unplaced (cross) contacts. -/
  cross : ℕ
  /-- Observed (certain) unplaced–unplaced contacts. -/
  ffLo : ℕ
  /-- Admissible upper bound on unplaced–unplaced contacts so far. -/
  ffHi : ℕ
  /-- Index of the next residue to place. -/
  idx : ℕ
  deriving DecidableEq, Repr

/-- First retained entry with index `i`. -/
def lookupH (h : List (Site × ℕ)) (i : ℕ) : Option (Site × ℕ) := h.find? (fun x => x.2 == i)

/-- Observed contacts of residue `j` at site `s` with retained unplaced residues. -/
def loAt (S : List Bool) (m : ℕ) (h : List (Site × ℕ)) (j : ℕ) (s : Site) : ℕ :=
  ∑ i ∈ range S.length, if (candidatePair S i j && decide (m ≤ i) &&
    ((lookupH h i).map (fun x => manhattan1 x.1 s)).getD false) = true then 1 else 0

/-- Candidate unplaced partners of residue `j` that are not retained. -/
def outAt (S : List Bool) (m : ℕ) (h : List (Site × ℕ)) (j : ℕ) : ℕ :=
  ∑ i ∈ range S.length, if (candidatePair S i j && decide (m ≤ i) &&
    (lookupH h i).isNone) = true then 1 else 0

/-- Admissible upper bound on the unplaced–unplaced contacts of residue `j`. -/
def upAt (S : List Bool) (m : ℕ) (h : List (Site × ℕ)) (j : ℕ) (s : Site) : ℕ :=
  min (slots S j) (loAt S m h j s + outAt S m h j)

/-- Admissible relaxed step: no reversal, no prefix site, no retained site. -/
def memAllowed (Wp : List Site) (st : MState) (e : Fin 4) : Bool :=
  decide (e ≠ opp st.dir) && !(Wp.contains (advance st.cur e)) &&
    !(st.hist.any (fun x => x.1 == advance st.cur e))

/-- Deterministic state update when stepping in direction `e`. -/
def memStep (S : List Bool) (Wp : List Site) (m k : ℕ) (st : MState) (e : Fin 4) : MState :=
  let s := advance st.cur e
  { cur := s
    dir := e
    hist := ((st.cur, st.idx - 1) :: st.hist).take k
    cross := st.cross + crossAt S Wp m st.idx s
    ffLo := st.ffLo + loAt S m st.hist st.idx s
    ffHi := st.ffHi + upAt S m st.hist st.idx s
    idx := st.idx + 1 }

/-- Relaxed successors. -/
def memNext (S : List Bool) (Wp : List Site) (m k : ℕ) (st : MState) : List MState :=
  (allDirs.filter (memAllowed Wp st)).map (memStep S Wp m k st)

/-- Initial state at prefix `p`. -/
def memInit (p : List (Fin 4)) : MState :=
  ⟨endFrom first p, lastDir p, [], 0, 0, 0, p.length + 2⟩

/-- Deterministic state of a genuine extension `p'` of `p`. -/
def gstate (S : List Bool) (k : ℕ) (p p' : List (Fin 4)) : MState :=
  (p'.drop p.length).foldl (memStep S (walkFrom p) (p.length + 2) k) (memInit p)

/-- Lexicographic comparator (efficiency only; soundness never depends on it). -/
def mstateKey (st : MState) : List ℤ :=
  [st.cur.1, st.cur.2, st.dir.val, st.cross, st.ffLo, st.ffHi, st.idx, st.hist.length] ++
    st.hist.flatMap (fun x => [x.1.1, x.1.2, (x.2 : ℤ)])

def lexLe : List ℤ → List ℤ → Bool
  | [], _ => true
  | _ :: _, [] => false
  | a :: as, b :: bs => decide (a < b) || (decide (a = b) && lexLe as bs)

def mstateLe (a b : MState × ℚ) : Bool := lexLe (mstateKey a.1) (mstateKey b.1)

/-! ### A stack-safe merged forward DP

`Transfer.dp` merges duplicate states with a non-tail-recursive merge sort,
which overflows the native stack on the large tables produced by memory
states.  `dpM` uses a bottom-up merge sort and a duplicate combiner that
are tail recursive *and* structurally recursive (so they also reduce in
the kernel).  Only the permutation property of the sort matters for
soundness (`wsum_dpM`); the comparator affects efficiency only. -/

section StackSafe

variable {α σ : Type*}

/-- Tail-recursive merge with fuel and accumulator. -/
def mergeAcc (le : α → α → Bool) : ℕ → List α → List α → List α → List α
  | 0, acc, a, b => acc.reverseAux (a ++ b)
  | _ + 1, acc, [], b => acc.reverseAux b
  | _ + 1, acc, a, [] => acc.reverseAux a
  | f + 1, acc, x :: xs, y :: ys =>
      if le x y then mergeAcc le f (x :: acc) xs (y :: ys)
      else mergeAcc le f (y :: acc) (x :: xs) ys

theorem mergeAcc_perm [DecidableEq α] (le : α → α → Bool) :
    ∀ f acc a b, (mergeAcc le f acc a b).Perm (acc ++ a ++ b) := by
  intro f
  induction f with
  | zero =>
      intro acc a b
      simp only [mergeAcc, List.reverseAux_eq, List.append_assoc]
      exact (List.reverse_perm acc).append_right _
  | succ f ih =>
      intro acc a b
      cases a with
      | nil =>
          simp only [mergeAcc, List.reverseAux_eq, List.append_nil]
          exact (List.reverse_perm acc).append_right _
      | cons x xs =>
          cases b with
          | nil =>
              simp only [mergeAcc, List.reverseAux_eq, List.append_nil]
              exact (List.reverse_perm acc).append_right _
          | cons y ys =>
              simp only [mergeAcc]
              split_ifs
              · refine (ih _ _ _).trans ?_
                rw [List.perm_iff_count]; intro z
                simp only [List.count_append, List.count_cons]; omega
              · refine (ih _ _ _).trans ?_
                rw [List.perm_iff_count]; intro z
                simp only [List.count_append, List.count_cons]; omega

/-- Merge adjacent runs (tail recursive). -/
def mergePairs (le : α → α → Bool) : List (List α) → List (List α) → List (List α)
  | acc, a :: b :: rest => mergePairs le (mergeAcc le (a.length + b.length) [] a b :: acc) rest
  | acc, [a] => a :: acc
  | acc, [] => acc

theorem mergePairs_perm [DecidableEq α] (le : α → α → Bool) :
    ∀ runs acc, (mergePairs le acc runs).flatten.Perm (acc.flatten ++ runs.flatten) := by
  intro runs
  induction runs using List.twoStepInduction with
  | nil => intro acc; simp [mergePairs]
  | singleton a =>
      intro acc
      simp only [mergePairs, List.flatten_cons, List.flatten_nil, List.append_nil]
      exact List.perm_append_comm
  | cons_cons a b rest ih _ =>
      intro acc
      simp only [mergePairs]
      refine (ih _).trans ?_
      have hm := List.perm_iff_count.1 (mergeAcc_perm le (a.length + b.length) [] a b)
      rw [List.perm_iff_count]; intro z
      have := hm z
      simp only [List.flatten_cons, List.count_append, List.nil_append] at this ⊢
      omega

/-- Repeated pairwise merging with fuel. -/
def sortRuns (le : α → α → Bool) : ℕ → List (List α) → List α
  | 0, runs => runs.flatten
  | _ + 1, [] => []
  | _ + 1, [r] => r
  | f + 1, r :: s :: rest => sortRuns le f (mergePairs le [] (r :: s :: rest))

theorem sortRuns_perm [DecidableEq α] (le : α → α → Bool) :
    ∀ f runs, (sortRuns le f runs).Perm runs.flatten := by
  intro f
  induction f with
  | zero => intro runs; exact List.Perm.refl _
  | succ f ih =>
      intro runs
      match runs with
      | [] => simp [sortRuns]
      | [r] => simp [sortRuns]
      | r :: s :: rest =>
          simp only [sortRuns]
          refine (ih _).trans ?_
          have := mergePairs_perm le (r :: s :: rest) []
          simpa using this

/-- Stack-safe bottom-up merge sort (a permutation of its input). -/
def tsort (le : α → α → Bool) (l : List α) : List α := sortRuns le l.length (l.map (fun x => [x]))

theorem tsort_perm [DecidableEq α] (le : α → α → Bool) (l : List α) : (tsort le l).Perm l := by
  unfold tsort
  refine (sortRuns_perm le _ _).trans ?_
  have : ∀ l : List α, (l.map (fun x => [x])).flatten = l := by
    intro l; induction l with
    | nil => rfl
    | cons x l ih => simp [ih]
  rw [this]

variable [DecidableEq σ]

/-- Tail-recursive combination of adjacent equal states (accumulator reversed). -/
def combineTR : List (σ × ℚ) → List (σ × ℚ) → List (σ × ℚ)
  | acc, [] => acc
  | [], x :: xs => combineTR [x] xs
  | a :: acc, x :: xs =>
      if a.1 = x.1 then combineTR ((a.1, a.2 + x.2) :: acc) xs else combineTR (x :: a :: acc) xs

theorem wsum_combineTR (F : σ → ℚ) :
    ∀ l acc, wsum F (combineTR acc l) = wsum F acc + wsum F l := by
  intro l
  induction l with
  | nil => intro acc; simp [combineTR, wsum_nil]
  | cons x xs ih =>
      intro acc
      cases acc with
      | nil => simp only [combineTR]; rw [ih, wsum_cons, wsum_cons, wsum_nil]; ring
      | cons a acc =>
          simp only [combineTR]
          split_ifs with h
          · rw [ih, wsum_cons, wsum_cons, wsum_cons, h]; ring
          · rw [ih, wsum_cons, wsum_cons, wsum_cons]; ring

omit [DecidableEq σ] in
theorem wsum_perm (F : σ → ℚ) {l l' : List (σ × ℚ)} (h : l.Perm l') : wsum F l = wsum F l' := by
  unfold wsum; exact (h.map _).sum_eq

/-- Stack-safe merge of duplicate states. -/
def mergeTR (le : σ × ℚ → σ × ℚ → Bool) (l : List (σ × ℚ)) : List (σ × ℚ) :=
  combineTR [] (tsort le l)

theorem wsum_mergeTR (F : σ → ℚ) (le : σ × ℚ → σ × ℚ → Bool) (l : List (σ × ℚ)) :
    wsum F (mergeTR le l) = wsum F l := by
  rw [mergeTR, wsum_combineTR, wsum_nil, zero_add, wsum_perm F (tsort_perm le l)]

/-- The stack-safe merged forward DP. -/
def dpM (next : σ → List σ) (le : σ × ℚ → σ × ℚ → Bool) : ℕ → List (σ × ℚ) → List (σ × ℚ)
  | 0, l => l
  | r + 1, l => dpM next le r (mergeTR le (expand next l))

/-- **Transfer-matrix soundness (stack-safe DP).** -/
theorem wsum_dpM (next : σ → List σ) (le : σ × ℚ → σ × ℚ → Bool) (G : σ → ℚ) (r : ℕ)
    (l : List (σ × ℚ)) : wsum G (dpM next le r l) = wsum (spec next G r) l := by
  induction r generalizing l with
  | zero => rfl
  | succ r ih => rw [dpM, ih, wsum_mergeTR, wsum_expand]

end StackSafe

/-- The memory-`k` table for prefix `p`. -/
def memTable (steps : ℕ) (S : List Bool) (k : ℕ) (p : List (Fin 4)) : List (MState × ℚ) :=
  dpM (memNext S (walkFrom p) (p.length + 2) k) mstateLe (steps - p.length) [(memInit p, 1)]

/-! ### Unplaced–unplaced contacts, residue by residue -/

/-- Contacts of residue `j` with earlier unplaced residues (`i ≥ m`). -/
def ffj (S : List Bool) (m : ℕ) (W : List Site) (j : ℕ) : ℕ :=
  ∑ i ∈ range S.length, if (candidatePair S i j && decide (m ≤ i) &&
    manhattan1 (siteAt W i) (siteAt W j)) = true then 1 else 0

/-- Unplaced–unplaced contacts whose larger index is below `J`. -/
def ffPart (S : List Bool) (m : ℕ) (W : List Site) (J : ℕ) : ℕ :=
  ∑ j ∈ range J, ffj S m W j

theorem ffPart_succ (S : List Bool) (m : ℕ) (W : List Site) (J : ℕ) :
    ffPart S m W (J + 1) = ffPart S m W J + ffj S m W J := by
  simp [ffPart, Finset.sum_range_succ]

theorem ffPart_base (S : List Bool) (m : ℕ) (W : List Site) : ffPart S m W m = 0 := by
  refine Finset.sum_eq_zero (fun j hj => Finset.sum_eq_zero (fun i _ => ?_))
  rw [Finset.mem_range] at hj
  by_cases hc : candidatePair S i j = true
  · have := candidatePair_lt hc
    have : ¬ m ≤ i := by omega
    simp [this]
  · simp [hc]

theorem ffContacts_eq_ffPart (S : List Bool) (m : ℕ) (W : List Site) (hW : W.length = S.length) :
    ffContacts S m W = ffPart S m W S.length := by
  rw [ffContacts, pairSum_comm, ffPart]
  refine Finset.sum_congr rfl (fun j hj => Finset.sum_congr rfl (fun i _ => ?_))
  rw [Finset.mem_range] at hj
  have : j < W.length := by omega
  simp [this]

theorem ffj_le_slots (S : List Bool) (m : ℕ) (l : List (Fin 4)) (hlen : S.length = l.length + 2)
    (hn : (walkFrom l).Nodup) (j : ℕ) (hj : j < S.length) :
    ffj S m (walkFrom l) j ≤ slots S j := by
  refine le_trans (Finset.sum_le_sum (fun i _ => ?_)) (contactsAt_le_slots S l hlen hn j hj)
  by_cases h : (candidatePair S i j && decide (m ≤ i) &&
      manhattan1 (siteAt (walkFrom l) i) (siteAt (walkFrom l) j)) = true
  · simp only [Bool.and_eq_true] at h
    simp [h]
  · simp [h]

theorem lookupH_spec {h : List (Site × ℕ)} {i : ℕ} {x : Site × ℕ} (hx : lookupH h i = some x) :
    x ∈ h ∧ x.2 = i := by
  unfold lookupH at hx
  exact ⟨List.mem_of_find?_eq_some hx, by simpa using List.find?_some hx⟩

/-- Genuine history: every retained entry records the actual site of its
index in the walk `W`. -/
def GenuineHist (W : List Site) (h : List (Site × ℕ)) : Prop :=
  ∀ x ∈ h, x.2 < W.length ∧ siteAt W x.2 = x.1

theorem loAt_le_ffj (S : List Bool) (m : ℕ) (W : List Site) (h : List (Site × ℕ)) (j : ℕ)
    (hg : GenuineHist W h) :
    loAt S m h j (siteAt W j) ≤ ffj S m W j := by
  refine Finset.sum_le_sum (fun i _ => ?_)
  cases hl : lookupH h i with
  | none => simp
  | some x =>
      obtain ⟨hx, hxi⟩ := lookupH_spec hl
      have hs := (hg x hx).2
      rw [hxi] at hs
      simp only [Option.map_some, Option.getD_some, ← hs]
      split_ifs <;> simp_all

theorem ffj_le_lo_add_out (S : List Bool) (m : ℕ) (W : List Site) (h : List (Site × ℕ)) (j : ℕ)
    (hg : GenuineHist W h) :
    ffj S m W j ≤ loAt S m h j (siteAt W j) + outAt S m h j := by
  unfold ffj loAt outAt
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_le_sum (fun i _ => ?_)
  cases hl : lookupH h i with
  | none => simp only [Option.map_none, Option.getD_none, Bool.and_false, Option.isNone_none,
      Bool.and_true]; split_ifs <;> simp_all
  | some x =>
      obtain ⟨hx, hxi⟩ := lookupH_spec hl
      have hs := (hg x hx).2
      rw [hxi] at hs
      simp only [Option.map_some, Option.getD_some, ← hs, Option.isNone_some, Bool.and_false]
      split_ifs <;> simp_all

/-! ### Invariants of genuine states -/

/-- Invariants of the state reached by a genuine partial conformation `p'`. -/
def Inv (steps : ℕ) (S : List Bool) (p p' : List (Fin 4)) (st : MState) : Prop :=
  st.cur = endFrom first p' ∧ st.dir = lastDir p' ∧
  st.cross = crossContacts S (p.length + 2) (walkFrom p') ∧
  st.idx = p'.length + 2 ∧ GenuineHist (walkFrom p') st.hist ∧
  ∀ l ∈ ext steps p', st.ffLo ≤ ffPart S (p.length + 2) (walkFrom l) (p'.length + 2) ∧
    ffPart S (p.length + 2) (walkFrom l) (p'.length + 2) ≤ st.ffHi

theorem siteAt_walkFrom_last (p' : List (Fin 4)) :
    siteAt (walkFrom p') (p'.length + 1) = endFrom first p' := by
  rcases List.eq_nil_or_concat p' with rfl | ⟨q, e, rfl⟩
  · rfl
  · rw [List.concat_eq_append, walkFrom_append_single, endFrom_append]
    have hl : (walkFrom q).length = q.length + 2 := walkFrom_length q
    simp only [siteAt, List.length_append, List.length_singleton]
    rw [List.getElem?_append_right (by omega)]
    simp [hl, endFrom]

theorem gstate_snoc (S : List Bool) (k : ℕ) {p p' : List (Fin 4)} (h : p <+: p') (d : Fin 4) :
    gstate S k p (p' ++ [d]) =
      memStep S (walkFrom p) (p.length + 2) k (gstate S k p p') d := by
  unfold gstate
  rw [List.drop_append_of_le_length h.length_le, List.foldl_append]
  rfl

theorem gstate_self (S : List Bool) (k : ℕ) (p : List (Fin 4)) :
    gstate S k p p = memInit p := by
  simp [gstate]

theorem memState_inv (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (k : ℕ)
    (p : List (Fin 4)) :
    ∀ s : List (Fin 4), (walkFrom (p ++ s)).Nodup → (p ++ s).length ≤ steps →
      Inv steps S p (p ++ s) (gstate S k p (p ++ s)) := by
  intro s
  induction s using List.reverseRecOn with
  | nil =>
      intro _ _
      rw [List.append_nil, gstate_self]
      refine ⟨rfl, rfl, (crossContacts_self S p).symm, rfl, fun x hx => by simp [memInit] at hx, ?_⟩
      intro l _
      simp only [memInit]
      rw [ffPart_base]
      exact ⟨le_refl _, le_refl _⟩
  | append_singleton s d ih =>
      intro hn hle
      rw [← List.append_assoc] at hn hle ⊢
      set p' := p ++ s with hp'
      have hpp' : p <+: p' := List.prefix_append p s
      have hn' : (walkFrom p').Nodup := hn.sublist (walkFrom_prefix (List.prefix_append _ _)).sublist
      have hle' : p'.length ≤ steps := by simp at hle ⊢; omega
      obtain ⟨hcur, hdir, hcross, hidx, hhist, hff⟩ := ih hn' hle'
      rw [gstate_snoc S k hpp' d]
      set st := gstate S k p p'
      have hWl : (walkFrom p').length = p'.length + 2 := walkFrom_length p'
      have hpre : walkFrom p <+: walkFrom p' := walkFrom_prefix hpp'
      have hW2 : walkFrom (p' ++ [d]) = walkFrom p' ++ [advance (endFrom first p') d] :=
        walkFrom_append_single p' d
      have hpre2 : walkFrom p' <+: walkFrom (p' ++ [d]) := by
        rw [hW2]; exact List.prefix_append _ _
      have hlast : siteAt (walkFrom (p' ++ [d])) (p'.length + 2) = advance (endFrom first p') d := by
        rw [hW2]; simp [siteAt, hWl]
      refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
      · simp only [memStep, hcur, endFrom_append, endFrom]
      · simp only [memStep, lastDir_append_single]
      · simp only [memStep, hcross, hidx, hcur]
        rw [hW2, crossContacts_snoc S p _ _ hpre (by rw [hWl]; simp at hle; omega), hWl]
      · simp only [memStep, hidx, List.length_append, List.length_singleton]
      · intro x hx
        simp only [memStep] at hx
        have hx' := List.mem_of_mem_take hx
        simp only [List.length_append, List.length_singleton, hW2]
        rcases List.mem_cons.1 hx' with rfl | hx'
        · refine ⟨by simp [hWl, hidx], ?_⟩
          simp only [hidx, hcur]
          have : p'.length + 2 - 1 = p'.length + 1 := by omega
          rw [this, ← siteAt_walkFrom_last p']
          exact (siteAt_prefix (List.prefix_append _ _) (by omega)).symm
        · obtain ⟨h1, h2⟩ := hhist x hx'
          refine ⟨by simp; omega, ?_⟩
          rw [← siteAt_prefix (List.prefix_append _ _) h1, h2]
      · intro l hl
        have hl' : l ∈ ext steps p' := by
          rw [mem_ext] at hl ⊢
          exact ⟨hl.1, hl.2.1, (List.prefix_append p' [d]).trans hl.2.2⟩
        obtain ⟨hlo, hhi⟩ := hff l hl'
        rw [mem_ext] at hl
        have hWlpre : walkFrom (p' ++ [d]) <+: walkFrom l := walkFrom_prefix hl.2.2
        have hlenl : (walkFrom l).length = steps + 2 := by rw [walkFrom_length, hl.1]
        have hJ : p'.length + 2 < (walkFrom l).length := by simp at hle; omega
        have hs : siteAt (walkFrom l) (p'.length + 2) = advance st.cur d := by
          rw [← siteAt_prefix hWlpre (by rw [hW2]; simp [hWl]), hlast, hcur]
        have hg : GenuineHist (walkFrom l) st.hist := by
          intro x hx
          obtain ⟨h1, h2⟩ := hhist x hx
          refine ⟨by omega, ?_⟩
          rw [← siteAt_prefix (hpre2.trans hWlpre) h1, h2]
        have e1 := loAt_le_ffj S (p.length + 2) (walkFrom l) st.hist (p'.length + 2) hg
        have e2 := ffj_le_lo_add_out S (p.length + 2) (walkFrom l) st.hist (p'.length + 2) hg
        have e3 := ffj_le_slots S (p.length + 2) l (by omega) hl.2.1 (p'.length + 2) (by omega)
        rw [hs] at e1 e2
        simp only [memStep, List.length_append, List.length_singleton, hidx]
        rw [show p'.length + 1 + 2 = p'.length + 2 + 1 by ring, ffPart_succ]
        refine ⟨by omega, ?_⟩
        unfold upAt
        have := le_min e3 e2
        omega

/-! ### Coverage -/

theorem sum_map_filter_eq' {α : Type*} (l : List α) (P : α → Bool) (f : α → ℚ) :
    ((l.filter P).map f).sum = (l.map (fun x => if P x then f x else 0)).sum :=
  sum_map_filter_eq l P f

/-- **Coverage / weighted domination.** For every genuine extension `p'`
of `p`, the sum over genuine completions of any nonnegative function of
the final memory state is bounded by the relaxed path-sum. -/
theorem memCover_spec (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (k : ℕ)
    (p : List (Fin 4)) (G : MState → ℚ) (hG : ∀ st, 0 ≤ G st) :
    ∀ r (p' : List (Fin 4)), p <+: p' → (walkFrom p').Nodup → p'.length ≤ steps →
      steps - p'.length = r →
      ∑ l ∈ ext steps p', G (gstate S k p l) ≤
        spec (memNext S (walkFrom p) (p.length + 2) k) G r (gstate S k p p') := by
  intro r
  induction r with
  | zero =>
      intro p' _ hn hle hr
      rw [ext_eq_singleton (by omega) hn, Finset.sum_singleton]
      rfl
  | succ r ih =>
      intro p' hpp' hn hle hr
      have hlt : p'.length < steps := by omega
      rw [sum_ext_split hlt, spec, memNext, List.map_map, sum_map_filter_eq, sum_allDirs]
      refine Finset.sum_le_sum (fun d _ => ?_)
      by_cases hnd : (walkFrom (p' ++ [d])).Nodup
      · obtain ⟨s, rfl⟩ := hpp'
        obtain ⟨hcur, hdir, _, _, hhist, _⟩ :=
          memState_inv steps S hS k p s hn (by omega)
        have hrev : d ≠ opp (lastDir (p ++ s)) := by
          rintro rfl; exact not_nodup_reverse _ hnd
        have hnew : advance (endFrom first (p ++ s)) d ∉ walkFrom (p ++ s) := by
          intro hm
          rw [walkFrom_append_single, List.nodup_append] at hnd
          exact hnd.2.2 _ hm _ (List.mem_singleton_self _) rfl
        have hpre : walkFrom p <+: walkFrom (p ++ s) := walkFrom_prefix (List.prefix_append _ _)
        have hall : memAllowed (walkFrom p) (gstate S k p (p ++ s)) d = true := by
          have hc : (walkFrom p).contains (advance (endFrom first (p ++ s)) d) = false := by
            simpa using fun hm => hnew (hpre.subset hm)
          have hh : (gstate S k p (p ++ s)).hist.any
              (fun x => x.1 == advance (endFrom first (p ++ s)) d) = false := by
            rw [List.any_eq_false]
            intro x hx hxe
            have hxe' : x.1 = advance (endFrom first (p ++ s)) d := by simpa using hxe
            obtain ⟨h1, h2⟩ := hhist x hx
            apply hnew
            rw [← hxe', ← h2, siteAt_eq_getElem h1]
            exact List.getElem_mem h1
          simp only [memAllowed, hcur, hdir, Bool.and_eq_true, decide_eq_true_eq,
            Bool.not_eq_true', hc, hh]
          exact ⟨⟨hrev, trivial⟩, trivial⟩
        rw [if_pos hall]
        have := ih ((p ++ s) ++ [d]) ((List.prefix_append p s).trans (List.prefix_append _ _)) hnd
          (by simp at hle hlt ⊢; omega) (by simp at hr ⊢; omega)
        rw [gstate_snoc S k (List.prefix_append p s)] at this
        exact this
      · rw [ext_eq_empty_of_not_nodup hnd, Finset.sum_empty]
        split_ifs
        · exact spec_nonneg _ _ hG _ _
        · exact le_refl _

/-- **Coverage theorem used by the checker.** -/
theorem memCover (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (k : ℕ)
    (p : List (Fin 4)) (hn : (walkFrom p).Nodup) (hle : p.length ≤ steps)
    (G : MState → ℚ) (hG : ∀ st, 0 ≤ G st) :
    ∑ l ∈ ext steps p, G (gstate S k p l) ≤ wsum G (memTable steps S k p) := by
  rw [memTable, wsum_dpM, wsum_cons, wsum_nil, one_mul, add_zero]
  have := memCover_spec steps S hS k p G hG _ p (List.prefix_refl p) hn hle rfl
  rwa [gstate_self] at this

/-! ### Contact bounds from final states -/

/-- Sites known to lie on a completion with final state `st`. -/
def knownSites (p : List (Fin 4)) (st : MState) : List Site :=
  walkFrom p ++ st.cur :: st.hist.map Prod.fst

/-- **Contact bounds for genuine completions.** -/
theorem memContacts_bounds (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (k : ℕ)
    {p l : List (Fin 4)} (hl : l ∈ ext steps p) :
    knownContacts S p + (gstate S k p l).cross + (gstate S k p l).ffLo ≤
      parityContactCount S (walkFrom l) ∧
    parityContactCount S (walkFrom l) ≤
      knownContacts S p + (gstate S k p l).cross + (gstate S k p l).ffHi ∧
    parityContactCount S (walkFrom l) ≤ geomCap (steps + 2) (knownSites p (gstate S k p l)) := by
  have hl' := hl
  rw [mem_ext] at hl'
  obtain ⟨hlen, hn, ⟨s, rfl⟩⟩ := hl'
  obtain ⟨hcur, _, hcross, _, hhist, hff⟩ :=
    memState_inv steps S hS k p s hn (by omega)
  have hself : p ++ s ∈ ext steps (p ++ s) := mem_ext.2 ⟨hlen, hn, List.prefix_refl _⟩
  obtain ⟨hlo, hhi⟩ := hff _ hself
  have hWlen : (walkFrom (p ++ s)).length = S.length := by rw [walkFrom_length]; omega
  have hfe := ffContacts_eq_ffPart S (p.length + 2) (walkFrom (p ++ s)) hWlen
  rw [show (p ++ s).length + 2 = S.length by omega] at hlo hhi
  have hdec := contacts_decompose S (List.prefix_append p s)
  refine ⟨by omega, by omega, ?_⟩
  apply parityContacts_le_geomCap S hS (FrontierCheck.mem_ext_valid hl)
  intro a ha
  simp only [knownSites, List.mem_append, List.mem_cons, List.mem_map] at ha
  rcases ha with ha | rfl | ⟨x, hx, rfl⟩
  · exact (walkFrom_prefix (List.prefix_append p s)).subset ha
  · rw [hcur]; exact endFrom_mem_walkFrom _
  · obtain ⟨h1, h2⟩ := hhist x hx
    rw [← h2, siteAt_eq_getElem h1]; exact List.getElem_mem h1

end MemTransfer
end FoldBounds
