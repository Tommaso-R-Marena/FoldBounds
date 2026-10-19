module

public import FoldBounds.LatticeContacts

/-!
# The conformation tree: completion sets, partition, and subtree bounds

A *prefix* is a word `p : List (Fin 4)` of directions after the fixed first
bond.  Its completion set `ext steps p` is the set of complete valid
conformations (encoded as direction words) having `p` as a prefix.  This
file proves the structural facts needed by the non-enumerative frontier
checker:

* `ext steps [] = validWords steps` (the root covers everything);
* for a nonterminal prefix the four one-step extensions *partition* the
  completion set (`ext_eq_biUnion`, `ext_pairwiseDisjoint`, `sum_ext_split`);
* invalid prefixes have empty completion sets;
* the counting bound `#ext steps p ≤ 3 ^ (steps - |p|)`
  (`card_ext_le_pow_three`), via the impossibility of immediate reversal;
* the exact contact decomposition `k = known + cross + futureFuture`
  for every completion, and the resulting contact bounds.

All sums here are *semantic*: the checker in `FrontierCheck` never
evaluates `validWords`.
-/
@[expose] public section

namespace FoldBounds
namespace Frontier

open LatticeModel Finset

/-! ### Valid complete words -/

/-- Valid complete direction words (semantic object; never enumerated by
the checker). -/
def validWords (steps : ℕ) : Finset (List (Fin 4)) :=
  Finset.univ.image (fun c : Conformation steps => List.ofFn c.1)

theorem mem_validWords {steps : ℕ} {l : List (Fin 4)} :
    l ∈ validWords steps ↔ l.length = steps ∧ (walkFrom l).Nodup := by
  constructor
  · intro h
    obtain ⟨c, _, rfl⟩ := Finset.mem_image.1 h
    exact ⟨by simp, c.2⟩
  · rintro ⟨hl, hn⟩
    subst hl
    refine Finset.mem_image.2 ⟨⟨fun i => l[i.1], ?_⟩, Finset.mem_univ _, ?_⟩
    · rw [List.ofFn_getElem]; exact hn
    · exact List.ofFn_getElem l

theorem ofFn_val_injective (steps : ℕ) :
    Function.Injective (fun c : Conformation steps => List.ofFn c.1) := by
  intro c c' h
  exact Subtype.ext (List.ofFn_injective h)

/-- Sums over conformations are sums over valid words. -/
theorem sum_conformation_eq_sum_validWords {M : Type*} [AddCommMonoid M]
    (steps : ℕ) (F : List (Fin 4) → M) (T : Finset (Conformation steps)) :
    ∑ c ∈ T, F (List.ofFn c.1) = ∑ l ∈ T.image (fun c => List.ofFn c.1), F l := by
  rw [Finset.sum_image (fun c _ c' _ h => ofFn_val_injective steps h)]

/-! ### Completion sets -/

/-- Completion set of a prefix: complete valid words extending it. -/
def ext (steps : ℕ) (p : List (Fin 4)) : Finset (List (Fin 4)) :=
  (validWords steps).filter (fun l => p <+: l)

theorem mem_ext {steps : ℕ} {p l : List (Fin 4)} :
    l ∈ ext steps p ↔ l.length = steps ∧ (walkFrom l).Nodup ∧ p <+: l := by
  simp [ext, mem_validWords, and_assoc]

theorem ext_nil (steps : ℕ) : ext steps [] = validWords steps := by
  ext l; simp [ext]

theorem ext_eq_empty_of_not_nodup {steps : ℕ} {p : List (Fin 4)}
    (h : ¬ (walkFrom p).Nodup) : ext steps p = ∅ := by
  ext l
  simp only [mem_ext, Finset.notMem_empty, iff_false, not_and]
  intro _ hn hp
  exact h (hn.sublist (walkFrom_prefix hp).sublist)

theorem ext_eq_empty_of_length_gt {steps : ℕ} {p : List (Fin 4)}
    (h : steps < p.length) : ext steps p = ∅ := by
  ext l
  simp only [mem_ext, Finset.notMem_empty, iff_false, not_and]
  intro hl _ hp
  have := hp.length_le
  omega

theorem ext_subset_singleton_of_length_eq {steps : ℕ} {p : List (Fin 4)}
    (h : p.length = steps) : ext steps p ⊆ {p} := by
  intro l hl
  rw [mem_ext] at hl
  rw [Finset.mem_singleton]
  exact (hl.2.2.eq_of_length (by omega)).symm

theorem ext_eq_singleton {steps : ℕ} {p : List (Fin 4)}
    (h : p.length = steps) (hn : (walkFrom p).Nodup) : ext steps p = {p} := by
  apply Finset.Subset.antisymm (ext_subset_singleton_of_length_eq h)
  rw [Finset.singleton_subset_iff, mem_ext]
  exact ⟨h, hn, List.prefix_refl p⟩

theorem prefix_getElem? {l₁ l₂ : List (Fin 4)} (h : l₁ <+: l₂) {n : ℕ}
    (hn : n < l₁.length) : l₁[n]? = l₂[n]? := by
  obtain ⟨t, rfl⟩ := h
  rw [List.getElem?_append_left hn]

theorem prefix_concat_of_prefix {p l : List (Fin 4)} (h : p <+: l)
    (hlen : p.length < l.length) : p ++ [l[p.length]] <+: l := by
  obtain ⟨t, rfl⟩ := h
  cases t with
  | nil => simp at hlen
  | cons a t =>
      simp only [List.getElem_append_right (le_refl _), Nat.sub_self,
        List.getElem_cons_zero]
      exact ⟨t, by simp⟩

/-- One-step extensions cover the completion set of a nonterminal prefix. -/
theorem ext_eq_biUnion {steps : ℕ} {p : List (Fin 4)} (h : p.length < steps) :
    ext steps p = (Finset.univ : Finset (Fin 4)).biUnion (fun d => ext steps (p ++ [d])) := by
  ext l
  simp only [Finset.mem_biUnion, Finset.mem_univ, true_and, mem_ext]
  constructor
  · rintro ⟨hl, hn, hp⟩
    exact ⟨l[p.length]'(by omega), hl, hn, prefix_concat_of_prefix hp (by omega)⟩
  · rintro ⟨d, hl, hn, hp⟩
    exact ⟨hl, hn, (List.prefix_append p [d]).trans hp⟩

/-- Distinct one-step extensions have disjoint completion sets. -/
theorem ext_pairwiseDisjoint (steps : ℕ) (p : List (Fin 4)) :
    (↑(Finset.univ : Finset (Fin 4)) : Set (Fin 4)).PairwiseDisjoint
      (fun d => ext steps (p ++ [d])) := by
  intro d _ e _ hde
  rw [Function.onFun, Finset.disjoint_left]
  intro l hd he
  rw [mem_ext] at hd he
  apply hde
  have h1 := prefix_getElem? hd.2.2 (n := p.length) (by simp)
  have h2 := prefix_getElem? he.2.2 (n := p.length) (by simp)
  simp only [List.getElem?_append_right (le_refl _), Nat.sub_self,
    List.getElem?_cons_zero] at h1 h2
  rw [← h1] at h2
  exact (Option.some_injective _ h2).symm

/-- **Partition of a completion sum** over the four one-step children. -/
theorem sum_ext_split {M : Type*} [AddCommMonoid M] {steps : ℕ} {p : List (Fin 4)}
    (h : p.length < steps) (F : List (Fin 4) → M) :
    ∑ l ∈ ext steps p, F l = ∑ d : Fin 4, ∑ l ∈ ext steps (p ++ [d]), F l := by
  rw [ext_eq_biUnion h, Finset.sum_biUnion (ext_pairwiseDisjoint steps p)]

/-! ### Immediate reversal is impossible -/

/-- The opposite cardinal direction. -/
def opp (d : Fin 4) : Fin 4 :=
  if d.val = 0 then 1 else if d.val = 1 then 0 else if d.val = 2 then 3 else 2

theorem advance_advance_opp (s : Site) (d : Fin 4) :
    advance (advance s d) (opp d) = s := by
  obtain ⟨a, b⟩ := s
  fin_cases d <;> simp [advance, opp]

/-- Direction of the last bond of the walk decoded from `p` (the fixed
first bond has direction `0`). -/
def lastDir (p : List (Fin 4)) : Fin 4 := p.getLastD 0

theorem lastDir_append_single (p : List (Fin 4)) (d : Fin 4) :
    lastDir (p ++ [d]) = d := by
  simp [lastDir]

theorem not_nodup_reverse (p : List (Fin 4)) :
    ¬ (walkFrom (p ++ [opp (lastDir p)])).Nodup := by
  rcases List.eq_nil_or_concat p with rfl | ⟨q, e, rfl⟩
  · decide
  · rw [List.concat_eq_append, lastDir_append_single, List.append_assoc,
      walkFrom_append]
    simp only [List.singleton_append, trace, advance_advance_opp]
    intro hn
    rw [List.nodup_append] at hn
    exact hn.2.2 _ (endFrom_mem_walkFrom q) _ (by simp) rfl

/-! ### Counting bound -/

theorem card_ext_le_one_of_length_ge {steps : ℕ} {p : List (Fin 4)}
    (h : steps ≤ p.length) : (ext steps p).card ≤ 1 := by
  rcases Nat.lt_or_eq_of_le h with h | h
  · rw [ext_eq_empty_of_length_gt h]; simp
  · exact (Finset.card_le_card (ext_subset_singleton_of_length_eq h.symm)).trans (by simp)

/-- **Non-enumerative counting bound.** A prefix with `|p|` free
directions has at most `3 ^ (steps - |p|)` valid completions. -/
theorem card_ext_le_pow_three (steps : ℕ) (p : List (Fin 4)) :
    (ext steps p).card ≤ 3 ^ (steps - p.length) := by
  induction h : steps - p.length generalizing p with
  | zero => exact card_ext_le_one_of_length_ge (by omega) |>.trans (by simp)
  | succ r ih =>
      have hlt : p.length < steps := by omega
      rw [ext_eq_biUnion hlt, Finset.card_biUnion (ext_pairwiseDisjoint steps p)]
      have hrev : ext steps (p ++ [opp (lastDir p)]) = ∅ :=
        ext_eq_empty_of_not_nodup (not_nodup_reverse p)
      rw [← Finset.add_sum_erase _ _ (Finset.mem_univ (opp (lastDir p))), hrev,
        Finset.card_empty, zero_add]
      have hchild : ∀ d ∈ (Finset.univ : Finset (Fin 4)).erase (opp (lastDir p)),
          (ext steps (p ++ [d])).card ≤ 3 ^ r := by
        intro d _
        exact ih (p ++ [d]) (by simp; omega)
      calc ∑ d ∈ (Finset.univ : Finset (Fin 4)).erase (opp (lastDir p)),
              (ext steps (p ++ [d])).card
          ≤ ((Finset.univ : Finset (Fin 4)).erase (opp (lastDir p))).card • 3 ^ r :=
            Finset.sum_le_card_nsmul _ _ _ hchild
        _ = 3 ^ (r + 1) := by
            rw [Finset.card_erase_of_mem (Finset.mem_univ _)]
            simp [pow_succ, mul_comm]

/-! ### Pair sums and the contact decomposition -/

/-- Number of index pairs `(i, j)` (both below `|S|`) satisfying `P`. -/
def pairSum (S : List Bool) (P : ℕ → ℕ → Bool) : ℕ :=
  ∑ i ∈ range S.length, ∑ j ∈ range S.length, if P i j then 1 else 0

theorem pairSum_mono (S : List Bool) {P Q : ℕ → ℕ → Bool}
    (h : ∀ i j, i < S.length → j < S.length → P i j = true → Q i j = true) :
    pairSum S P ≤ pairSum S Q := by
  unfold pairSum
  refine Finset.sum_le_sum (fun i hi => Finset.sum_le_sum (fun j hj => ?_))
  simp only [Finset.mem_range] at hi hj
  by_cases hp : P i j = true
  · simp [hp, h i j hi hj hp]
  · simp [hp]

theorem pairSum_congr (S : List Bool) {P Q : ℕ → ℕ → Bool}
    (h : ∀ i j, i < S.length → j < S.length → P i j = Q i j) :
    pairSum S P = pairSum S Q := by
  unfold pairSum
  refine Finset.sum_congr rfl (fun i hi => Finset.sum_congr rfl (fun j hj => ?_))
  simp only [Finset.mem_range] at hi hj
  rw [h i j hi hj]

theorem pairSum_split (S : List Bool) {P Q R : ℕ → ℕ → Bool}
    (h : ∀ i j, i < S.length → j < S.length →
      (P i j = true ↔ (Q i j = true ∨ R i j = true)) ∧ ¬ (Q i j = true ∧ R i j = true)) :
    pairSum S P = pairSum S Q + pairSum S R := by
  unfold pairSum
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun j hj => ?_)
  simp only [Finset.mem_range] at hi hj
  obtain ⟨h1, h2⟩ := h i j hi hj
  by_cases hq : Q i j = true <;> by_cases hr : R i j = true <;> simp_all

/-- Site lookup with the same default convention as `contactCount`. -/
def siteAt (W : List Site) (i : ℕ) : Site := (W[i]?).getD origin

theorem parityContactCount_eq_pairSum (S : List Bool) (W : List Site) :
    parityContactCount S W = pairSum S (fun i j =>
      candidatePair S i j && decide (j < W.length) && manhattan1 (siteAt W i) (siteAt W j)) :=
  rfl

theorem candidatePair_lt {S : List Bool} {i j : ℕ} (h : candidatePair S i j = true) :
    i + 3 ≤ j := by
  simp only [candidatePair, Bool.and_eq_true, decide_eq_true_eq] at h
  exact h.1.1.1

/-- Contacts already fixed by the prefix (Python: `_contacts(path, pairs)`). -/
def knownContacts (S : List Bool) (p : List (Fin 4)) : ℕ :=
  parityContactCount S (walkFrom p)

/-- Candidate pairs whose larger index is not yet placed (`m` = number of
placed sites). Python: `remaining_potential` in `_score_range`. -/
def futurePairs (S : List Bool) (m : ℕ) : ℕ :=
  pairSum S (fun i j => candidatePair S i j && decide (m ≤ j))

/-- Candidate pairs with one end placed and one end unplaced. -/
def crossPairs (S : List Bool) (m : ℕ) : ℕ :=
  pairSum S (fun i j => candidatePair S i j && decide (i < m) && decide (m ≤ j))

/-- Candidate pairs with both ends unplaced. -/
def ffPairs (S : List Bool) (m : ℕ) : ℕ :=
  pairSum S (fun i j => candidatePair S i j && decide (m ≤ i))

/-- Realised placed–unplaced contacts of a walk `W` relative to `m`
placed prefix sites. -/
def crossContacts (S : List Bool) (m : ℕ) (W : List Site) : ℕ :=
  pairSum S (fun i j => candidatePair S i j && decide (i < m) && decide (m ≤ j) &&
    decide (j < W.length) && manhattan1 (siteAt W i) (siteAt W j))

/-- Realised unplaced–unplaced contacts of a walk `W`. -/
def ffContacts (S : List Bool) (m : ℕ) (W : List Site) : ℕ :=
  pairSum S (fun i j => candidatePair S i j && decide (m ≤ i) &&
    decide (j < W.length) && manhattan1 (siteAt W i) (siteAt W j))

theorem futurePairs_eq (S : List Bool) (m : ℕ) :
    futurePairs S m = crossPairs S m + ffPairs S m := by
  apply pairSum_split
  intro i j _ _
  have := @candidatePair_lt S i j
  by_cases hc : candidatePair S i j = true
  · have := this hc
    simp only [hc, Bool.true_and, Bool.and_eq_true, decide_eq_true_eq]
    omega
  · simp [hc]

theorem crossContacts_le (S : List Bool) (m : ℕ) (W : List Site) :
    crossContacts S m W ≤ crossPairs S m := by
  apply pairSum_mono
  intro i j _ _ h
  simp only [Bool.and_eq_true] at h ⊢
  tauto

theorem ffContacts_le (S : List Bool) (m : ℕ) (W : List Site) :
    ffContacts S m W ≤ ffPairs S m := by
  apply pairSum_mono
  intro i j _ _ h
  simp only [Bool.and_eq_true] at h ⊢
  tauto

theorem siteAt_prefix {W V : List Site} (h : W <+: V) {i : ℕ} (hi : i < W.length) :
    siteAt W i = siteAt V i := by
  obtain ⟨t, rfl⟩ := h
  simp [siteAt, List.getElem?_append_left hi]

/-- **Exact contact decomposition.** For every word `l` extending `p`,
the Python-style contact count of `l` splits into the contacts fixed by
the prefix, placed–unplaced contacts, and unplaced–unplaced contacts. -/
theorem contacts_decompose (S : List Bool) {p l : List (Fin 4)} (h : p <+: l) :
    parityContactCount S (walkFrom l) =
      knownContacts S p + crossContacts S (p.length + 2) (walkFrom l) +
        ffContacts S (p.length + 2) (walkFrom l) := by
  have hpre := walkFrom_prefix h
  have hplen : (walkFrom p).length = p.length + 2 := walkFrom_length p
  have hle : p.length + 2 ≤ (walkFrom l).length := by
    have := hpre.length_le; omega
  rw [knownContacts, crossContacts, ffContacts, parityContactCount_eq_pairSum,
    parityContactCount_eq_pairSum]
  rw [pairSum_split S (Q := fun i j => candidatePair S i j && decide (j < (walkFrom p).length) &&
      manhattan1 (siteAt (walkFrom p) i) (siteAt (walkFrom p) j))
    (R := fun i j => (candidatePair S i j && decide (i < p.length + 2) &&
      decide (p.length + 2 ≤ j) && decide (j < (walkFrom l).length) &&
      manhattan1 (siteAt (walkFrom l) i) (siteAt (walkFrom l) j)) ||
      (candidatePair S i j && decide (p.length + 2 ≤ i) &&
      decide (j < (walkFrom l).length) &&
      manhattan1 (siteAt (walkFrom l) i) (siteAt (walkFrom l) j)))]
  · rw [add_assoc]
    congr 1
    apply pairSum_split
    intro i j _ _
    simp only [Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq]
    refine ⟨by simp, ?_⟩
    rintro ⟨⟨⟨⟨⟨_, h1⟩, _⟩, _⟩, _⟩, ⟨⟨_, h2⟩, _⟩, _⟩
    omega
  · intro i j _ _
    by_cases hc : candidatePair S i j = true
    · have hij := candidatePair_lt hc
      by_cases hj : j < p.length + 2
      · rw [siteAt_prefix hpre (show i < (walkFrom p).length by omega),
          siteAt_prefix hpre (show j < (walkFrom p).length by omega)]
        simp only [hc, Bool.true_and, Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq]
        constructor
        · constructor
          · rintro ⟨_, hm⟩; exact Or.inl ⟨by omega, hm⟩
          · rintro (⟨_, hm⟩ | h' | h')
            · exact ⟨by omega, hm⟩
            · exact absurd h'.1.1.2 (by omega)
            · exact absurd h'.1.1 (by omega)
        · rintro ⟨⟨_, _⟩, _ | _⟩ <;> omega
      · simp only [hc, Bool.true_and, Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq]
        constructor
        · constructor
          · rintro ⟨hj', hm⟩
            by_cases hi : i < p.length + 2
            · exact Or.inr (Or.inl ⟨⟨⟨hi, by omega⟩, hj'⟩, hm⟩)
            · exact Or.inr (Or.inr ⟨⟨by omega, hj'⟩, hm⟩)
          · rintro (⟨h1, _⟩ | ⟨⟨_, hj'⟩, hm⟩ | ⟨⟨_, hj'⟩, hm⟩)
            · omega
            · exact ⟨hj', hm⟩
            · exact ⟨hj', hm⟩
        · rintro ⟨⟨h1, _⟩, _⟩; omega
    · simp [hc]

/-- Contacts of a completion are at least the contacts fixed by its prefix. -/
theorem knownContacts_le (S : List Bool) {p l : List (Fin 4)} (h : p <+: l) :
    knownContacts S p ≤ parityContactCount S (walkFrom l) := by
  rw [contacts_decompose S h]; omega

/-- **Basic contact bound** `k(c) ≤ K_max(p)`: known contacts plus the
number of candidate pairs whose larger index is still unplaced. -/
theorem contacts_le_kmax (S : List Bool) {p l : List (Fin 4)} (h : p <+: l) :
    parityContactCount S (walkFrom l) ≤ knownContacts S p + futurePairs S (p.length + 2) := by
  rw [contacts_decompose S h, futurePairs_eq]
  have := crossContacts_le S (p.length + 2) (walkFrom l)
  have := ffContacts_le S (p.length + 2) (walkFrom l)
  omega

/-- Refined two-sided contact bounds used by the transfer relaxation. -/
theorem contacts_le_cross_bound (S : List Bool) {p l : List (Fin 4)} (h : p <+: l) :
    knownContacts S p + crossContacts S (p.length + 2) (walkFrom l) ≤
      parityContactCount S (walkFrom l) ∧
    parityContactCount S (walkFrom l) ≤
      knownContacts S p + crossContacts S (p.length + 2) (walkFrom l) +
        ffPairs S (p.length + 2) := by
  rw [contacts_decompose S h]
  have := ffContacts_le S (p.length + 2) (walkFrom l)
  omega

/-! ### Semantic contact count on words -/

/-- Contact count of a complete word, via the Lean reference scorer. -/
def wordContacts (S : List Bool) (l : List (Fin 4)) : ℕ := contactCount S (walkFrom l)

theorem wordContacts_eq_parity {steps : ℕ} (S : List Bool) (hS : S.length = steps + 2)
    {l : List (Fin 4)} (hl : l ∈ validWords steps) :
    wordContacts S l = parityContactCount S (walkFrom l) := by
  rw [mem_validWords] at hl
  exact contactCount_eq_parityContactCount S l (by omega)

end Frontier
end FoldBounds
