module

public import FoldBounds.Frontier

/-!
# Geometric degree caps on future contacts

On the square lattice a site has four neighbours.  For residue `j` of a
self-avoiding walk, the chain neighbours `j - 1` and (if it exists) `j + 1`
occupy two (resp. one) of them, so at most `slots j = 2` (resp. `3` for the
last residue) earlier residues can be in contact with `j`
(`contactsAt_le_slots`).  Assigning every future contact to its larger
index gives the caps

* `futureCap S m = ∑_{j ≥ m} min (slots j) (#candidates (·, j))`,
* `ffCap S m = ∑_{j ≥ m} min (slots j) (#candidates (i, j) with i ≥ m)`,

which replace the quadratic candidate-pair counts `futurePairs`/`ffPairs`
by quantities linear in the number of unplaced residues
(`contacts_le_known_futureCap`, `ffContacts_le_ffCap`,
`futureCap_le_futurePairs`).
-/
@[expose] public section

namespace FoldBounds
namespace Frontier

open LatticeModel Finset

theorem adjacent_symm {a b : Site} (h : adjacent a b = true) : adjacent b a = true := by
  obtain ⟨a1, a2⟩ := a
  obtain ⟨b1, b2⟩ := b
  simp only [adjacent, decide_eq_true_eq] at h ⊢
  omega

theorem adjacentChain_getElem {W : List Site} (h : adjacentChain W = true) {k : ℕ}
    (hk : k + 1 < W.length) : adjacent W[k] W[k + 1] = true := by
  induction W generalizing k with
  | nil => simp at hk
  | cons a W ih =>
      cases W with
      | nil => simp at hk
      | cons b W =>
          simp only [adjacentChain, Bool.and_eq_true] at h
          cases k with
          | zero => simpa using h.1
          | succ k =>
              simp only [List.getElem_cons_succ]
              exact ih h.2 (by simpa using hk)

/-- The four lattice neighbours of a site. -/
def nbrs (s : Site) : Finset Site := (Finset.univ : Finset (Fin 4)).image (advance s)

theorem card_nbrs_le (s : Site) : (nbrs s).card ≤ 4 :=
  (Finset.card_image_le).trans (by simp)

theorem mem_nbrs {s a : Site} (h : adjacent s a = true) : a ∈ nbrs s :=
  Finset.mem_image.2 ⟨stepDir s a, Finset.mem_univ _, advance_stepDir h⟩

/-- Free contact slots of residue `j`: two for interior residues, three for
the last residue. -/
def slots (S : List Bool) (j : ℕ) : ℕ := if j + 1 < S.length then 2 else 3

theorem siteAt_eq_getElem {W : List Site} {i : ℕ} (h : i < W.length) : siteAt W i = W[i] := by
  simp [siteAt, List.getElem?_eq_getElem h]

/-- **Degree lemma.** In a complete self-avoiding walk, residue `j` is in
contact with at most `slots S j` earlier candidate residues. -/
theorem contactsAt_le_slots (S : List Bool) (l : List (Fin 4)) (hlen : S.length = l.length + 2)
    (hn : (walkFrom l).Nodup) (j : ℕ) (hj : j < S.length) :
    ∑ i ∈ range S.length, (if (candidatePair S i j &&
      manhattan1 (siteAt (walkFrom l) i) (siteAt (walkFrom l) j)) = true then 1 else 0) ≤
      slots S j := by
  set W := walkFrom l with hWdef
  have hW : W.length = S.length := by rw [hWdef, walkFrom_length]; omega
  have hadj : adjacentChain W = true := walkFrom_adjacentChain l
  rw [Finset.sum_boole]
  simp only [Nat.cast_id]
  by_cases hj3 : j < 3
  · have : ((range S.length).filter (fun i => (candidatePair S i j &&
        manhattan1 (siteAt W i) (siteAt W j)) = true)) = ∅ := by
      apply Finset.filter_false_of_mem
      intro i _ hc
      simp only [Bool.and_eq_true] at hc
      have := candidatePair_lt hc.1
      omega
    rw [this]; simp
  · have hjW : j < W.length := by omega
    have hj1 : j - 1 < W.length := by omega
    set s := W[j] with hs
    set a := W[j - 1] with ha
    have hadj_a : adjacent s a = true := by
      have := adjacentChain_getElem hadj (k := j - 1) (by omega)
      have e : j - 1 + 1 = j := by omega
      simp only [e] at this
      exact adjacent_symm this
    have ha_mem : a ∈ nbrs s := mem_nbrs hadj_a
    -- target set of admissible partner sites
    let N : Finset Site :=
      if h : j + 1 < W.length then ((nbrs s).erase a).erase W[j + 1] else (nbrs s).erase a
    have hN : N.card ≤ slots S j := by
      unfold slots
      by_cases h : j + 1 < W.length
      · have hb : W[j + 1] ∈ (nbrs s).erase a := by
          refine Finset.mem_erase.2 ⟨?_, mem_nbrs (adjacentChain_getElem hadj h)⟩
          intro he
          have := (List.Nodup.getElem_inj_iff hn).1 he
          omega
        simp only [N, h, dif_pos]
        rw [Finset.card_erase_of_mem hb, Finset.card_erase_of_mem ha_mem]
        have := card_nbrs_le s
        rw [if_pos (by omega)]
        omega
      · simp only [N, h, dif_neg, not_false_eq_true]
        rw [Finset.card_erase_of_mem ha_mem]
        have := card_nbrs_le s
        split_ifs <;> omega
    refine le_trans (Finset.card_le_card_of_injOn (fun i => siteAt W i) ?_ ?_) hN
    · intro i hi
      simp only [Finset.coe_filter, Finset.mem_range, Bool.and_eq_true, Set.mem_setOf_eq] at hi
      obtain ⟨hiS, hc, hm⟩ := hi
      have hij := candidatePair_lt hc
      have hiW : i < W.length := by omega
      rw [siteAt_eq_getElem hiW, siteAt_eq_getElem hjW, ← adjacent_eq_manhattan1] at hm
      have hmem : W[i] ∈ nbrs s := mem_nbrs (adjacent_symm hm)
      have hne_a : W[i] ≠ a := by
        intro he
        have := (List.Nodup.getElem_inj_iff hn).1 he
        omega
      simp only [Finset.mem_coe, N]
      rw [siteAt_eq_getElem hiW]
      split_ifs with h
      · refine Finset.mem_erase.2 ⟨?_, Finset.mem_erase.2 ⟨hne_a, hmem⟩⟩
        intro he
        have := (List.Nodup.getElem_inj_iff hn).1 he
        omega
      · exact Finset.mem_erase.2 ⟨hne_a, hmem⟩
    · intro i hi i' hi' he
      simp only [Finset.coe_filter, Finset.mem_range, Set.mem_setOf_eq] at hi hi'
      have h1 : i < W.length := by omega
      have h2 : i' < W.length := by omega
      simp only [siteAt_eq_getElem h1, siteAt_eq_getElem h2] at he
      exact (List.Nodup.getElem_inj_iff hn).1 he

/-- Degree-capped number of future contacts (`m` placed sites). -/
def futureCap (S : List Bool) (m : ℕ) : ℕ :=
  ∑ j ∈ range S.length, if m ≤ j then
    min (slots S j) (∑ i ∈ range S.length, if candidatePair S i j = true then 1 else 0) else 0

/-- Degree-capped number of unplaced–unplaced contacts. -/
def ffCap (S : List Bool) (m : ℕ) : ℕ :=
  ∑ j ∈ range S.length, if m ≤ j then
    min (slots S j)
      (∑ i ∈ range S.length, if (candidatePair S i j && decide (m ≤ i)) = true then 1 else 0)
  else 0

/-- Realised future contacts (larger index unplaced). -/
def futureContacts (S : List Bool) (m : ℕ) (W : List Site) : ℕ :=
  pairSum S (fun i j => candidatePair S i j && decide (m ≤ j) && decide (j < W.length) &&
    manhattan1 (siteAt W i) (siteAt W j))

theorem cross_add_ff (S : List Bool) (m : ℕ) (W : List Site) :
    crossContacts S m W + ffContacts S m W = futureContacts S m W := by
  symm
  apply pairSum_split
  intro i j _ _
  by_cases hc : candidatePair S i j = true
  · have := candidatePair_lt hc
    simp only [hc, Bool.true_and, Bool.and_eq_true, decide_eq_true_eq]
    constructor
    · constructor
      · rintro ⟨⟨h1, h2⟩, h3⟩
        by_cases hi : i < m
        · exact Or.inl ⟨⟨⟨hi, h1⟩, h2⟩, h3⟩
        · exact Or.inr ⟨⟨by omega, h2⟩, h3⟩
      · rintro (⟨⟨⟨_, h1⟩, h2⟩, h3⟩ | ⟨⟨h1, h2⟩, h3⟩)
        · exact ⟨⟨h1, h2⟩, h3⟩
        · exact ⟨⟨by omega, h2⟩, h3⟩
    · rintro ⟨⟨⟨⟨h1, _⟩, _⟩, _⟩, ⟨⟨h2, _⟩, _⟩⟩; omega
  · simp [hc]

theorem pairSum_comm (S : List Bool) (P : ℕ → ℕ → Bool) :
    pairSum S P = ∑ j ∈ range S.length, ∑ i ∈ range S.length, if P i j then 1 else 0 :=
  Finset.sum_comm

theorem futureContacts_le_futureCap (S : List Bool) (l : List (Fin 4))
    (hlen : S.length = l.length + 2) (hn : (walkFrom l).Nodup) (m : ℕ) :
    futureContacts S m (walkFrom l) ≤ futureCap S m := by
  rw [futureContacts, pairSum_comm, futureCap]
  refine Finset.sum_le_sum (fun j hj => ?_)
  rw [Finset.mem_range] at hj
  split_ifs with hm
  · refine le_min ?_ ?_
    · refine le_trans (Finset.sum_le_sum (fun i _ => ?_)) (contactsAt_le_slots S l hlen hn j hj)
      by_cases h : (candidatePair S i j && decide (m ≤ j) && decide (j < (walkFrom l).length) &&
        manhattan1 (siteAt (walkFrom l) i) (siteAt (walkFrom l) j)) = true
      · simp only [Bool.and_eq_true] at h
        simp [h]
      · simp [h]
    · refine Finset.sum_le_sum (fun i _ => ?_)
      by_cases h : candidatePair S i j = true
      · simp only [h, if_true]; split_ifs <;> omega
      · simp [h]
  · refine le_of_eq (Finset.sum_eq_zero (fun i _ => ?_))
    simp [hm]

theorem ffContacts_le_ffCap (S : List Bool) (l : List (Fin 4))
    (hlen : S.length = l.length + 2) (hn : (walkFrom l).Nodup) (m : ℕ) :
    ffContacts S m (walkFrom l) ≤ ffCap S m := by
  rw [ffContacts, pairSum_comm, ffCap]
  refine Finset.sum_le_sum (fun j hj => ?_)
  rw [Finset.mem_range] at hj
  split_ifs with hm
  · refine le_min ?_ ?_
    · refine le_trans (Finset.sum_le_sum (fun i _ => ?_)) (contactsAt_le_slots S l hlen hn j hj)
      by_cases h : (candidatePair S i j && decide (m ≤ i) && decide (j < (walkFrom l).length) &&
        manhattan1 (siteAt (walkFrom l) i) (siteAt (walkFrom l) j)) = true
      · simp only [Bool.and_eq_true] at h
        simp [h]
      · simp [h]
    · refine Finset.sum_le_sum (fun i _ => ?_)
      by_cases h : (candidatePair S i j && decide (m ≤ i)) = true
      · simp only [Bool.and_eq_true] at h
        simp only [h, Bool.and_self, if_true, Bool.true_and]
        split_ifs <;> omega
      · simp only [Bool.and_eq_true, not_and] at h
        by_cases hc : candidatePair S i j = true
        · have := h hc
          simp [hc, this]
        · simp [hc]
  · refine le_of_eq (Finset.sum_eq_zero (fun i _ => ?_))
    by_cases hc : candidatePair S i j = true
    · have := candidatePair_lt hc
      have : ¬ m ≤ i := by omega
      simp [this]
    · simp [hc]

/-- **Degree-capped contact bounds** for every complete valid completion `l`
of a prefix `p`. -/
theorem contacts_le_degree_caps (S : List Bool) {p l : List (Fin 4)} (h : p <+: l)
    (hlen : S.length = l.length + 2) (hn : (walkFrom l).Nodup) :
    parityContactCount S (walkFrom l) ≤ knownContacts S p + futureCap S (p.length + 2) ∧
    parityContactCount S (walkFrom l) ≤
      knownContacts S p + crossContacts S (p.length + 2) (walkFrom l) +
        ffCap S (p.length + 2) := by
  rw [contacts_decompose S h]
  have h1 := cross_add_ff S (p.length + 2) (walkFrom l)
  have h2 := futureContacts_le_futureCap S l hlen hn (p.length + 2)
  have h3 := ffContacts_le_ffCap S l hlen hn (p.length + 2)
  omega

/-- The degree cap never exceeds the Python candidate-pair cap. -/
theorem futureCap_le_futurePairs (S : List Bool) (m : ℕ) :
    futureCap S m ≤ futurePairs S m := by
  rw [futurePairs, pairSum_comm, futureCap]
  refine Finset.sum_le_sum (fun j _ => ?_)
  split_ifs with hm
  · refine (min_le_right _ _).trans (Finset.sum_le_sum (fun i _ => ?_))
    by_cases hc : candidatePair S i j = true <;> simp [hc, hm]
  · exact Nat.zero_le _

end Frontier
end FoldBounds
