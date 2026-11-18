module

public import FoldBounds.FrontierKernel
public import FoldBounds.FrontierSoundness
public import FoldBounds.VerifiedFamilies
public import FoldBounds.VerifiedGeometry

/-!
# Soundness of the executable frontier kernel

`FoldBounds.Kernel` (module `FoldBounds/FrontierKernel.lean`) is an
executable, Mathlib-free checker for `FOLDBOUNDS-KERNEL 1` certificates.
This file proves it sound against the Lean HP lattice model, reusing the
semantic library of the frontier checker (`Frontier.ext`, the split
partition, the counting bound, degree and geometric contact caps, and the
partially directed and diagonal families).

Structure of the proof:

* `Rep S p st`: the incremental state `st` represents the direction prefix
  `p` (end site, reversed tagged walk, number of placed sites, known
  contacts).
* `stepState_none` / `stepState_some`: one step of the kernel agrees with
  appending a direction to the prefix; an occupied target is exactly a
  collision of the extended walk.
* `encEntry_sound`: a bound entry `(a, K, cnt, fam)` whose four claims hold
  semantically yields a sound enclosure for every `q ≥ 1` and every
  threshold; `boundEntry_ok` shows the recomputed entry has that property.
* `enumGo_sound`, `runTree_sound`: the enumerator and the certificate
  interpreter add a sound enclosure of the current subtree to the
  accumulator.
* `verify_sound`: every report returned by the string-level function
  `Kernel.verify` is sound for the echoed sequence and each echoed query.
-/
@[expose] public section

namespace FoldBounds
namespace Kernel

open LatticeModel Frontier FrontierCheck Finset

/-! ### The kernel primitives are the model primitives -/

theorem step_eq_advance (p : Site) (d : Fin 4) : step p d = advance p d := rfl

theorem manhattan1_eq (p r : Site) :
    Kernel.manhattan1 p r = LatticeModel.manhattan1 p r := rfl

theorem getD_toArray (S : List Bool) (i : ℕ) :
    S.toArray.getD i false = isHydrophobic S i := by
  simp [isHydrophobic]

theorem candPair_eq (S : List Bool) (i j : ℕ) :
    candPair S.toArray i j = candidatePair S i j := by
  simp [candPair, candidatePair, isHydrophobic]

/-! ### The reversed tagged walk -/

/-- Tag the sites of a walk with the hydrophobicity of their residues. -/
def tagFrom (S : List Bool) : ℕ → List Site → List (Site × Bool)
  | _, [] => []
  | i, s :: W => (s, isHydrophobic S i) :: tagFrom S (i + 1) W

theorem tagFrom_append (S : List Bool) (i : ℕ) (W V : List Site) :
    tagFrom S i (W ++ V) = tagFrom S i W ++ tagFrom S (i + W.length) V := by
  induction W generalizing i with
  | nil => simp [tagFrom]
  | cons s W ih =>
      simp only [List.cons_append, tagFrom, ih, List.length_cons]
      rw [show i + 1 + W.length = i + (W.length + 1) by omega]

theorem mem_tagFrom (S : List Bool) (i : ℕ) (W : List Site) (s : Site) :
    (∃ e ∈ tagFrom S i W, e.1 = s) ↔ s ∈ W := by
  induction W generalizing i with
  | nil => simp [tagFrom]
  | cons a W ih =>
      simp only [tagFrom, List.mem_cons, exists_eq_or_imp, ih]
      constructor
      · rintro (h | h)
        · exact Or.inl h.symm
        · exact Or.inr h
      · rintro (h | h)
        · exact Or.inl h.symm
        · exact Or.inr h

/-- The kernel's representation of a walk: most recent site first. -/
def revOf (S : List Bool) (W : List Site) : List (Site × Bool) := (tagFrom S 0 W).reverse

theorem revOf_snoc (S : List Bool) (W : List Site) (s : Site) :
    revOf S (W ++ [s]) = (s, isHydrophobic S W.length) :: revOf S W := by
  simp [revOf, tagFrom_append, tagFrom]

theorem mem_revOf (S : List Bool) (W : List Site) (s : Site) :
    (∃ e ∈ revOf S W, e.1 = s) ↔ s ∈ W := by
  simp only [revOf, List.mem_reverse]
  exact mem_tagFrom S 0 W s

theorem forall_mem_revOf (S : List Bool) (W : List Site) (P : Site → Prop) :
    (∀ e ∈ revOf S W, P e.1) ↔ ∀ s ∈ W, P s := by
  constructor
  · intro h s hs
    obtain ⟨e, he, rfl⟩ := (mem_revOf S W s).2 hs
    exact h e he
  · intro h e he
    exact h e.1 ((mem_revOf S W e.1).1 ⟨e, he, rfl⟩)

/-! ### The occupancy-and-contact scan -/

/-- Number of contacts counted by `scan` (specification). -/
def scanCount (s : Site) (hNew : Bool) : ℕ → List (Site × Bool) → ℕ
  | _, [] => 0
  | g, e :: rest =>
      (if hNew && e.2 && decide (3 ≤ g) && decide (g % 2 = 1) && Kernel.manhattan1 e.1 s
        then 1 else 0) + scanCount s hNew (g + 1) rest

theorem scan_none_iff (s : Site) (hNew : Bool) (g : ℕ) (R : List (Site × Bool)) (acc : ℕ) :
    scan s hNew g R acc = none ↔ ∃ e ∈ R, e.1 = s := by
  induction R generalizing g acc with
  | nil => simp [scan]
  | cons e rest ih =>
      by_cases he : e.1 = s
      · simp [scan, he]
      · simp only [scan, he, if_false, List.mem_cons, exists_eq_or_imp, false_or]
        exact ih _ _

theorem scan_some (s : Site) (hNew : Bool) (g : ℕ) (R : List (Site × Bool)) (acc c : ℕ)
    (h : scan s hNew g R acc = some c) : c = acc + scanCount s hNew g R := by
  induction R generalizing g acc with
  | nil => simp [scan] at h; simp [scanCount, h]
  | cons e rest ih =>
      by_cases he : e.1 = s
      · simp [scan, he] at h
      · simp only [scan, he, if_false] at h
        have := ih _ _ h
        rw [this, scanCount]
        split_ifs <;> omega

/-- Appending a site to the placed walk adds its own cross term. -/
theorem crossAt_snoc (S : List Bool) (W : List Site) (a : Site) (j : ℕ) (s : Site)
    (hW : W.length < S.length) :
    Transfer.crossAt S (W ++ [a]) (W.length + 1) j s =
      Transfer.crossAt S W W.length j s +
        (if (candidatePair S W.length j && LatticeModel.manhattan1 a s) = true then 1 else 0) := by
  unfold Transfer.crossAt
  have hsingle : (if (candidatePair S W.length j && LatticeModel.manhattan1 a s) = true
      then 1 else 0) = ∑ i ∈ range S.length, (if i = W.length then
        (if (candidatePair S W.length j && LatticeModel.manhattan1 a s) = true then 1 else 0)
        else 0) := by
    rw [Finset.sum_ite_eq' (range S.length) W.length]
    simp [hW]
  rw [hsingle, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rcases Nat.lt_trichotomy i W.length with hi | hi | hi
  · have h1 : siteAt (W ++ [a]) i = siteAt W i := by
      simp [siteAt, List.getElem?_append_left hi]
    have h2 : i < W.length + 1 := by omega
    have h3 : i ≠ W.length := by omega
    simp [h1, hi, h2, h3]
  · subst hi
    have h1 : siteAt (W ++ [a]) W.length = a := by simp [siteAt]
    simp [h1]
  · have h2 : ¬ i < W.length + 1 := by omega
    have h3 : ¬ i < W.length := by omega
    have h4 : i ≠ W.length := by omega
    simp [h2, h3, h4]

/-- The scan over the reversed tagged walk counts exactly the cross term of
the model. -/
theorem scanCount_revOf (S : List Bool) (s : Site) (j : ℕ) (W : List Site)
    (hW : W.length ≤ S.length) (hj : W.length ≤ j) :
    scanCount s (isHydrophobic S j) (j + 1 - W.length) (revOf S W) =
      Transfer.crossAt S W W.length j s := by
  induction W using List.reverseRecOn with
  | nil =>
      simp [revOf, tagFrom, scanCount, Transfer.crossAt]
  | append_singleton W a ih =>
      have hlen : (W ++ [a]).length = W.length + 1 := by simp
      rw [hlen] at hW hj ⊢
      rw [revOf_snoc, scanCount, crossAt_snoc S W a j s (by omega)]
      have hg : j + 1 - (W.length + 1) + 1 = j + 1 - W.length := by omega
      rw [hg, ih (by omega) (by omega), add_comm]
      congr 1
      have hcond : (isHydrophobic S j && isHydrophobic S W.length &&
          decide (3 ≤ j + 1 - (W.length + 1)) && decide ((j + 1 - (W.length + 1)) % 2 = 1) &&
          Kernel.manhattan1 a s) = (candidatePair S W.length j && LatticeModel.manhattan1 a s) := by
        have e1 : j + 1 - (W.length + 1) = j - W.length := by omega
        rw [e1, manhattan1_eq]
        unfold candidatePair
        have e2 : decide (3 ≤ j - W.length) = decide (W.length + 3 ≤ j) := by
          apply Bool.decide_congr; omega
        rw [e2]
        cases isHydrophobic S j <;> cases isHydrophobic S W.length <;>
          cases decide (W.length + 3 ≤ j) <;> cases decide ((j - W.length) % 2 = 1) <;> simp
      simp only [hcond]

/-! ### Known contacts after one more step -/

theorem ffContacts_snoc_zero (S : List Bool) (W : List Site) (s : Site) :
    ffContacts S W.length (W ++ [s]) = 0 := by
  unfold ffContacts pairSum
  refine Finset.sum_eq_zero (fun i _ => Finset.sum_eq_zero (fun j _ => ?_))
  by_cases hc : candidatePair S i j = true
  · have := candidatePair_lt hc
    have h2 : ¬ (W.length ≤ i ∧ j < W.length + 1) := by omega
    simp only [List.length_append, List.length_singleton, Bool.and_eq_true, decide_eq_true_eq,
      ite_eq_right_iff, one_ne_zero, imp_false]
    tauto
  · simp [hc]

/-- Placing one more residue adds exactly its contacts with the placed walk. -/
theorem knownContacts_snoc (S : List Bool) (p : List (Fin 4)) (d : Fin 4)
    (hlt : p.length + 2 < S.length) :
    knownContacts S (p ++ [d]) =
      knownContacts S p +
        Transfer.crossAt S (walkFrom p) (p.length + 2) (p.length + 2)
          (advance (endFrom first p) d) := by
  have hlen : (walkFrom p).length = p.length + 2 := walkFrom_length p
  have hdec := contacts_decompose S (List.prefix_append p [d])
  rw [walkFrom_append_single] at hdec
  have hcross := Transfer.crossContacts_snoc S p (walkFrom p) (advance (endFrom first p) d)
    (List.prefix_refl _) (by omega)
  rw [Transfer.crossContacts_self, zero_add, hlen] at hcross
  have hff := ffContacts_snoc_zero S (walkFrom p) (advance (endFrom first p) d)
  rw [hlen] at hff
  unfold knownContacts
  rw [walkFrom_append_single, hdec, hcross, hff]
  rfl

/-! ### The representation invariant -/

/-- The kernel state `st` represents the direction prefix `p`. -/
structure Rep (S : List Bool) (p : List (Fin 4)) (st : State) : Prop where
  cur : st.cur = endFrom first p
  rev : st.rev = revOf S (walkFrom p)
  placed : st.placed = p.length + 2
  known : st.known = knownContacts S p

theorem rep_root (S : List Bool) : Rep S [] (rootState S.toArray) := by
  refine ⟨rfl, ?_, rfl, ?_⟩
  · simp [rootState, revOf, tagFrom, walkFrom, trace, isHydrophobic, site0, site1,
      LatticeModel.first, LatticeModel.origin]
  · simp [rootState, knownContacts_nil]

theorem nodup_snoc_iff (p : List (Fin 4)) (d : Fin 4) :
    (walkFrom (p ++ [d])).Nodup ↔
      (walkFrom p).Nodup ∧ advance (endFrom first p) d ∉ walkFrom p := by
  rw [walkFrom_append_single, List.nodup_append]
  simp only [List.nodup_singleton, List.mem_singleton, true_and]
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨h1, fun hm => h2 _ hm _ rfl rfl⟩
  · rintro ⟨h1, h2⟩
    refine ⟨h1, ?_⟩
    rintro a ha b rfl rfl
    exact h2 ha

/-- The kernel refuses a step exactly when the target site is occupied. -/
theorem stepState_none {S : List Bool} {p : List (Fin 4)} {st : State} (hrep : Rep S p st)
    (d : Fin 4) :
    stepState S.toArray st d = none ↔ advance (endFrom first p) d ∈ walkFrom p := by
  unfold stepState
  simp only
  rw [← mem_revOf S (walkFrom p), ← hrep.rev, ← hrep.cur, ← step_eq_advance,
    ← scan_none_iff (step st.cur d) (S.toArray.getD st.placed false) 1 st.rev 0]
  cases scan (step st.cur d) (S.toArray.getD st.placed false) 1 st.rev 0 <;> simp

/-- An accepted step yields the state of the extended prefix. -/
theorem stepState_some {S : List Bool} {p : List (Fin 4)} {st st' : State} (hrep : Rep S p st)
    (hlt : p.length + 2 < S.length) {d : Fin 4} (h : stepState S.toArray st d = some st') :
    Rep S (p ++ [d]) st' ∧ advance (endFrom first p) d ∉ walkFrom p := by
  have hnot : advance (endFrom first p) d ∉ walkFrom p := by
    intro hm
    rw [(stepState_none hrep d).2 hm] at h
    exact absurd h (by simp)
  refine ⟨?_, hnot⟩
  unfold stepState at h
  simp only at h
  cases hs : scan (step st.cur d) (S.toArray.getD st.placed false) 1 st.rev 0 with
  | none => rw [hs] at h; exact absurd h (by simp)
  | some c =>
      rw [hs] at h
      simp only [Option.some.injEq] at h
      subst h
      have hc := scan_some _ _ _ _ _ _ hs
      have hlen : (walkFrom p).length = p.length + 2 := walkFrom_length p
      rw [zero_add, hrep.rev, hrep.placed, getD_toArray, hrep.cur, step_eq_advance] at hc
      have hcount := scanCount_revOf S (advance (endFrom first p) d) (p.length + 2) (walkFrom p)
        (by omega) (by omega)
      rw [hlen] at hcount
      have hg : p.length + 2 + 1 - (p.length + 2) = 1 := by omega
      rw [hg] at hcount
      refine ⟨?_, ?_, ?_, ?_⟩
      · simp [hrep.cur, step_eq_advance, endFrom_append, endFrom]
      · simp only
        rw [walkFrom_append_single, revOf_snoc, hlen, hrep.rev, hrep.placed, getD_toArray,
          hrep.cur, step_eq_advance]
      · simp [hrep.placed]
      · simp only
        rw [knownContacts_snoc S p d hlt, hc, hcount, hrep.known]

/-! ### Enclosure algebra -/

theorem enc_ext {a b : Enclosure} (h1 : a.zLo = b.zLo) (h2 : a.zHi = b.zHi) (h3 : a.eLo = b.eLo)
    (h4 : a.eHi = b.eHi) (h5 : a.nLo = b.nLo) (h6 : a.nHi = b.nHi) : a = b := by
  cases a; cases b; simp_all

@[simp] theorem add_zLo (a b : Enclosure) : (a + b).zLo = a.zLo + b.zLo := rfl
@[simp] theorem add_zHi (a b : Enclosure) : (a + b).zHi = a.zHi + b.zHi := rfl
@[simp] theorem add_eLo (a b : Enclosure) : (a + b).eLo = a.eLo + b.eLo := rfl
@[simp] theorem add_eHi (a b : Enclosure) : (a + b).eHi = a.eHi + b.eHi := rfl
@[simp] theorem add_nLo (a b : Enclosure) : (a + b).nLo = a.nLo + b.nLo := rfl
@[simp] theorem add_nHi (a b : Enclosure) : (a + b).nHi = a.nHi + b.nHi := rfl
@[simp] theorem zero_zLo : (0 : Enclosure).zLo = 0 := rfl
@[simp] theorem zero_zHi : (0 : Enclosure).zHi = 0 := rfl
@[simp] theorem zero_eLo : (0 : Enclosure).eLo = 0 := rfl
@[simp] theorem zero_eHi : (0 : Enclosure).eHi = 0 := rfl
@[simp] theorem zero_nLo : (0 : Enclosure).nLo = 0 := rfl
@[simp] theorem zero_nHi : (0 : Enclosure).nHi = 0 := rfl

theorem enc_add_zero (a : Enclosure) : a + 0 = a := by apply enc_ext <;> simp

theorem enc_zero_add (a : Enclosure) : 0 + a = a := by apply enc_ext <;> simp

theorem enc_add_assoc (a b c : Enclosure) : a + b + c = a + (b + c) := by
  apply enc_ext <;> simp [add_assoc]

/-! ### Evaluation of the accumulator -/

/-- The six accumulated endpoints for `(q, t)`, as an `Enclosure`. -/
def encOf (q : ℚ) (t : ℕ) (ex : Array ℕ) (es : List Entry) : Enclosure :=
  ⟨histE q t ex + histN q t ex + (es.map (entZLo q)).sum,
   histE q t ex + histN q t ex + (es.map (entZHi q)).sum,
   histE q t ex + (es.map (entELo q t)).sum,
   histE q t ex + (es.map (entEHi q t)).sum,
   histN q t ex + (es.map (entNLo q t)).sum,
   histN q t ex + (es.map (entNHi q t)).sum⟩

/-- Enclosure of one conformation with exactly `k` contacts. -/
def encExact (q : ℚ) (t k : ℕ) : Enclosure :=
  if t ≤ k then ⟨q ^ k, q ^ k, q ^ k, q ^ k, 0, 0⟩ else ⟨q ^ k, q ^ k, 0, 0, q ^ k, q ^ k⟩

/-- Enclosure contributed by one bound entry. -/
def encEntry (q : ℚ) (t : ℕ) (e : Entry) : Enclosure :=
  ⟨entZLo q e, entZHi q e, entELo q t e, entEHi q t e, entNLo q t e, entNHi q t e⟩

theorem list_sum_range {M : Type*} [AddCommMonoid M] (n : ℕ) (g : ℕ → M) :
    ((List.range n).map g).sum = ∑ i ∈ Finset.range n, g i := by
  induction n with
  | zero => simp
  | succ n ih => rw [List.sum_range_succ, Finset.sum_range_succ, ih]

/-- Adding one conformation with `k` contacts changes a weighted histogram
sum by exactly the weight of `k`. -/
theorem hist_bump (w : ℕ → ℚ) {ex ex' : Array ℕ} {k : ℕ} (h : bump ex k = some ex') :
    ((List.range ex'.size).map (fun i => (ex'.getD i 0 : ℚ) * w i)).sum =
      ((List.range ex.size).map (fun i => (ex.getD i 0 : ℚ) * w i)).sum + w k := by
  unfold bump at h
  split at h
  · rename_i hk
    simp only [Option.some.injEq] at h
    subst h
    rw [list_sum_range, list_sum_range, Array.size_set]
    have hsingle : w k = ∑ i ∈ Finset.range ex.size, (if i = k then w k else 0) := by
      rw [Finset.sum_ite_eq' (Finset.range ex.size) k]
      simp [hk]
    rw [hsingle, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    rw [Finset.mem_range] at hi
    by_cases hik : i = k
    · subst hik
      simp [Array.getD_eq_getD_getElem?, hi]
      ring
    · have hki : k ≠ i := fun h => hik h.symm
      simp [Array.getD_eq_getD_getElem?, hki, hik]
  · exact absurd h (by simp)

theorem histE_eq (q : ℚ) (t : ℕ) (ex : Array ℕ) :
    histE q t ex = ((List.range ex.size).map
      (fun i => (ex.getD i 0 : ℚ) * (if t ≤ i then q ^ i else 0))).sum := by
  unfold histE
  congr 1
  apply List.map_congr_left
  intro i _
  split_ifs <;> simp

theorem histN_eq (q : ℚ) (t : ℕ) (ex : Array ℕ) :
    histN q t ex = ((List.range ex.size).map
      (fun i => (ex.getD i 0 : ℚ) * (if i < t then q ^ i else 0))).sum := by
  unfold histN
  congr 1
  apply List.map_congr_left
  intro i _
  split_ifs <;> simp

theorem histE_bump (q : ℚ) (t : ℕ) {ex ex' : Array ℕ} {k : ℕ} (h : bump ex k = some ex') :
    histE q t ex' = histE q t ex + (if t ≤ k then q ^ k else 0) := by
  rw [histE_eq, histE_eq, hist_bump (fun i => if t ≤ i then q ^ i else 0) h]

theorem histN_bump (q : ℚ) (t : ℕ) {ex ex' : Array ℕ} {k : ℕ} (h : bump ex k = some ex') :
    histN q t ex' = histN q t ex + (if k < t then q ^ k else 0) := by
  rw [histN_eq, histN_eq, hist_bump (fun i => if i < t then q ^ i else 0) h]

theorem encOf_bump (q : ℚ) (t : ℕ) (es : List Entry) {ex ex' : Array ℕ} {k : ℕ}
    (h : bump ex k = some ex') :
    encOf q t ex' es = encOf q t ex es + encExact q t k := by
  have hE := histE_bump q t h
  have hN := histN_bump q t h
  by_cases ht : t ≤ k
  · have hnt : ¬ k < t := by omega
    rw [if_pos ht] at hE
    rw [if_neg hnt] at hN
    apply enc_ext <;> simp [encOf, encExact, ht, hE, hN] <;> ring
  · have hnt : k < t := by omega
    rw [if_neg ht] at hE
    rw [if_pos hnt] at hN
    apply enc_ext <;> simp [encOf, encExact, ht, hE, hN] <;> ring

theorem encOf_cons (q : ℚ) (t : ℕ) (ex : Array ℕ) (e : Entry) (es : List Entry) :
    encOf q t ex (e :: es) = encOf q t ex es + encEntry q t e := by
  apply enc_ext <;> simp [encOf, encEntry] <;> ring

theorem encOf_empty (q : ℚ) (t : ℕ) (N : ℕ) : encOf q t (Array.replicate N 0) [] = 0 := by
  have hz : ∀ i, ((Array.replicate N 0 : Array ℕ).getD i 0 : ℚ) = 0 := by
    intro i
    simp only [Array.getD_eq_getD_getElem?, Array.getElem?_replicate]
    split_ifs <;> simp
  have hE : histE q t (Array.replicate N 0) = 0 := by
    rw [histE_eq, list_sum_range]
    exact Finset.sum_eq_zero (fun i _ => by rw [hz]; simp)
  have hN : histN q t (Array.replicate N 0) = 0 := by
    rw [histN_eq, list_sum_range]
    exact Finset.sum_eq_zero (fun i _ => by rw [hz]; simp)
  apply enc_ext <;> simp [encOf, hE, hN]

/-! ### Bound entries: semantic validity implies a sound enclosure -/

/-- The four claims of a bound entry at prefix `p`. -/
def EntryOK (steps : ℕ) (S : List Bool) (p : List (Fin 4)) (e : Entry) : Prop :=
  (∀ l ∈ ext steps p, e.a ≤ parityContactCount S (walkFrom l) ∧
      parityContactCount S (walkFrom l) ≤ e.K) ∧
  (ext steps p).card ≤ e.cnt ∧ e.fam ≤ (ext steps p).card

theorem sum_le_cnt_mul {T : Finset (List (Fin 4))} (G : List (Fin 4) → ℚ) (C : ℚ) (cnt : ℕ)
    (hC : 0 ≤ C) (hG : ∀ l ∈ T, G l ≤ C) (hcard : T.card ≤ cnt) :
    ∑ l ∈ T, G l ≤ (cnt : ℚ) * C := by
  calc ∑ l ∈ T, G l ≤ T.card • C := Finset.sum_le_card_nsmul _ _ _ hG
    _ = (T.card : ℚ) * C := by rw [nsmul_eq_mul]
    _ ≤ (cnt : ℚ) * C := by
        apply mul_le_mul_of_nonneg_right _ hC
        exact_mod_cast hcard

theorem fam_mul_le_sum {T : Finset (List (Fin 4))} (F : List (Fin 4) → ℚ) (c : ℚ) (fam : ℕ)
    (hc : 0 ≤ c) (hF : ∀ l ∈ T, c ≤ F l) (hcard : fam ≤ T.card) :
    (fam : ℚ) * c ≤ ∑ l ∈ T, F l := by
  calc (fam : ℚ) * c ≤ (T.card : ℚ) * c := by
        apply mul_le_mul_of_nonneg_right _ hc
        exact_mod_cast hcard
    _ = T.card • c := by rw [nsmul_eq_mul]
    _ ≤ ∑ l ∈ T, F l := Finset.card_nsmul_le_sum _ _ _ hF

/-- **Bound-entry soundness, uniformly in `q` and the threshold.** If every
completion of `p` has between `a` and `K` contacts and the number of
completions is between `fam` and `cnt`, then the entry's six endpoints
enclose the subtree masses for every `q ≥ 1` and every threshold `t`.  The
complement upper endpoint uses the exponent `min K (t - 1)`. -/
theorem encEntry_sound {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) {q : ℚ}
    (hq : 1 ≤ q) (t : ℕ) {p : List (Fin 4)} {e : Entry} (h : EntryOK steps S p e) :
    Encloses steps S q t p (encEntry q t e) := by
  obtain ⟨hk, hcnt, hfam⟩ := h
  have hq0 : (0 : ℚ) ≤ q := by linarith
  have hpos : ∀ k : ℕ, (0 : ℚ) ≤ q ^ k := fun k => pow_nonneg hq0 k
  have mZ := massZ_eq steps S hS q p
  have mA := massA_eq steps S hS q t p
  have mN := massN_eq steps S hS q t p
  have hlo : ∀ l ∈ ext steps p, q ^ e.a ≤ q ^ parityContactCount S (walkFrom l) :=
    fun l hl => qpow_mono hq (hk l hl).1
  have hhi : ∀ l ∈ ext steps p, q ^ parityContactCount S (walkFrom l) ≤ q ^ e.K :=
    fun l hl => qpow_mono hq (hk l hl).2
  have hfamZ : (e.fam : ℚ) * q ^ e.a ≤
      ∑ l ∈ ext steps p, q ^ parityContactCount S (walkFrom l) :=
    fam_mul_le_sum _ _ _ (hpos _) hlo hfam
  have hcntZ : ∑ l ∈ ext steps p, q ^ parityContactCount S (walkFrom l) ≤
      (e.cnt : ℚ) * q ^ e.K :=
    sum_le_cnt_mul _ _ _ (hpos _) hhi hcnt
  have hA0 : (0 : ℚ) ≤ ∑ l ∈ ext steps p,
      (if t ≤ parityContactCount S (walkFrom l) then q ^ parityContactCount S (walkFrom l)
        else 0) :=
    Finset.sum_nonneg (fun l _ => by split_ifs <;> simp [hpos])
  have hN0 : (0 : ℚ) ≤ ∑ l ∈ ext steps p,
      (if parityContactCount S (walkFrom l) < t then q ^ parityContactCount S (walkFrom l)
        else 0) :=
    Finset.sum_nonneg (fun l _ => by split_ifs <;> simp [hpos])
  rw [Encloses, mZ, mA, mN]
  refine ⟨hfamZ, hcntZ, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · show entELo q t e ≤ _
    unfold entELo
    split_ifs with ht
    · have hsum : ∑ l ∈ ext steps p,
          (if t ≤ parityContactCount S (walkFrom l) then
            q ^ parityContactCount S (walkFrom l) else 0) =
          ∑ l ∈ ext steps p, q ^ parityContactCount S (walkFrom l) :=
        Finset.sum_congr rfl (fun l hl => by rw [if_pos (ht.trans (hk l hl).1)])
      rw [hsum]
      exact hfamZ
    · exact hA0
  · show _ ≤ entEHi q t e
    unfold entEHi
    split_ifs with ht
    · refine sum_le_cnt_mul _ _ _ (hpos _) (fun l hl => ?_) hcnt
      split_ifs
      · exact hhi l hl
      · exact hpos _
    · refine le_of_eq (Finset.sum_eq_zero (fun l hl => ?_))
      have := (hk l hl).2
      rw [if_neg (by omega)]
  · show entNLo q t e ≤ _
    unfold entNLo
    split_ifs with ht
    · have hsum : ∑ l ∈ ext steps p,
          (if parityContactCount S (walkFrom l) < t then
            q ^ parityContactCount S (walkFrom l) else 0) =
          ∑ l ∈ ext steps p, q ^ parityContactCount S (walkFrom l) :=
        Finset.sum_congr rfl (fun l hl => by rw [if_pos (lt_of_le_of_lt (hk l hl).2 ht)])
      rw [hsum]
      exact hfamZ
    · exact hN0
  · show _ ≤ entNHi q t e
    unfold entNHi
    split_ifs with ht
    · refine sum_le_cnt_mul _ _ _ (hpos _) (fun l hl => ?_) hcnt
      split_ifs with h1
      · have := (hk l hl).2
        exact qpow_mono hq (le_min this (by omega))
      · exact hpos _
    · refine le_of_eq (Finset.sum_eq_zero (fun l hl => ?_))
      have := (hk l hl).1
      rw [if_neg (by omega)]
  · show (0 : ℚ) ≤ entZLo q e
    unfold entZLo
    exact mul_nonneg (Nat.cast_nonneg _) (hpos _)
  · show (0 : ℚ) ≤ entELo q t e
    unfold entELo
    split_ifs
    · exact mul_nonneg (Nat.cast_nonneg _) (hpos _)
    · exact le_refl _
  · show (0 : ℚ) ≤ entNLo q t e
    unfold entNLo
    split_ifs
    · exact mul_nonneg (Nat.cast_nonneg _) (hpos _)
    · exact le_refl _

/-! ### The recomputed bound entry is semantically valid -/

/-- The cap table returns the model's degree cap. -/
def CapTabOK (S : List Bool) (capTab : Array ℕ) : Prop :=
  ∀ m fc, capTab[m]? = some fc → fc = Frontier.futureCap S m

theorem futureCap_eq (S : List Bool) (m : ℕ) :
    Kernel.futureCap S.toArray m = Frontier.futureCap S m := by
  unfold Kernel.futureCap Frontier.futureCap
  rw [list_sum_range]
  simp only [List.size_toArray]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  have hc : candCount S.toArray j =
      ∑ i ∈ Finset.range S.length, if candidatePair S i j = true then 1 else 0 := by
    unfold candCount
    rw [list_sum_range]
    simp only [List.size_toArray, candPair_eq]
  rw [hc]
  rfl

theorem capTable_ok (S : List Bool) : CapTabOK S (capTable S.toArray) := by
  intro m fc h
  unfold capTable at h
  rw [Array.getElem?_ofFn] at h
  split at h
  · simp only [Option.some.injEq] at h
    rw [← h]
    exact futureCap_eq S m
  · exact absurd h (by simp)

theorem pdTriple_eq (r : ℕ) : Kernel.pdTriple r = Families.pdTriple r := by
  induction r with
  | zero => rfl
  | succ r ih => simp only [Kernel.pdTriple, Families.pdTriple, ih]

theorem pdCount_eq (r : ℕ) : Kernel.pdCount r = Families.pdCount r := by
  unfold Kernel.pdCount Families.pdCount
  rw [pdTriple_eq]

theorem pdTest_imp {S : List Bool} {p : List (Fin 4)} {st : State} (hrep : Rep S p st)
    {m : Fin 4} (h : Kernel.pdTest m st.cur st.rev = true) : Families.pdOK m p = true := by
  unfold Families.pdOK
  unfold Kernel.pdTest at h
  rw [List.all_eq_true] at h ⊢
  rw [hrep.rev, hrep.cur] at h
  exact (forall_mem_revOf S (walkFrom p) (fun s =>
    (s == endFrom first p ||
      decide (Families.oF m s < Families.oF m (endFrom first p))) = true)).1 h

theorem diagTest_imp {S : List Bool} {p : List (Fin 4)} {st : State} (hrep : Rep S p st)
    {a b : Bool} (h : Kernel.diagTest a b st.cur st.rev = true) : Frontier.diagOK a b p = true := by
  unfold Frontier.diagOK
  unfold Kernel.diagTest at h
  rw [List.all_eq_true] at h ⊢
  rw [hrep.rev, hrep.cur] at h
  exact (forall_mem_revOf S (walkFrom p) (fun s =>
    decide (Frontier.pot a b s ≤ Frontier.pot a b (endFrom first p)) = true)).1 h

/-- The certified family count never exceeds the number of completions. -/
theorem famCount_le {steps : ℕ} {S : List Bool} {p : List (Fin 4)} {st : State}
    (hrep : Rep S p st) (hn : (walkFrom p).Nodup) (hle : p.length ≤ steps) :
    famCount st (steps - p.length) ≤ (ext steps p).card := by
  unfold famCount
  split_ifs with h1 h2
  · have key : ∀ m, Kernel.pdTest m st.cur st.rev = true →
        Kernel.pdCount (steps - p.length) ≤ (ext steps p).card := by
      intro m hm
      rw [pdCount_eq, ← Families.card_pdFamily m p]
      exact Finset.card_le_card (Families.pdFamily_subset_ext hn (pdTest_imp hrep hm) hle)
    simp only [Bool.or_eq_true] at h1
    rcases h1 with ((h | h) | h) | h <;> exact key _ h
  · have key : ∀ a b, Kernel.diagTest a b st.cur st.rev = true →
        2 ^ (steps - p.length) ≤ (ext steps p).card := by
      intro a b hab
      rw [← Frontier.card_family a b p]
      exact Finset.card_le_card (Frontier.family_subset_ext hn (diagTest_imp hrep hab) hle)
    simp only [Bool.or_eq_true] at h2
    rcases h2 with ((h | h) | h) | h <;> exact key _ _ h
  · exact Nat.zero_le _

theorem minX_mem : ∀ (R : List (Site × Bool)), R ≠ [] → ∃ e ∈ R, e.1.1 = minX R
  | [], h => absurd rfl h
  | [e], _ => ⟨e, by simp, rfl⟩
  | e :: e' :: rest, _ => by
      obtain ⟨x, hx, hxe⟩ := minX_mem (e' :: rest) (by simp)
      simp only [minX]
      rcases min_choice e.1.1 (minX (e' :: rest)) with h | h
      · exact ⟨e, by simp, h.symm⟩
      · exact ⟨x, List.mem_cons_of_mem _ hx, by rw [h, hxe]⟩

theorem maxX_mem : ∀ (R : List (Site × Bool)), R ≠ [] → ∃ e ∈ R, e.1.1 = maxX R
  | [], h => absurd rfl h
  | [e], _ => ⟨e, by simp, rfl⟩
  | e :: e' :: rest, _ => by
      obtain ⟨x, hx, hxe⟩ := maxX_mem (e' :: rest) (by simp)
      simp only [maxX]
      rcases max_choice e.1.1 (maxX (e' :: rest)) with h | h
      · exact ⟨e, by simp, h.symm⟩
      · exact ⟨x, List.mem_cons_of_mem _ hx, by rw [h, hxe]⟩

theorem minY_mem : ∀ (R : List (Site × Bool)), R ≠ [] → ∃ e ∈ R, e.1.2 = minY R
  | [], h => absurd rfl h
  | [e], _ => ⟨e, by simp, rfl⟩
  | e :: e' :: rest, _ => by
      obtain ⟨x, hx, hxe⟩ := minY_mem (e' :: rest) (by simp)
      simp only [minY]
      rcases min_choice e.1.2 (minY (e' :: rest)) with h | h
      · exact ⟨e, by simp, h.symm⟩
      · exact ⟨x, List.mem_cons_of_mem _ hx, by rw [h, hxe]⟩

theorem maxY_mem : ∀ (R : List (Site × Bool)), R ≠ [] → ∃ e ∈ R, e.1.2 = maxY R
  | [], h => absurd rfl h
  | [e], _ => ⟨e, by simp, rfl⟩
  | e :: e' :: rest, _ => by
      obtain ⟨x, hx, hxe⟩ := maxY_mem (e' :: rest) (by simp)
      simp only [maxY]
      rcases max_choice e.1.2 (maxY (e' :: rest)) with h | h
      · exact ⟨e, by simp, h.symm⟩
      · exact ⟨x, List.mem_cons_of_mem _ hx, by rw [h, hxe]⟩

theorem revOf_ne_nil (S : List Bool) {W : List Site} (h : W ≠ []) : revOf S W ≠ [] := by
  cases W with
  | nil => exact absurd rfl h
  | cons a W => simp [revOf, tagFrom]

/-- **Geometric cap from the bounding box of the placed sites.** Every
completion of `p` has at most `n + 1 - minPerim n w h` contacts, where
`w × h` is the bounding box of the placed walk. -/
theorem contacts_le_geom {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2)
    {p l : List (Fin 4)} (hl : l ∈ ext steps p) :
    parityContactCount S (walkFrom l) ≤
      S.length + 1 - Geometry.minPerim S.length
        (maxX (revOf S (walkFrom p)) - minX (revOf S (walkFrom p)) + 1).toNat
        (maxY (revOf S (walkFrom p)) - minY (revOf S (walkFrom p)) + 1).toNat := by
  obtain ⟨hlen, hn, hpre⟩ := mem_ext.1 hl
  have hWlen : (walkFrom l).length = steps + 2 := by rw [walkFrom_length, hlen]
  have hc : adjacentChain (walkFrom l) = true := walkFrom_adjacentChain l
  have hsub : ∀ a ∈ walkFrom p, a ∈ walkFrom l := fun a ha => (walkFrom_prefix hpre).subset ha
  have hne : walkFrom p ≠ [] := by
    intro h0
    have := walkFrom_length p
    rw [h0] at this
    simp at this
  have hR := revOf_ne_nil S hne
  obtain ⟨a, ha, hax⟩ := minX_mem _ hR
  obtain ⟨b, hb, hbx⟩ := maxX_mem _ hR
  obtain ⟨c, hc', hcy⟩ := minY_mem _ hR
  obtain ⟨d, hd, hdy⟩ := maxY_mem _ hR
  have hmem : ∀ e ∈ revOf S (walkFrom p), e.1 ∈ walkFrom l :=
    fun e he => hsub _ ((mem_revOf S (walkFrom p) e.1).1 ⟨e, he, rfl⟩)
  have hw : (maxX (revOf S (walkFrom p)) - minX (revOf S (walkFrom p)) + 1).toNat ≤
      (Geometry.cols (walkFrom l)).card := by
    rw [← hax, ← hbx]
    exact Geometry.cols_ge_span hc (hmem a ha) (hmem b hb)
  have hh : (maxY (revOf S (walkFrom p)) - minY (revOf S (walkFrom p)) + 1).toNat ≤
      (Geometry.rows (walkFrom l)).card := by
    rw [← hcy, ← hdy]
    exact Geometry.rows_ge_span hc (hmem c hc') (hmem d hd)
  have h1 := Geometry.parityContactCount_le_contactPairs S (walkFrom l) (by omega)
  have h2 := Geometry.contacts_le_rows_cols hn hc (by omega)
  have hcols : 0 < (Geometry.cols (walkFrom l)).card :=
    Finset.card_pos.2 ⟨a.1.1, Finset.mem_image.2 ⟨a.1, List.mem_toFinset.2 (hmem a ha), rfl⟩⟩
  have hp := Geometry.minPerim_le (walkFrom l).length _ _ (Geometry.cols (walkFrom l)).card
    (Geometry.rows (walkFrom l)).card hw Geometry.cols_le_length hh
    (Geometry.length_le_rows_mul_cols hn) hcols
  rw [hWlen] at h2 hp
  rw [hS]
  omega

/-- The entry recomputed by the kernel satisfies its four semantic claims. -/
theorem boundEntry_ok {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2)
    {capTab : Array ℕ} (hcap : CapTabOK S capTab) {p : List (Fin 4)} {st : State}
    (hrep : Rep S p st) (hn : (walkFrom p).Nodup) (hle : p.length ≤ steps) {e : Entry}
    (h : boundEntry S.toArray capTab st = some e) : EntryOK steps S p e := by
  unfold boundEntry at h
  cases hfc : capTab[st.placed]? with
  | none => rw [hfc] at h; exact absurd h (by simp)
  | some fc =>
      rw [hfc] at h
      simp only [Option.some.injEq] at h
      subst h
      have hfc' := hcap _ _ hfc
      have hr : steps + 2 - (p.length + 2) = steps - p.length := by omega
      refine ⟨fun l hl => ⟨?_, ?_⟩, ?_, ?_⟩
      · simp only
        rw [hrep.known]
        exact knownContacts_le S (mem_ext.1 hl).2.2
      · simp only [List.size_toArray]
        refine le_min ?_ ?_
        · have hl' := mem_ext.1 hl
          have := (contacts_le_degree_caps S hl'.2.2 (by omega) hl'.2.1).1
          rw [hrep.known, hfc', hrep.placed]
          exact this
        · have := contacts_le_geom hS hl
          rw [hrep.rev]
          exact this
      · simp only [List.size_toArray]
        rw [hrep.placed, hS, hr]
        exact card_ext_le_pow_three steps p
      · simp only [List.size_toArray]
        rw [hrep.placed, hS, hr]
        exact famCount_le hrep hn hle

/-! ### Exhaustive enumeration -/

/-- A direction whose target site is occupied has no completion. -/
theorem encloses_blocked {steps : ℕ} {S : List Bool} {q : ℚ} {t : ℕ} {p : List (Fin 4)}
    {d : Fin 4} (h : advance (endFrom first p) d ∈ walkFrom p) :
    Encloses steps S q t (p ++ [d]) 0 :=
  masses_of_empty (ext_eq_empty_of_not_nodup (fun hn => ((nodup_snoc_iff p d).1 hn).2 h))

theorem encExact_sound {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) {q : ℚ}
    (hq : 1 ≤ q) (t : ℕ) {p : List (Fin 4)} (hlen : p.length = steps) (hn : (walkFrom p).Nodup) :
    Encloses steps S q t p (encExact q t (knownContacts S p)) :=
  exactEnc_sound steps S hS q hq t p hlen hn

theorem childWith_sound {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) {q : ℚ}
    (t : ℕ) (es : List Entry) {r : ℕ} {k : State → Array ℕ → Option (Array ℕ)}
    (ih : ∀ (p : List (Fin 4)) (st : State) (ex ex' : Array ℕ), Rep S p st →
      (walkFrom p).Nodup → p.length + r = steps → k st ex = some ex' →
      ∃ e, Encloses steps S q t p e ∧ encOf q t ex' es = encOf q t ex es + e)
    {p : List (Fin 4)} {st : State} (hrep : Rep S p st) (hn : (walkFrom p).Nodup)
    (hlen : p.length + (r + 1) = steps) (d : Fin 4) {ex ex' : Array ℕ}
    (h : childWith S.toArray st d ex k = some ex') :
    ∃ e, Encloses steps S q t (p ++ [d]) e ∧ encOf q t ex' es = encOf q t ex es + e := by
  unfold childWith at h
  cases hs : stepState S.toArray st d with
  | none =>
      rw [hs] at h
      simp only [Option.some.injEq] at h
      subst h
      exact ⟨0, encloses_blocked ((stepState_none hrep d).1 hs), (enc_add_zero _).symm⟩
  | some st' =>
      rw [hs] at h
      simp only at h
      obtain ⟨hrep', hfree⟩ := stepState_some hrep (by omega) hs
      exact ih (p ++ [d]) st' ex ex' hrep' ((nodup_snoc_iff p d).2 ⟨hn, hfree⟩)
        (by simp; omega) h

/-- **Soundness of the enumerator.** The histogram returned by `enumGo`
differs from the input histogram by a sound (indeed exact) enclosure of the
subtree of `p`, for every `q ≥ 1` and every threshold. -/
theorem enumGo_sound {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) {q : ℚ}
    (hq : 1 ≤ q) (t : ℕ) (es : List Entry) :
    ∀ (r : ℕ) (p : List (Fin 4)) (st : State) (ex ex' : Array ℕ), Rep S p st →
      (walkFrom p).Nodup → p.length + r = steps → enumGo S.toArray r st ex = some ex' →
      ∃ e, Encloses steps S q t p e ∧ encOf q t ex' es = encOf q t ex es + e := by
  intro r
  induction r with
  | zero =>
      intro p st ex ex' hrep hn hlen h
      rw [enumGo] at h
      refine ⟨encExact q t st.known, ?_, encOf_bump q t es h⟩
      rw [hrep.known]
      exact encExact_sound hS hq t (by omega) hn
  | succ r ih =>
      intro p st ex ex' hrep hn hlen h
      rw [enumGo] at h
      cases h0 : childWith S.toArray st 0 ex (enumGo S.toArray r) with
      | none => rw [h0] at h; exact absurd h (by simp)
      | some ex0 =>
        rw [h0] at h
        simp only at h
        cases h1 : childWith S.toArray st 1 ex0 (enumGo S.toArray r) with
        | none => rw [h1] at h; exact absurd h (by simp)
        | some ex1 =>
          rw [h1] at h
          simp only at h
          cases h2 : childWith S.toArray st 2 ex1 (enumGo S.toArray r) with
          | none => rw [h2] at h; exact absurd h (by simp)
          | some ex2 =>
            rw [h2] at h
            simp only at h
            obtain ⟨e0, he0, ha0⟩ := childWith_sound hS t es ih hrep hn hlen 0 h0
            obtain ⟨e1, he1, ha1⟩ := childWith_sound hS t es ih hrep hn hlen 1 h1
            obtain ⟨e2, he2, ha2⟩ := childWith_sound hS t es ih hrep hn hlen 2 h2
            obtain ⟨e3, he3, ha3⟩ := childWith_sound hS t es ih hrep hn hlen 3 h
            refine ⟨e0 + e1 + e2 + e3, Encloses_node (by omega) he0 he1 he2 he3, ?_⟩
            rw [ha3, ha2, ha1, ha0]
            simp only [enc_add_assoc]

/-! ### The certificate interpreter -/

/-- The enclosure accumulated so far, for `(q, t)`. -/
def accEnc (q : ℚ) (t : ℕ) (acc : Acc) : Enclosure := encOf q t acc.ex acc.entries

theorem runNodeBody_sound {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) {q : ℚ}
    (hq : 1 ≤ q) (t : ℕ) {capTab : Array ℕ} (hcap : CapTabOK S capTab) (cert : ByteArray)
    {r : ℕ} {p : List (Fin 4)} {st : State} (hrep : Rep S p st) (hn : (walkFrom p).Nodup)
    (hlen : p.length + r = steps) {pos : ℕ} {acc acc' : Acc} {pos' : ℕ}
    {onSplit : Acc → Except Err (Acc × ℕ)}
    (hsplit : ∀ acc₁ acc₂ pos₂, onSplit acc₁ = .ok (acc₂, pos₂) →
      ∃ e, Encloses steps S q t p e ∧ accEnc q t acc₂ = accEnc q t acc₁ + e)
    (h : runNodeBody S.toArray capTab cert r st pos acc onSplit = .ok (acc', pos')) :
    ∃ e, Encloses steps S q t p e ∧ accEnc q t acc' = accEnc q t acc + e := by
  unfold runNodeBody at h
  by_cases hpos : pos < cert.size
  · rw [dif_pos hpos] at h
    simp only at h
    by_cases hE : cert[pos] = 69
    · rw [if_pos hE] at h
      by_cases hr : r = 0
      · rw [if_pos hr] at h
        obtain ⟨ex, entries, s1, s2, s3, s4⟩ := acc
        simp only at h
        cases hb : bump ex st.known with
        | none => rw [hb] at h; exact absurd h (by simp)
        | some ex' =>
            rw [hb] at h
            simp only [Except.ok.injEq, Prod.mk.injEq] at h
            obtain ⟨rfl, rfl⟩ := h
            refine ⟨encExact q t st.known, ?_, encOf_bump q t entries hb⟩
            rw [hrep.known]
            exact encExact_sound hS hq t (by omega) hn
      · rw [if_neg hr] at h
        exact absurd h (by simp)
    · rw [if_neg hE] at h
      by_cases hX : cert[pos] = 88
      · rw [if_pos hX] at h
        obtain ⟨ex, entries, s1, s2, s3, s4⟩ := acc
        simp only at h
        cases hb : enumGo S.toArray r st ex with
        | none => rw [hb] at h; exact absurd h (by simp)
        | some ex' =>
            rw [hb] at h
            simp only [Except.ok.injEq, Prod.mk.injEq] at h
            obtain ⟨rfl, rfl⟩ := h
            exact enumGo_sound hS hq t entries r p st ex ex' hrep hn hlen hb
      · rw [if_neg hX] at h
        by_cases hB : cert[pos] = 66
        · rw [if_pos hB] at h
          cases hb : boundEntry S.toArray capTab st with
          | none => rw [hb] at h; exact absurd h (by simp)
          | some e =>
              rw [hb] at h
              obtain ⟨ex, entries, s1, s2, s3, s4⟩ := acc
              simp only [Except.ok.injEq, Prod.mk.injEq] at h
              obtain ⟨rfl, rfl⟩ := h
              exact ⟨encEntry q t e,
                encEntry_sound hS hq t (boundEntry_ok hS hcap hrep hn (by omega) hb),
                encOf_cons q t ex e entries⟩
        · rw [if_neg hB] at h
          by_cases hN : cert[pos] = 78
          · rw [if_pos hN] at h
            exact hsplit _ _ _ h
          · rw [if_neg hN] at h
            exact absurd h (by simp)
  · rw [dif_neg hpos] at h
    exact absurd h (by simp)

theorem runChild_sound {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) {q : ℚ}
    (t : ℕ) (cert : ByteArray) {r : ℕ} {k : State → ℕ → Acc → Except Err (Acc × ℕ)}
    (ih : ∀ (p : List (Fin 4)) (st : State) (pos : ℕ) (acc acc' : Acc) (pos' : ℕ),
      Rep S p st → (walkFrom p).Nodup → p.length + r = steps →
      k st pos acc = .ok (acc', pos') →
      ∃ e, Encloses steps S q t p e ∧ accEnc q t acc' = accEnc q t acc + e)
    {p : List (Fin 4)} {st : State} (hrep : Rep S p st) (hn : (walkFrom p).Nodup)
    (hlen : p.length + (r + 1) = steps) (d : Fin 4) {pos : ℕ} {acc acc' : Acc} {pos' : ℕ}
    (h : runChild S.toArray cert st d pos acc k = .ok (acc', pos')) :
    ∃ e, Encloses steps S q t (p ++ [d]) e ∧ accEnc q t acc' = accEnc q t acc + e := by
  unfold runChild at h
  cases hs : stepState S.toArray st d with
  | none =>
      rw [hs] at h
      simp only at h
      by_cases hpos : pos < cert.size
      · rw [dif_pos hpos] at h
        by_cases hm : cert[pos] = 45
        · rw [if_pos hm] at h
          simp only [Except.ok.injEq, Prod.mk.injEq] at h
          obtain ⟨rfl, rfl⟩ := h
          exact ⟨0, encloses_blocked ((stepState_none hrep d).1 hs), (enc_add_zero _).symm⟩
        · rw [if_neg hm] at h
          exact absurd h (by simp)
      · rw [dif_neg hpos] at h
        exact absurd h (by simp)
  | some st' =>
      rw [hs] at h
      simp only at h
      obtain ⟨hrep', hfree⟩ := stepState_some hrep (by omega) hs
      exact ih (p ++ [d]) st' pos acc acc' pos' hrep' ((nodup_snoc_iff p d).2 ⟨hn, hfree⟩)
        (by simp; omega) h

/-- **Soundness of the certificate interpreter.** Whatever bytes the
certificate contains, if `runTree` accepts the subtree at prefix `p`, the
accumulator grows by a sound enclosure of the subtree masses of `p`, for
every `q ≥ 1` and every threshold `t`. -/
theorem runTree_sound {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) {q : ℚ}
    (hq : 1 ≤ q) (t : ℕ) {capTab : Array ℕ} (hcap : CapTabOK S capTab) (cert : ByteArray) :
    ∀ (r : ℕ) (p : List (Fin 4)) (st : State) (pos : ℕ) (acc acc' : Acc) (pos' : ℕ),
      Rep S p st → (walkFrom p).Nodup → p.length + r = steps →
      runTree S.toArray capTab cert r st pos acc = .ok (acc', pos') →
      ∃ e, Encloses steps S q t p e ∧ accEnc q t acc' = accEnc q t acc + e := by
  intro r
  induction r with
  | zero =>
      intro p st pos acc acc' pos' hrep hn hlen h
      rw [runTree] at h
      exact runNodeBody_sound hS hq t hcap cert hrep hn hlen
        (fun _ _ _ h' => absurd h' (by simp)) h
  | succ r ih =>
      intro p st pos acc acc' pos' hrep hn hlen h
      rw [runTree] at h
      refine runNodeBody_sound hS hq t hcap cert hrep hn hlen ?_ h
      intro acc₁ acc₂ pos₂ h'
      simp only at h'
      cases h0 : runChild S.toArray cert st 0 (pos + 1)
          { acc₁ with splits := acc₁.splits + 1 } (runTree S.toArray capTab cert r) with
      | error e => rw [h0] at h'; exact absurd h' (by simp)
      | ok res0 =>
        obtain ⟨a0, p0⟩ := res0
        rw [h0] at h'
        simp only at h'
        cases h1 : runChild S.toArray cert st 1 p0 a0 (runTree S.toArray capTab cert r) with
        | error e => rw [h1] at h'; exact absurd h' (by simp)
        | ok res1 =>
          obtain ⟨a1, p1⟩ := res1
          rw [h1] at h'
          simp only at h'
          cases h2 : runChild S.toArray cert st 2 p1 a1 (runTree S.toArray capTab cert r) with
          | error e => rw [h2] at h'; exact absurd h' (by simp)
          | ok res2 =>
            obtain ⟨a2, p2⟩ := res2
            rw [h2] at h'
            simp only at h'
            obtain ⟨e0, he0, hc0⟩ := runChild_sound hS t cert ih hrep hn hlen 0 h0
            obtain ⟨e1, he1, hc1⟩ := runChild_sound hS t cert ih hrep hn hlen 1 h1
            obtain ⟨e2, he2, hc2⟩ := runChild_sound hS t cert ih hrep hn hlen 2 h2
            obtain ⟨e3, he3, hc3⟩ := runChild_sound hS t cert ih hrep hn hlen 3 h'
            refine ⟨e0 + e1 + e2 + e3, Encloses_node (by omega) he0 he1 he2 he3, ?_⟩
            rw [hc3, hc2, hc1, hc0]
            simp only [enc_add_assoc]
            rfl

/-- **Soundness of a whole tree.** An accepted tree yields, for every
`q ≥ 1` and every threshold, a sound enclosure of the root masses. -/
theorem runRoot_sound {S : List Bool} {cert : ByteArray} {start : ℕ} {acc : Acc} {stop : ℕ}
    (h : runRoot S.toArray cert start = .ok (acc, stop)) :
    ∃ steps, S.length = steps + 2 ∧
      ∀ q : ℚ, 1 ≤ q → ∀ t : ℕ, Encloses steps S q t [] (accEnc q t acc) := by
  unfold runRoot at h
  by_cases hsz : S.toArray.size < 2
  · rw [if_pos hsz] at h
    exact absurd h (by simp)
  · rw [if_neg hsz] at h
    simp only [List.size_toArray, not_lt] at hsz
    have hS : S.length = (S.length - 2) + 2 := by omega
    refine ⟨S.length - 2, hS, fun q hq t => ?_⟩
    simp only [List.size_toArray] at h
    obtain ⟨e, he, hacc⟩ := runTree_sound hS hq t (capTable_ok S) cert (S.length - 2) []
      (rootState S.toArray) start (emptyAcc S.length) acc stop (rep_root S) (by decide)
      (by simp) h
    have h0 : accEnc q t (emptyAcc S.length) = 0 := encOf_empty q t _
    rw [hacc, h0, enc_zero_add]
    exact he

/-! ### Query evaluation and the string-level theorem -/

theorem evalQuery_fields (acc : Acc) (q : ℚ) (t : ℕ) :
    (evalQuery acc q t).q = q ∧ (evalQuery acc q t).t = t ∧
    (evalQuery acc q t).zLo = (accEnc q t acc).zLo ∧
    (evalQuery acc q t).zHi = (accEnc q t acc).zHi ∧
    (evalQuery acc q t).eLo = (accEnc q t acc).eLo ∧
    (evalQuery acc q t).eHi = (accEnc q t acc).eHi ∧
    (evalQuery acc q t).nLo = (accEnc q t acc).nLo ∧
    (evalQuery acc q t).nHi = (accEnc q t acc).nHi ∧
    (evalQuery acc q t).pLo =
      (accEnc q t acc).eLo / ((accEnc q t acc).eLo + (accEnc q t acc).nHi) ∧
    (evalQuery acc q t).pHi =
      (accEnc q t acc).eHi / ((accEnc q t acc).eHi + (accEnc q t acc).nLo) := by
  refine ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, ?_, ?_⟩
  · show (if (accEnc q t acc).eLo = 0 then 0
      else (accEnc q t acc).eLo / ((accEnc q t acc).eLo + (accEnc q t acc).nHi)) = _
    split_ifs with h0
    · rw [h0, zero_div]
    · rfl
  · show (if (accEnc q t acc).eHi = 0 then 0 else if (accEnc q t acc).nLo = 0 then 1
      else (accEnc q t acc).eHi / ((accEnc q t acc).eHi + (accEnc q t acc).nLo)) = _
    split_ifs with h0 h1
    · rw [h0, zero_div]
    · rw [h1, add_zero, div_self h0]
    · rfl

/-- Soundness of one reported result for the Lean HP model of `seq`: the
echoed `q` is at least one, and the eight reported rationals enclose the
exact partition function, event mass, complement mass and probability
`P(k ≥ t)` for the echoed `q` and threshold. -/
def ResultSound (steps : ℕ) (seq : Fin (steps + 2) → Bool) (r : QueryResult) : Prop :=
  1 ≤ r.q ∧
  Executable.massQ (weightQ steps r.q seq) univ ∈ Set.Icc r.zLo r.zHi ∧
  Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t) ∈
    Set.Icc r.eLo r.eHi ∧
  Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t)ᶜ ∈
    Set.Icc r.nLo r.nHi ∧
  Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t) /
    Executable.massQ (weightQ steps r.q seq) univ ∈ Set.Icc r.pLo r.pHi

theorem evalQuery_sound (steps : ℕ) (seq : Fin (steps + 2) → Bool) (acc : Acc) (q : ℚ)
    (hq : 1 ≤ q) (t : ℕ) (h : Encloses steps (List.ofFn seq) q t [] (accEnc q t acc)) :
    ResultSound steps seq (evalQuery acc q t) := by
  obtain ⟨z1, z2, a1, a2, n1, n2, _, a0, n0⟩ := h
  obtain ⟨hq', ht', hz1, hz2, he1, he2, hn1, hn2, hp1, hp2⟩ := evalQuery_fields acc q t
  unfold ResultSound
  rw [hq', ht', hz1, hz2, he1, he2, hn1, hn2, hp1, hp2, massQ_univ_eq, massQ_event_eq,
    massQ_compl_eq]
  refine ⟨hq, ⟨z1, z2⟩, ⟨a1, a2⟩, ⟨n1, n2⟩, ?_⟩
  rw [massZ_eq_add steps (List.ofFn seq) q t []]
  have hpos := massZ_root_pos steps (List.ofFn seq) (q := q) (by linarith)
  rw [massZ_eq_add steps (List.ofFn seq) q t []] at hpos
  exact ratio_mem_coupled_rat a0 n0 ⟨a1, a2⟩ ⟨n1, n2⟩ hpos

/-- Soundness of a whole report: the echoed sequence is the H/P sequence of
some chain length `steps + 2`, and every reported result is sound for it. -/
def ReportSound (rep : Report) : Prop :=
  ∃ (steps : ℕ) (seq : Fin (steps + 2) → Bool),
    rep.seq = List.ofFn seq ∧ ∀ r ∈ rep.results, ResultSound steps seq r

/-- Inversion of an accepting `checkDoc`: the report is the evaluation of an
accepted tree for the document's sequence, and every query has `q ≥ 1`. -/
theorem checkDoc_ok (doc : Doc) (rep : Report) (h : checkDoc doc = .ok rep) :
    ∃ (acc : Acc) (stop : ℕ),
      runRoot doc.seq.toList.toArray doc.bytes doc.treeStart = .ok (acc, stop) ∧
      rep.seq = doc.seq.toList ∧
      rep.results = doc.queries.map (fun qt => evalQuery acc qt.1 qt.2) ∧
      rep.entries = acc.entries.length ∧
      ∀ qt ∈ doc.queries, (1 : ℚ) ≤ qt.1 := by
  unfold checkDoc at h
  by_cases hq1 : doc.queries.any (fun qt => decide (qt.1 < 1)) = true
  · rw [if_pos hq1] at h
    exact absurd h (by simp)
  · rw [if_neg hq1] at h
    cases hr : runRoot doc.seq doc.bytes doc.treeStart with
    | error e => rw [hr] at h; exact absurd h (by simp)
    | ok res =>
      obtain ⟨acc, stop⟩ := res
      rw [hr] at h
      simp only at h
      by_cases htr : expectLit doc.bytes (lit "\nend\n") stop = some doc.bytes.size
      · rw [if_pos htr] at h
        simp only [Except.ok.injEq] at h
        subst h
        refine ⟨acc, stop, rfl, rfl, rfl, rfl, ?_⟩
        intro qt hqt
        by_contra hlt
        apply hq1
        rw [List.any_eq_true]
        exact ⟨qt, hqt, by simpa using hlt⟩
      · rw [if_neg htr] at h
        exact absurd h (by simp)

/-- The echoed sequence as a function on residue indices. -/
theorem exists_seq_fn (S : List Bool) (steps : ℕ) (hS : S.length = steps + 2) :
    ∃ seq : Fin (steps + 2) → Bool, S = List.ofFn seq := by
  refine ⟨fun i => S[i.1]'(by omega), ?_⟩
  apply List.ext_getElem
  · simp [hS]
  · intro i h1 h2
    rw [List.getElem_ofFn]

theorem checkDoc_sound (doc : Doc) (rep : Report) (h : checkDoc doc = .ok rep) :
    ReportSound rep := by
  obtain ⟨acc, stop, hr, hseq, hres, _, hq⟩ := checkDoc_ok doc rep h
  obtain ⟨steps, hS, hall⟩ := runRoot_sound hr
  obtain ⟨seq, hofFn⟩ := exists_seq_fn doc.seq.toList steps hS
  refine ⟨steps, seq, hseq.trans hofFn, ?_⟩
  intro r hr'
  rw [hres, List.mem_map] at hr'
  obtain ⟨qt, hqt, rfl⟩ := hr'
  apply evalQuery_sound steps seq acc qt.1 (hq qt hqt) qt.2
  rw [← hofFn]
  exact hall qt.1 (hq qt hqt) qt.2

/-- **Soundness of the checker as a function of the raw input bytes.**
No hypothesis about the input: whatever bytes are supplied, every report
returned by `verifyBytes` is sound for the model it echoes. -/
theorem verifyBytes_sound (b : ByteArray) (rep : Report) (h : verifyBytes b = .ok rep) :
    ReportSound rep := by
  unfold verifyBytes at h
  cases hp : parseHeader b with
  | error e => rw [hp] at h; exact absurd h (by simp)
  | ok doc =>
      rw [hp] at h
      exact checkDoc_sound doc rep h

/-- **Soundness of the checker as a function of the raw input string.** -/
theorem verify_sound (input : String) (rep : Report) (h : verify input = .ok rep) :
    ReportSound rep :=
  verifyBytes_sound input.toUTF8 rep h

/-! ### Exactness: no bound entry means point intervals equal to the truth -/

/-- With no bound entries, every reported interval is a single point. -/
theorem evalQuery_point (acc : Acc) (hent : acc.entries = []) (q : ℚ) (t : ℕ) :
    (evalQuery acc q t).zLo = (evalQuery acc q t).zHi ∧
    (evalQuery acc q t).eLo = (evalQuery acc q t).eHi ∧
    (evalQuery acc q t).nLo = (evalQuery acc q t).nHi ∧
    (evalQuery acc q t).pLo = (evalQuery acc q t).pHi := by
  obtain ⟨_, _, hz1, hz2, he1, he2, hn1, hn2, hp1, hp2⟩ := evalQuery_fields acc q t
  have hE : (accEnc q t acc).eLo = (accEnc q t acc).eHi := by simp [accEnc, encOf, hent]
  have hN : (accEnc q t acc).nLo = (accEnc q t acc).nHi := by simp [accEnc, encOf, hent]
  have hZ : (accEnc q t acc).zLo = (accEnc q t acc).zHi := by simp [accEnc, encOf, hent]
  refine ⟨by rw [hz1, hz2, hZ], by rw [he1, he2, hE], by rw [hn1, hn2, hN], ?_⟩
  rw [hp1, hp2, hE, hN]

/-- A reported result that equals the exact model values. -/
def ResultExact (steps : ℕ) (seq : Fin (steps + 2) → Bool) (r : QueryResult) : Prop :=
  1 ≤ r.q ∧
  r.zLo = Executable.massQ (weightQ steps r.q seq) univ ∧
  r.zHi = Executable.massQ (weightQ steps r.q seq) univ ∧
  r.eLo = Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t) ∧
  r.eHi = Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t) ∧
  r.nLo = Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t)ᶜ ∧
  r.nHi = Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t)ᶜ ∧
  r.pLo = Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t) /
    Executable.massQ (weightQ steps r.q seq) univ ∧
  r.pHi = Executable.massQ (weightQ steps r.q seq) (contactEvent steps seq r.t) /
    Executable.massQ (weightQ steps r.q seq) univ

theorem resultExact_of_point {steps : ℕ} {seq : Fin (steps + 2) → Bool} {r : QueryResult}
    (hs : ResultSound steps seq r) (hz : r.zLo = r.zHi) (he : r.eLo = r.eHi)
    (hn : r.nLo = r.nHi) (hp : r.pLo = r.pHi) : ResultExact steps seq r := by
  obtain ⟨hq, ⟨z1, z2⟩, ⟨e1, e2⟩, ⟨n1, n2⟩, ⟨p1, p2⟩⟩ := hs
  refine ⟨hq, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact le_antisymm z1 (hz ▸ z2)
  · exact le_antisymm (hz ▸ z1) z2
  · exact le_antisymm e1 (he ▸ e2)
  · exact le_antisymm (he ▸ e1) e2
  · exact le_antisymm n1 (hn ▸ n2)
  · exact le_antisymm (hn ▸ n1) n2
  · exact le_antisymm p1 (hp ▸ p2)
  · exact le_antisymm (hp ▸ p1) p2

/-- **Exactness.** If an accepted report contains no bound entry (the tree
uses only `N`, `E`, `X` and blocked markers), then every reported endpoint
*equals* the exact value of the Lean HP model: lower and upper endpoints
coincide with the true partition function, event mass, complement mass and
probability.  The fast kernel loses no accuracy relative to exhaustive
enumeration. -/
theorem checkDoc_exact (doc : Doc) (rep : Report) (h : checkDoc doc = .ok rep)
    (h0 : rep.entries = 0) :
    ∃ (steps : ℕ) (seq : Fin (steps + 2) → Bool),
      rep.seq = List.ofFn seq ∧ ∀ r ∈ rep.results, ResultExact steps seq r := by
  obtain ⟨acc, stop, hr, hseq, hres, hent, hq⟩ := checkDoc_ok doc rep h
  obtain ⟨steps, hS, hall⟩ := runRoot_sound hr
  obtain ⟨seq, hofFn⟩ := exists_seq_fn doc.seq.toList steps hS
  have hnil : acc.entries = [] := List.eq_nil_of_length_eq_zero (hent ▸ h0)
  refine ⟨steps, seq, hseq.trans hofFn, ?_⟩
  intro r hr'
  rw [hres, List.mem_map] at hr'
  obtain ⟨qt, hqt, rfl⟩ := hr'
  obtain ⟨hz, he, hn, hp⟩ := evalQuery_point acc hnil qt.1 qt.2
  refine resultExact_of_point ?_ hz he hn hp
  apply evalQuery_sound steps seq acc qt.1 (hq qt hqt) qt.2
  rw [← hofFn]
  exact hall qt.1 (hq qt hqt) qt.2

theorem verifyBytes_exact (b : ByteArray) (rep : Report) (h : verifyBytes b = .ok rep)
    (h0 : rep.entries = 0) :
    ∃ (steps : ℕ) (seq : Fin (steps + 2) → Bool),
      rep.seq = List.ofFn seq ∧ ∀ r ∈ rep.results, ResultExact steps seq r := by
  unfold verifyBytes at h
  cases hp : parseHeader b with
  | error e => rw [hp] at h; exact absurd h (by simp)
  | ok doc =>
      rw [hp] at h
      exact checkDoc_exact doc rep h h0

/-! ### Dominance over the version-1 basic bound leaf -/

/-- **The kernel's bound leaf is never looser than the version-1 basic
bound leaf** (`FrontierCheck.boundEnc false []`, the rule of the Python
`verify_certificate`): at the same prefix, for every `q ≥ 1` and every
threshold, all three upper endpoints are at most the version-1 ones and all
three lower endpoints are at least the version-1 ones. -/
theorem encEntry_le_basic {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2)
    {capTab : Array ℕ} (hcap : CapTabOK S capTab) {p : List (Fin 4)} {st : State}
    (hrep : Rep S p st) {e : Entry} (h : boundEntry S.toArray capTab st = some e)
    {q : ℚ} (hq : 1 ≤ q) (t : ℕ) :
    (encEntry q t e).zHi ≤ (boundEnc false [] steps S q t p).zHi ∧
    (encEntry q t e).eHi ≤ (boundEnc false [] steps S q t p).eHi ∧
    (encEntry q t e).nHi ≤ (boundEnc false [] steps S q t p).nHi ∧
    (boundEnc false [] steps S q t p).zLo ≤ (encEntry q t e).zLo ∧
    (boundEnc false [] steps S q t p).eLo ≤ (encEntry q t e).eLo ∧
    (boundEnc false [] steps S q t p).nLo ≤ (encEntry q t e).nLo := by
  have hq0 : (0 : ℚ) ≤ q := by linarith
  have hpos : ∀ k : ℕ, (0 : ℚ) ≤ q ^ k := fun k => pow_nonneg hq0 k
  unfold boundEntry at h
  cases hfc : capTab[st.placed]? with
  | none => rw [hfc] at h; exact absurd h (by simp)
  | some fc =>
      rw [hfc] at h
      simp only [Option.some.injEq] at h
      subst h
      have hfc' := hcap _ _ hfc
      have hr : S.toArray.size - st.placed = steps - p.length := by
        simp only [List.size_toArray]
        rw [hrep.placed, hS]
        omega
      have hK : min (st.known + fc) (S.toArray.size + 1 - Kernel.minPerim S.toArray.size
          (maxX st.rev - minX st.rev + 1).toNat (maxY st.rev - minY st.rev + 1).toNat) ≤
          kmax S p := by
        refine (min_le_left _ _).trans ?_
        rw [hrep.known, hfc', hrep.placed]
        unfold kmax
        have := futureCap_le_futurePairs S (p.length + 2)
        omega
      have hcnt : ((3 ^ (S.toArray.size - st.placed) : ℕ) : ℚ) = 3 ^ (steps - p.length) := by
        rw [hr]
        push_cast
        rfl
      have hbase : ∀ k : ℕ, k ≤ kmax S p →
          ((3 ^ (S.toArray.size - st.placed) : ℕ) : ℚ) * q ^ k ≤
            3 ^ (steps - p.length) * q ^ kmax S p := by
        intro k hk
        rw [hcnt]
        exact mul_le_mul_of_nonneg_left (qpow_mono hq hk) (by positivity)
      have hbz : (0 : ℚ) ≤ 3 ^ (steps - p.length) * q ^ kmax S p := by positivity
      simp only [encEntry, entZLo, entZHi, entELo, entEHi, entNLo, entNHi, boundEnc,
        envelopeUppers, basicZ, basicA, basicN, kcap, List.map_nil, List.sum_nil, sub_zero,
        Bool.false_eq_true, if_false]
      rw [hrep.known] at hK ⊢
      refine ⟨hbase _ hK, ?_, ?_, ?_, ?_, ?_⟩
      · split_ifs with h1 h2
        · exact hbase _ hK
        · exact absurd (h1.trans hK) h2
        · exact hbz
        · exact le_refl _
      · split_ifs with h1
        · exact hbase _ ((min_le_left _ _).trans hK)
        · exact le_refl _
      · exact mul_nonneg (Nat.cast_nonneg _) (hpos _)
      · split_ifs
        · exact mul_nonneg (Nat.cast_nonneg _) (hpos _)
        · exact le_refl _
      · split_ifs
        · exact mul_nonneg (Nat.cast_nonneg _) (hpos _)
        · exact le_refl _

end Kernel
end FoldBounds
