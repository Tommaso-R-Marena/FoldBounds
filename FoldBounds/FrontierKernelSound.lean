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

theorem scan_none_iff (sx sy : ℤ) (hNew : Bool) (g : ℕ) (R : List (Site × Bool)) (acc : ℕ) :
    scan sx sy hNew g R acc = none ↔ ∃ e ∈ R, e.1 = (sx, sy) := by
  induction R generalizing g acc with
  | nil => simp [scan]
  | cons e rest ih =>
      by_cases he : e.1.1 = sx ∧ e.1.2 = sy
      · have he' : e.1 = (sx, sy) := Prod.ext he.1 he.2
        simp [scan, he']
      · have he' : ¬ e.1 = (sx, sy) := fun h => he ⟨congrArg Prod.fst h, congrArg Prod.snd h⟩
        rw [scan, if_neg he]
        simp only [List.mem_cons, exists_eq_or_imp, he', false_or]
        exact ih _ _

theorem scan_some (sx sy : ℤ) (hNew : Bool) (g : ℕ) (R : List (Site × Bool)) (acc c : ℕ)
    (h : scan sx sy hNew g R acc = some c) : c = acc + scanCount (sx, sy) hNew g R := by
  induction R generalizing g acc with
  | nil => simp [scan] at h; simp [scanCount, h]
  | cons e rest ih =>
      by_cases he : e.1.1 = sx ∧ e.1.2 = sy
      · simp [scan, he] at h
      · rw [scan, if_neg he] at h
        have := ih _ _ h
        rw [this, scanCount]
        have hm : Kernel.manhattan1 e.1 (sx, sy) =
            decide ((e.1.1 - sx).natAbs + (e.1.2 - sy).natAbs = 1) := rfl
        rw [hm]
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
    show (∃ e ∈ st.rev, e.1 = step st.cur d) ↔
      (∃ e ∈ st.rev, e.1 = ((step st.cur d).1, (step st.cur d).2)) from Iff.rfl,
    ← scan_none_iff (step st.cur d).1 (step st.cur d).2 (S.toArray.getD st.placed false) 1
      st.rev 0]
  cases scan (step st.cur d).1 (step st.cur d).2 (S.toArray.getD st.placed false) 1 st.rev 0 <;>
    simp

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
  cases hs : scan (step st.cur d).1 (step st.cur d).2 (S.toArray.getD st.placed false) 1
      st.rev 0 with
  | none => rw [hs] at h; exact absurd h (by simp)
  | some c =>
      rw [hs] at h
      simp only [Option.some.injEq] at h
      subst h
      have hc := scan_some _ _ _ _ _ _ _ hs
      have hlen : (walkFrom p).length = p.length + 2 := walkFrom_length p
      rw [zero_add, hrep.rev, hrep.placed, getD_toArray, hrep.cur, step_eq_advance,
        show ((advance (endFrom first p) d).1, (advance (endFrom first p) d).2) =
          advance (endFrom first p) d from rfl] at hc
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

/-! ### The grid enumerator returns what the list enumerator returns -/

theorem getD_setIfInBounds (g : Array ℕ) (i j v d : ℕ) :
    (g.setIfInBounds i v).getD j d = if i = j ∧ j < g.size then v else g.getD j d := by
  simp only [Array.getD_eq_getD_getElem?, Array.getElem?_setIfInBounds]
  by_cases h : i = j
  · subst h
    by_cases h2 : i < g.size <;> simp [h2]
  · simp [h]

theorem getD_modify (g : Array ℕ) (i j d : ℕ) (f : ℕ → ℕ) :
    (g.modify i f).getD j d = if i = j ∧ j < g.size then f (g.getD j d) else g.getD j d := by
  simp only [Array.getD_eq_getD_getElem?, Array.getElem?_modify]
  by_cases h : i = j
  · subst h
    by_cases h2 : i < g.size <;> simp [h2]
  · simp [h]

theorem getD_default (g : Array ℕ) {j : ℕ} (h : j < g.size) (d d' : ℕ) :
    g.getD j d = g.getD j d' := by
  simp [Array.getD_eq_getD_getElem?, h]

theorem getD_eq_getElem (g : Array ℕ) {j : ℕ} (h : j < g.size) : g.getD j 0 = g[j] := by
  simp [Array.getD_eq_getD_getElem?, h]

/-- Sites within distance `n` of the anchor in each coordinate. -/
def InBox (n : ℕ) (s : Site) : Prop :=
  -(n : ℤ) ≤ s.1 ∧ s.1 ≤ n ∧ -(n : ℤ) ≤ s.2 ∧ s.2 ≤ n

theorem InBox.mono {m n : ℕ} {s : Site} (h : InBox m s) (hmn : m ≤ n) : InBox n s := by
  unfold InBox at *
  omega

theorem InBox.advance {m : ℕ} {s : Site} (h : InBox m s) (d : Fin 4) :
    InBox (m + 1) (advance s d) := by
  obtain ⟨x, y⟩ := s
  unfold InBox at h ⊢
  have hx := h.1; have hx2 := h.2.1; have hy := h.2.2.1; have hy2 := h.2.2.2
  -- `simp [advance]` does not fire under the theorem name `InBox.advance`, so case
  -- on the concrete Fin value and expand the coordinates by hand.
  match d with
  | ⟨0, _⟩ =>
      change (_ : ℤ) ≤ x + 1 ∧ x + 1 ≤ _ ∧ _ ≤ y ∧ y ≤ _
      exact ⟨by omega, by omega, by omega, by omega⟩
  | ⟨1, _⟩ =>
      change (_ : ℤ) ≤ x - 1 ∧ x - 1 ≤ _ ∧ _ ≤ y ∧ y ≤ _
      exact ⟨by omega, by omega, by omega, by omega⟩
  | ⟨2, _⟩ =>
      change (_ : ℤ) ≤ x ∧ x ≤ _ ∧ _ ≤ y + 1 ∧ y + 1 ≤ _
      exact ⟨by omega, by omega, by omega, by omega⟩
  | ⟨3, _⟩ =>
      change (_ : ℤ) ≤ x ∧ x ≤ _ ∧ _ ≤ y - 1 ∧ y - 1 ≤ _
      exact ⟨by omega, by omega, by omega, by omega⟩
  | ⟨n + 4, hn⟩ => omega

/-- Every site of a decoded walk lies in the box of radius `|p| + 1`. -/
theorem walk_box (p : List (Fin 4)) : ∀ a ∈ walkFrom p, InBox (p.length + 1) a := by
  induction p using List.reverseRecOn with
  | nil =>
      intro a ha
      simp only [walkFrom, trace, List.append_nil, List.mem_cons, List.not_mem_nil, or_false] at ha
      rcases ha with rfl | rfl <;> simp [InBox, origin, first]
  | append_singleton p d ih =>
      intro a ha
      rw [walkFrom_append_single, List.mem_append, List.mem_singleton] at ha
      rw [List.length_append, List.length_singleton]
      rcases ha with ha | rfl
      · exact (ih a ha).mono (by omega)
      · exact (ih _ (endFrom_mem_walkFrom p)).advance d

theorem tagFrom_map_fst (S : List Bool) (i : ℕ) (W : List Site) :
    (tagFrom S i W).map Prod.fst = W := by
  induction W generalizing i with
  | nil => rfl
  | cons a W ih => simp [tagFrom, ih]

theorem revOf_map_fst (S : List Bool) (W : List Site) :
    (revOf S W).map Prod.fst = W.reverse := by
  simp [revOf, List.map_reverse, tagFrom_map_fst]

theorem revOf_length (S : List Bool) (W : List Site) : (revOf S W).length = W.length := by
  have := congrArg List.length (revOf_map_fst S W)
  simpa using this

theorem radix_inj {W a b c d : ℕ} (hb : b < W) (hd : d < W) (h : a * W + b = c * W + d) :
    a = c ∧ b = d := by
  have hW : 0 < W := by omega
  have h1 : (a * W + b) / W = a := by
    rw [Nat.add_comm, Nat.add_mul_div_right _ _ hW, Nat.div_eq_of_lt hb, Nat.zero_add]
  have h2 : (c * W + d) / W = c := by
    rw [Nat.add_comm, Nat.add_mul_div_right _ _ hW, Nat.div_eq_of_lt hd, Nat.zero_add]
  have hac : a = c := by rw [← h1, ← h2, h]
  subst hac
  exact ⟨rfl, Nat.add_left_cancel h⟩

theorem codeOf_lt {n : ℕ} {s : Site} (hs : InBox n s) :
    codeOf n s < (2 * n + 1) * (2 * n + 1) := by
  unfold InBox at hs
  unfold codeOf
  have ha : (s.1 + n).toNat ≤ 2 * n := by omega
  have hb : (s.2 + n).toNat < 2 * n + 1 := by omega
  have h1 : (s.1 + n).toNat * (2 * n + 1) ≤ 2 * n * (2 * n + 1) := Nat.mul_le_mul_right _ ha
  have h2 : (2 * n + 1) * (2 * n + 1) = 2 * n * (2 * n + 1) + (2 * n + 1) := by ring
  omega

theorem codeOf_inj {n : ℕ} {a b : Site} (ha : InBox n a) (hb : InBox n b)
    (h : codeOf n a = codeOf n b) : a = b := by
  unfold InBox at ha hb
  unfold codeOf at h
  obtain ⟨h1, h2⟩ := radix_inj (by omega) (by omega) h
  obtain ⟨a1, a2⟩ := a
  obtain ⟨b1, b2⟩ := b
  simp only at h1 h2 ha hb
  have e1 : a1 = b1 := by omega
  have e2 : a2 = b2 := by omega
  rw [e1, e2]

/-- Cell codes of the four neighbours of an interior site. -/
theorem codeOf_advance {n : ℕ} {s : Site} (hx1 : -(n : ℤ) < s.1) (hx2 : s.1 < n)
    (hy1 : -(n : ℤ) < s.2) (hy2 : s.2 < n) :
    codeOf n (advance s 0) = codeOf n s + (2 * n + 1) ∧
    codeOf n (advance s 1) + (2 * n + 1) = codeOf n s ∧
    codeOf n (advance s 2) = codeOf n s + 1 ∧
    codeOf n (advance s 3) + 1 = codeOf n s := by
  obtain ⟨x, y⟩ := s
  simp only at hx1 hx2 hy1 hy2
  refine ⟨?_, ?_, ?_, ?_⟩
  · show codeOf n (x + 1, y) = _
    unfold codeOf
    have e : (x + 1 + (n : ℤ)).toNat = (x + n).toNat + 1 := by omega
    simp only [e]
    ring
  · show codeOf n (x - 1, y) + _ = _
    unfold codeOf
    have e : (x + (n : ℤ)).toNat = (x - 1 + n).toNat + 1 := by omega
    simp only [e]
    ring
  · show codeOf n (x, y + 1) = _
    unfold codeOf
    have e : (y + 1 + (n : ℤ)).toNat = (y + n).toNat + 1 := by omega
    simp only [e]
    ring
  · show codeOf n (x, y - 1) + 1 = _
    unfold codeOf
    have e : (y + (n : ℤ)).toNat = (y - 1 + n).toNat + 1 := by omega
    simp only [e]
    ring

/-- Grid contents determined by the reversed tagged walk (specification). -/
def cellOf (n : ℕ) : List (Site × Bool) → ℕ → ℕ
  | [], _ => 0
  | e :: rest, t => if codeOf n e.1 = t then cellVal rest.length e.2 else cellOf n rest t

theorem cellVal_ne_zero (i : ℕ) (h : Bool) : cellVal i h ≠ 0 := by
  unfold cellVal
  omega

theorem size_loadCells (n : ℕ) (rev : List (Site × Bool)) (g : Array ℕ) :
    (loadCells n rev g).size = g.size := by
  induction rev with
  | nil => rfl
  | cons e rest ih => simp [loadCells, ih]

theorem getD_loadCells (n N : ℕ) :
    ∀ (rev : List (Site × Bool)) (t : ℕ), t < N →
      (loadCells n rev (Array.replicate N 0)).getD t 0 = cellOf n rev t := by
  intro rev
  induction rev with
  | nil =>
      intro t ht
      simp [loadCells, cellOf, Array.getD_eq_getD_getElem?, ht]
  | cons e rest ih =>
      intro t ht
      rw [loadCells, getD_setIfInBounds, size_loadCells, Array.size_replicate, cellOf, ih t ht]
      by_cases hc : codeOf n e.1 = t
      · simp [hc, ht]
      · simp [hc]

theorem cellOf_eq_zero_iff {n : ℕ} {s : Site} (hs : InBox n s) :
    ∀ (rev : List (Site × Bool)), (∀ e ∈ rev, InBox n e.1) →
      (cellOf n rev (codeOf n s) = 0 ↔ ∀ e ∈ rev, e.1 ≠ s) := by
  intro rev
  induction rev with
  | nil => intro _; simp [cellOf]
  | cons e rest ih =>
      intro hbox
      have hrest : ∀ e' ∈ rest, InBox n e'.1 := fun e' he' => hbox e' (List.mem_cons_of_mem _ he')
      rw [cellOf]
      by_cases hc : codeOf n e.1 = codeOf n s
      · have he : e.1 = s := codeOf_inj (hbox e (by simp)) hs hc
        simp [cellVal_ne_zero, he]
      · have he : e.1 ≠ s := fun h => hc (by rw [h])
        simp [hc, he, ih hrest]

theorem cellContact_zero (hNew : Bool) (j : ℕ) : cellContact hNew j 0 = 0 := by
  simp [cellContact]

theorem cellContact_le_one (hNew : Bool) (j v : ℕ) : cellContact hNew j v ≤ 1 := by
  unfold cellContact
  split_ifs <;> omega

theorem cellContact_cellVal (hNew h : Bool) (j i : ℕ) :
    cellContact hNew j (cellVal i h) =
      if (hNew && h && decide (i + 3 ≤ j) && decide ((j - i) % 2 = 1)) = true then 1 else 0 := by
  unfold cellContact cellVal
  have e1 : ((i + 1) * 2 + (if h = true then 1 else 0)) % 2 = if h = true then 1 else 0 := by
    split_ifs <;> omega
  have e2 : ((i + 1) * 2 + (if h = true then 1 else 0)) / 2 = i + 1 := by
    split_ifs <;> omega
  have e3 : j + 1 - (i + 1) = j - i := by omega
  have e4 : decide (i + 1 + 2 ≤ j) = decide (i + 3 ≤ j) := by
    apply Bool.decide_congr
    omega
  rw [e1, e2, e3, e4]
  cases hNew <;> cases h <;> simp

/-- Exactly one direction leads from `s` to a lattice neighbour `a`. -/
theorem sum_dir_eq (s a : Site) (C : ℕ) :
    (∑ d : Fin 4, if a = advance s d then C else 0) =
      if LatticeModel.manhattan1 a s = true then C else 0 := by
  by_cases hm : LatticeModel.manhattan1 a s = true
  · have hadj : adjacent s a = true := by
      rw [adjacent_eq_manhattan1]
      exact Geometry.manhattan1_symm hm
    have ha : advance s (stepDir s a) = a := advance_stepDir hadj
    rw [if_pos hm, Finset.sum_eq_single (stepDir s a)]
    · rw [if_pos ha.symm]
    · intro d _ hd
      rw [if_neg]
      intro h
      exact hd (advance_injective s (h.symm.trans ha.symm))
    · intro h
      exact absurd (Finset.mem_univ _) h
  · rw [if_neg hm]
    refine Finset.sum_eq_zero (fun d _ => ?_)
    rw [if_neg]
    intro h
    apply hm
    rw [h]
    have := advance_adjacent s d
    rw [adjacent_eq_manhattan1] at this
    exact Geometry.manhattan1_symm this

/-- **The four neighbouring cells see exactly the contacts found by the
list scan.** -/
theorem grid_contacts {n : ℕ} (hNew : Bool) (j : ℕ) {s : Site}
    (hs : ∀ d : Fin 4, InBox n (advance s d)) :
    ∀ (rev : List (Site × Bool)), (∀ e ∈ rev, InBox n e.1) → (rev.map Prod.fst).Nodup →
      rev.length ≤ j →
      (∑ d : Fin 4, cellContact hNew j (cellOf n rev (codeOf n (advance s d)))) =
        scanCount s hNew (j + 1 - rev.length) rev := by
  intro rev
  induction rev with
  | nil => intro _ _ _; simp [cellOf, cellContact_zero, scanCount]
  | cons e rest ih =>
      intro hbox hnd hlen
      obtain ⟨a, h⟩ := e
      have hrest : ∀ e' ∈ rest, InBox n e'.1 := fun e' he' => hbox e' (List.mem_cons_of_mem _ he')
      simp only [List.map_cons, List.nodup_cons, List.mem_map, not_exists, not_and] at hnd
      obtain ⟨hnot, hnd'⟩ := hnd
      simp only [List.length_cons] at hlen
      have hterm : ∀ d : Fin 4,
          cellContact hNew j (cellOf n ((a, h) :: rest) (codeOf n (advance s d))) =
            (if a = advance s d then cellContact hNew j (cellVal rest.length h) else 0) +
              cellContact hNew j (cellOf n rest (codeOf n (advance s d))) := by
        intro d
        by_cases had : a = advance s d
        · have hz : cellOf n rest (codeOf n (advance s d)) = 0 :=
            (cellOf_eq_zero_iff (hs d) rest hrest).2 (fun e' he' h' => hnot e' he' (h'.trans had.symm))
          rw [if_pos had, hz, cellContact_zero, add_zero, cellOf, if_pos (by rw [had])]
        · have hc : codeOf n a ≠ codeOf n (advance s d) :=
            fun hc => had (codeOf_inj (hbox (a, h) (by simp)) (hs d) hc)
          rw [if_neg had, zero_add, cellOf, if_neg hc]
      rw [Finset.sum_congr rfl (fun d _ => hterm d), Finset.sum_add_distrib, sum_dir_eq,
        ih hrest hnd' (by omega), scanCount, List.length_cons, cellContact_cellVal]
      have hg : j + 1 - (rest.length + 1) + 1 = j + 1 - rest.length := by omega
      rw [hg]
      congr 1
      have e1 : j + 1 - (rest.length + 1) = j - rest.length := by omega
      have e2 : decide (3 ≤ j - rest.length) = decide (rest.length + 3 ≤ j) := by
        apply Bool.decide_congr
        omega
      rw [e1, e2, manhattan1_eq]
      cases hNew <;> cases h <;> cases decide (rest.length + 3 ≤ j) <;>
        cases decide ((j - rest.length) % 2 = 1) <;> cases LatticeModel.manhattan1 a s <;> simp

/-- The grid part and the size of an array are unchanged. -/
def Frame (base : ℕ) (g g' : Array ℕ) : Prop :=
  g'.size = g.size ∧ ∀ u, u < base → g'.getD u 0 = g.getD u 0

theorem Frame.refl (base : ℕ) (g : Array ℕ) : Frame base g g := ⟨rfl, fun _ _ => rfl⟩

theorem Frame.trans {base : ℕ} {g g' g'' : Array ℕ} (h1 : Frame base g g')
    (h2 : Frame base g' g'') : Frame base g g'' :=
  ⟨h2.1.trans h1.1, fun u hu => (h2.2 u hu).trans (h1.2 u hu)⟩

theorem gridChild_frame (seq : Array Bool) (W base t j known : ℕ) (g : Array ℕ)
    (k : ℕ → ℕ → ℕ → Array ℕ → Array ℕ)
    (hk : ∀ t' j' kn' g', Frame base g' (k t' j' kn' g')) :
    Frame base g (gridChild seq W t j known g k) := by
  unfold gridChild
  by_cases h0 : g.getD t 1 = 0
  · rw [if_pos h0]
    simp only
    refine ⟨by rw [Array.size_setIfInBounds, (hk _ _ _ _).1, Array.size_setIfInBounds],
      fun u hu => ?_⟩
    rw [getD_setIfInBounds, (hk _ _ _ _).1, Array.size_setIfInBounds]
    by_cases hut : t = u
    · subst hut
      by_cases hts : t < g.size
      · have hz : g.getD t 0 = 0 := by rw [getD_default g hts 0 1]; exact h0
        rw [if_pos ⟨rfl, hts⟩, hz]
      · exfalso
        simp [Array.getD_eq_getD_getElem?, hts] at h0
    · rw [if_neg (fun h => hut h.1), (hk _ _ _ _).2 u hu, getD_setIfInBounds,
        if_neg (fun h => hut h.1)]
  · rw [if_neg h0]
    exact Frame.refl base g

theorem enumGrid_frame (seq : Array Bool) (W base : ℕ) :
    ∀ (r code j known : ℕ) (g : Array ℕ),
      Frame base g (enumGrid seq W base r code j known g) := by
  intro r
  induction r with
  | zero =>
      intro code j known g
      rw [enumGrid]
      refine ⟨Array.size_modify, fun u hu => ?_⟩
      rw [getD_modify, if_neg]
      intro h
      omega
  | succ r ih =>
      intro code j known g
      rw [enumGrid]
      exact (((gridChild_frame seq W base _ j known g _ ih).trans
        (gridChild_frame seq W base _ j known _ _ ih)).trans
        (gridChild_frame seq W base _ j known _ _ ih)).trans
        (gridChild_frame seq W base _ j known _ _ ih)

/-- The array's grid part holds the cells of the placed walk. -/
def GridRep (n : ℕ) (g : Array ℕ) (rev : List (Site × Bool)) : Prop :=
  ∀ u, u < (2 * n + 1) * (2 * n + 1) → g.getD u 0 = cellOf n rev u

/-- The array's histogram part holds the histogram `ex`. -/
def HistRep (base : ℕ) (g ex : Array ℕ) : Prop :=
  g.size = base + ex.size ∧ ∀ k, k < ex.size → g.getD (base + k) 0 = ex.getD k 0

theorem gridChild_refines {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) {r : ℕ}
    (ih : ∀ (p : List (Fin 4)) (st : State) (g ex : Array ℕ), Rep S p st →
      (walkFrom p).Nodup → p.length + r = steps → GridRep S.length g st.rev →
      HistRep ((2 * S.length + 1) * (2 * S.length + 1)) g ex → st.known + 4 * r < ex.size →
      ∃ ex', enumGo S.toArray r st ex = some ex' ∧
        HistRep ((2 * S.length + 1) * (2 * S.length + 1))
          (enumGrid S.toArray (2 * S.length + 1) ((2 * S.length + 1) * (2 * S.length + 1)) r
            (codeOf S.length st.cur) st.placed st.known g) ex')
    {p : List (Fin 4)} {st : State} (hrep : Rep S p st) (hn : (walkFrom p).Nodup)
    (hlen : p.length + (r + 1) = steps) {g ex : Array ℕ} (hg : GridRep S.length g st.rev)
    (hh : HistRep ((2 * S.length + 1) * (2 * S.length + 1)) g ex)
    (hk : st.known + 4 * (r + 1) < ex.size) (d : Fin 4) :
    ∃ ex', childWith S.toArray st d ex (enumGo S.toArray r) = some ex' ∧
      ex'.size = ex.size ∧
      HistRep ((2 * S.length + 1) * (2 * S.length + 1))
        (gridChild S.toArray (2 * S.length + 1) (codeOf S.length (advance st.cur d)) st.placed
          st.known g (enumGrid S.toArray (2 * S.length + 1)
            ((2 * S.length + 1) * (2 * S.length + 1)) r)) ex' ∧
      GridRep S.length
        (gridChild S.toArray (2 * S.length + 1) (codeOf S.length (advance st.cur d)) st.placed
          st.known g (enumGrid S.toArray (2 * S.length + 1)
            ((2 * S.length + 1) * (2 * S.length + 1)) r)) st.rev := by
  have hwlen : (walkFrom p).length = p.length + 2 := walkFrom_length p
  have hrevlen : st.rev.length = st.placed := by
    rw [hrep.rev, revOf_length, hwlen, hrep.placed]
  have hcurbox : InBox (p.length + 1) st.cur := by
    rw [hrep.cur]
    exact walk_box p _ (endFrom_mem_walkFrom p)
  have hsbox : InBox (p.length + 2) (advance st.cur d) := hcurbox.advance d
  have hnbr : ∀ d' : Fin 4, InBox S.length (advance (advance st.cur d) d') :=
    fun d' => (hsbox.advance d').mono (by omega)
  have hsbox' : InBox S.length (advance st.cur d) := hsbox.mono (by omega)
  have hrevbox : ∀ e ∈ st.rev, InBox S.length e.1 := by
    intro e he
    rw [hrep.rev] at he
    exact (walk_box p _ ((mem_revOf S (walkFrom p) e.1).1 ⟨e, he, rfl⟩)).mono (by omega)
  have hrevnd : (st.rev.map Prod.fst).Nodup := by
    rw [hrep.rev, revOf_map_fst]
    exact List.nodup_reverse.2 hn
  have hcode_lt := codeOf_lt hsbox'
  have hsz : (2 * S.length + 1) * (2 * S.length + 1) ≤ g.size := by rw [hh.1]; omega
  unfold childWith gridChild
  cases hs : stepState S.toArray st d with
  | none =>
      have hocc : advance (endFrom first p) d ∈ walkFrom p := (stepState_none hrep d).1 hs
      have hne : cellOf S.length st.rev (codeOf S.length (advance st.cur d)) ≠ 0 := by
        intro hz
        have := (cellOf_eq_zero_iff hsbox' st.rev hrevbox).1 hz
        rw [hrep.rev, hrep.cur] at this
        obtain ⟨e, he, hes⟩ := (mem_revOf S (walkFrom p) _).2 hocc
        exact this e he hes
      have hg1 : ¬ g.getD (codeOf S.length (advance st.cur d)) 1 = 0 := by
        rw [getD_default g (by omega) 1 0, hg _ hcode_lt]
        exact hne
      rw [if_neg hg1]
      exact ⟨ex, rfl, rfl, hh, hg⟩
  | some st' =>
      obtain ⟨hrep', hfree⟩ := stepState_some hrep (by omega) hs
      have hz : cellOf S.length st.rev (codeOf S.length (advance st.cur d)) = 0 := by
        refine (cellOf_eq_zero_iff hsbox' st.rev hrevbox).2 (fun e he hes => hfree ?_)
        rw [hrep.rev] at he
        rw [hrep.cur] at hes
        exact (mem_revOf S (walkFrom p) _).1 ⟨e, he, hes⟩
      have hg0 : g.getD (codeOf S.length (advance st.cur d)) 1 = 0 := by
        rw [getD_default g (by omega) 1 0, hg _ hcode_lt]
        exact hz
      rw [if_pos hg0]
      simp only
      -- the four neighbouring cells
      have hint : -(S.length : ℤ) < (advance st.cur d).1 ∧ (advance st.cur d).1 < S.length ∧
          -(S.length : ℤ) < (advance st.cur d).2 ∧ (advance st.cur d).2 < S.length := by
        unfold InBox at hsbox
        omega
      obtain ⟨c0, c1, c2, c3⟩ := codeOf_advance hint.1 hint.2.1 hint.2.2.1 hint.2.2.2
      have hc1 : codeOf S.length (advance st.cur d) - (2 * S.length + 1) =
          codeOf S.length (advance (advance st.cur d) 1) := by omega
      have hc3 : codeOf S.length (advance st.cur d) - 1 =
          codeOf S.length (advance (advance st.cur d) 3) := by omega
      have hcell : ∀ d' : Fin 4, g.getD (codeOf S.length (advance (advance st.cur d) d')) 0 =
          cellOf S.length st.rev (codeOf S.length (advance (advance st.cur d) d')) :=
        fun d' => hg _ (codeOf_lt (hnbr d'))
      have hsum : cellContact (S.toArray.getD st.placed false) st.placed
            (g.getD (codeOf S.length (advance st.cur d) + (2 * S.length + 1)) 0) +
          cellContact (S.toArray.getD st.placed false) st.placed
            (g.getD (codeOf S.length (advance st.cur d) - (2 * S.length + 1)) 0) +
          cellContact (S.toArray.getD st.placed false) st.placed
            (g.getD (codeOf S.length (advance st.cur d) + 1) 0) +
          cellContact (S.toArray.getD st.placed false) st.placed
            (g.getD (codeOf S.length (advance st.cur d) - 1) 0) =
          scanCount (advance st.cur d) (S.toArray.getD st.placed false) 1 st.rev := by
        rw [← c0, hc1, ← c2, hc3, hcell 0, hcell 1, hcell 2, hcell 3]
        have := grid_contacts (S.toArray.getD st.placed false) st.placed hnbr st.rev hrevbox
          hrevnd (by omega)
        rw [Fin.sum_univ_four, hrevlen] at this
        rw [this]
        congr 1
        omega
      -- the list step computes the same state
      -- Parenthesize the structure literal: bare `{ field := ... }` inside
      -- `have h : _ = {..} :=` is parsed as ending the `have` at the first `:=`.
      have hst' : st' = ({
          cur := advance st.cur d
          rev := (advance st.cur d, S.toArray.getD st.placed false) :: st.rev
          placed := st.placed + 1
          known := st.known +
            scanCount (advance st.cur d) (S.toArray.getD st.placed false) 1 st.rev }) := by
        unfold stepState at hs
        simp only at hs
        cases hsc : scan (step st.cur d).1 (step st.cur d).2 (S.toArray.getD st.placed false) 1
            st.rev 0 with
        | none => rw [hsc] at hs; exact absurd hs (by simp)
        | some c =>
            rw [hsc] at hs
            simp only [Option.some.injEq] at hs
            have hc := scan_some _ _ _ _ _ _ _ hsc
            rw [zero_add] at hc
            rw [← hs, hc]
            rfl
      rw [hsum]
      have hle4 : scanCount (advance st.cur d) (S.toArray.getD st.placed false) 1 st.rev ≤ 4 := by
        rw [← hsum]
        have h1 := cellContact_le_one (S.toArray.getD st.placed false) st.placed
          (g.getD (codeOf S.length (advance st.cur d) + (2 * S.length + 1)) 0)
        have h2 := cellContact_le_one (S.toArray.getD st.placed false) st.placed
          (g.getD (codeOf S.length (advance st.cur d) - (2 * S.length + 1)) 0)
        have h3 := cellContact_le_one (S.toArray.getD st.placed false) st.placed
          (g.getD (codeOf S.length (advance st.cur d) + 1) 0)
        have h4 := cellContact_le_one (S.toArray.getD st.placed false) st.placed
          (g.getD (codeOf S.length (advance st.cur d) - 1) 0)
        omega
      -- the array after writing the new cell represents the extended walk
      have hg1 : GridRep S.length
          (g.setIfInBounds (codeOf S.length (advance st.cur d))
            (cellVal st.placed (S.toArray.getD st.placed false))) st'.rev := by
        intro u hu
        rw [getD_setIfInBounds, hst', cellOf, hrevlen]
        by_cases hcu : codeOf S.length (advance st.cur d) = u
        · rw [if_pos ⟨hcu, by omega⟩, if_pos hcu]
        · rw [if_neg (fun h => hcu h.1), if_neg hcu]
          exact hg u hu
      have hh1 : HistRep ((2 * S.length + 1) * (2 * S.length + 1))
          (g.setIfInBounds (codeOf S.length (advance st.cur d))
            (cellVal st.placed (S.toArray.getD st.placed false))) ex := by
        refine ⟨by rw [Array.size_setIfInBounds]; exact hh.1, fun k hk' => ?_⟩
        rw [getD_setIfInBounds, if_neg (fun h => by omega)]
        exact hh.2 k hk'
      have hkn : st'.known + 4 * r < ex.size := by
        rw [hst']
        simp only
        omega
      obtain ⟨ex', he, hhe⟩ := ih (p ++ [d]) st' _ ex hrep' ((nodup_snoc_iff p d).2 ⟨hn, hfree⟩)
        (by simp; omega) hg1 hh1 hkn
      have hfr := enumGrid_frame S.toArray (2 * S.length + 1)
        ((2 * S.length + 1) * (2 * S.length + 1)) r (codeOf S.length st'.cur) st'.placed st'.known
        (g.setIfInBounds (codeOf S.length (advance st.cur d))
          (cellVal st.placed (S.toArray.getD st.placed false)))
      have hcur' : st'.cur = advance st.cur d := by rw [hst']
      have hplaced' : st'.placed = st.placed + 1 := by rw [hst']
      have hknown' : st'.known = st.known +
          scanCount (advance st.cur d) (S.toArray.getD st.placed false) 1 st.rev := by rw [hst']
      rw [hcur', hplaced', hknown'] at hhe hfr
      have hexsz : ex'.size = ex.size := by
        have h1 := hhe.1
        have h2 := hfr.1
        rw [Array.size_setIfInBounds] at h2
        have h3 := hh.1
        omega
      refine ⟨ex', he, hexsz, ⟨?_, fun k hk' => ?_⟩, fun u hu => ?_⟩
      · rw [Array.size_setIfInBounds]
        exact hhe.1
      · rw [getD_setIfInBounds, if_neg (fun h => by omega)]
        exact hhe.2 k hk'
      · rw [getD_setIfInBounds]
        by_cases hcu : codeOf S.length (advance st.cur d) = u
        · rw [if_pos ⟨hcu, by rw [hfr.1, Array.size_setIfInBounds]; omega⟩, ← hcu, hz]
        · rw [if_neg (fun h => hcu h.1), hfr.2 u hu, getD_setIfInBounds,
            if_neg (fun h => hcu h.1)]
          exact hg u hu

/-- **Refinement.** Started on an array that represents the walk and the
histogram, the grid enumerator leaves in the histogram part exactly the
histogram returned by the list enumerator. -/
theorem enumGrid_refines {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2) :
    ∀ (r : ℕ) (p : List (Fin 4)) (st : State) (g ex : Array ℕ), Rep S p st →
      (walkFrom p).Nodup → p.length + r = steps → GridRep S.length g st.rev →
      HistRep ((2 * S.length + 1) * (2 * S.length + 1)) g ex → st.known + 4 * r < ex.size →
      ∃ ex', enumGo S.toArray r st ex = some ex' ∧
        HistRep ((2 * S.length + 1) * (2 * S.length + 1))
          (enumGrid S.toArray (2 * S.length + 1) ((2 * S.length + 1) * (2 * S.length + 1)) r
            (codeOf S.length st.cur) st.placed st.known g) ex' := by
  intro r
  induction r with
  | zero =>
      intro p st g ex hrep hn hlen hg hh hk
      have hkn : st.known < ex.size := by omega
      rw [enumGo, enumGrid]
      refine ⟨ex.set st.known (ex[st.known] + 1), by simp [bump, hkn], ?_, fun k hk' => ?_⟩
      · rw [Array.size_modify, Array.size_set]
        exact hh.1
      · rw [Array.size_set] at hk'
        rw [getD_modify]
        by_cases hkk : st.known = k
        · subst hkk
          rw [if_pos ⟨rfl, by rw [hh.1]; omega⟩, hh.2 _ hk', getD_eq_getElem ex hkn]
          simp [Array.getD_eq_getD_getElem?, hkn]
        · rw [if_neg (fun h => hkk (by omega)), hh.2 k hk']
          simp [Array.getD_eq_getD_getElem?, hkk]
  | succ r ih =>
      intro p st g ex hrep hn hlen hg hh hk
      have hcurbox : InBox (p.length + 1) st.cur := by
        rw [hrep.cur]
        exact walk_box p _ (endFrom_mem_walkFrom p)
      have hint : -(S.length : ℤ) < st.cur.1 ∧ st.cur.1 < S.length ∧
          -(S.length : ℤ) < st.cur.2 ∧ st.cur.2 < S.length := by
        unfold InBox at hcurbox
        omega
      obtain ⟨c0, c1, c2, c3⟩ := codeOf_advance hint.1 hint.2.1 hint.2.2.1 hint.2.2.2
      have hc1 : codeOf S.length st.cur - (2 * S.length + 1) =
          codeOf S.length (advance st.cur 1) := by omega
      have hc3 : codeOf S.length st.cur - 1 = codeOf S.length (advance st.cur 3) := by omega
      -- Unfold the succ clauses by definitional equality, then rewrite neighbour
      -- codes into `advance` form so they match `gridChild_refines`.
      rw [show enumGo S.toArray (r + 1) st ex =
          (match childWith S.toArray st 0 ex (enumGo S.toArray r) with
            | none => none
            | some ex =>
              match childWith S.toArray st 1 ex (enumGo S.toArray r) with
              | none => none
              | some ex =>
                match childWith S.toArray st 2 ex (enumGo S.toArray r) with
                | none => none
                | some ex => childWith S.toArray st 3 ex (enumGo S.toArray r)) from rfl]
      rw [show enumGrid S.toArray (2 * S.length + 1)
            ((2 * S.length + 1) * (2 * S.length + 1)) (r + 1)
            (codeOf S.length st.cur) st.placed st.known g =
          (let g := gridChild S.toArray (2 * S.length + 1)
              (codeOf S.length st.cur + (2 * S.length + 1)) st.placed st.known g
              (enumGrid S.toArray (2 * S.length + 1)
                ((2 * S.length + 1) * (2 * S.length + 1)) r)
            let g := gridChild S.toArray (2 * S.length + 1)
              (codeOf S.length st.cur - (2 * S.length + 1)) st.placed st.known g
              (enumGrid S.toArray (2 * S.length + 1)
                ((2 * S.length + 1) * (2 * S.length + 1)) r)
            let g := gridChild S.toArray (2 * S.length + 1)
              (codeOf S.length st.cur + 1) st.placed st.known g
              (enumGrid S.toArray (2 * S.length + 1)
                ((2 * S.length + 1) * (2 * S.length + 1)) r)
            gridChild S.toArray (2 * S.length + 1)
              (codeOf S.length st.cur - 1) st.placed st.known g
              (enumGrid S.toArray (2 * S.length + 1)
                ((2 * S.length + 1) * (2 * S.length + 1)) r)) from rfl]
      rw [← c0, hc1, ← c2, hc3]
      obtain ⟨e0, h0, s0, hh0, hg0⟩ := gridChild_refines hS ih hrep hn hlen hg hh hk 0
      obtain ⟨e1, h1, s1, hh1, hg1⟩ := gridChild_refines hS ih hrep hn hlen hg0 hh0
        (by rw [s0]; exact hk) 1
      obtain ⟨e2, h2, s2, hh2, hg2⟩ := gridChild_refines hS ih hrep hn hlen hg1 hh1
        (by rw [s1, s0]; exact hk) 2
      obtain ⟨e3, h3, _s3, hh3, _hg3⟩ := gridChild_refines hS ih hrep hn hlen hg2 hh2
        (by rw [s2, s1, s0]; exact hk) 3
      -- Reduce the nested `match some _` after each successful child.
      simp only [h0, h1, h2, h3]
      exact ⟨e3, rfl, hh3⟩

/-- **The fast enumerator is the list enumerator.** For every state that
represents a self-avoiding prefix, `enumFast` returns exactly the result of
`enumGo`. -/
theorem enumFast_eq {steps : ℕ} {S : List Bool} (hS : S.length = steps + 2)
    {p : List (Fin 4)} {st : State} (hrep : Rep S p st) (hn : (walkFrom p).Nodup) {r : ℕ}
    (hlen : p.length + r = steps) (ex : Array ℕ) :
    enumFast S.toArray r st ex = enumGo S.toArray r st ex := by
  unfold enumFast
  simp only [List.size_toArray]
  by_cases hk : st.known + 4 * r < ex.size
  · rw [if_pos hk]
    have hgrid : GridRep S.length
        (loadCells S.length st.rev
          (Array.replicate ((2 * S.length + 1) * (2 * S.length + 1)) 0) ++ ex) st.rev := by
      intro u hu
      have hlt : u < (loadCells S.length st.rev
          (Array.replicate ((2 * S.length + 1) * (2 * S.length + 1)) 0)).size := by
        rw [size_loadCells, Array.size_replicate]
        exact hu
      rw [← getD_loadCells S.length _ st.rev u hu]
      simp [Array.getD_eq_getD_getElem?, Array.getElem?_append_left hlt]
    have hhist : HistRep ((2 * S.length + 1) * (2 * S.length + 1))
        (loadCells S.length st.rev
          (Array.replicate ((2 * S.length + 1) * (2 * S.length + 1)) 0) ++ ex) ex := by
      refine ⟨by rw [Array.size_append, size_loadCells, Array.size_replicate], fun k hk' => ?_⟩
      set baseArr := loadCells S.length st.rev
          (Array.replicate ((2 * S.length + 1) * (2 * S.length + 1)) 0) with hbaseArr
      have hsz : baseArr.size = (2 * S.length + 1) * (2 * S.length + 1) := by
        rw [hbaseArr, size_loadCells, Array.size_replicate]
      -- Rewrite only the index `(2n+1)^2` to `baseArr.size`, not the array itself.
      conv => lhs; arg 2; rw [← hsz]
      have hle : baseArr.size ≤ baseArr.size + k := Nat.le_add_right _ _
      simp only [Array.getD_eq_getD_getElem?, Array.getElem?_append_right hle]
      rw [Nat.add_sub_cancel_left]
    obtain ⟨ex', he, hsz, hh⟩ := enumGrid_refines hS r p st _ ex hrep hn hlen hgrid hhist hk
    have hfr := enumGrid_frame S.toArray (2 * S.length + 1)
      ((2 * S.length + 1) * (2 * S.length + 1)) r (codeOf S.length st.cur) st.placed st.known
      (loadCells S.length st.rev
        (Array.replicate ((2 * S.length + 1) * (2 * S.length + 1)) 0) ++ ex)
    have hexsz : ex'.size = ex.size := by
      have h2 := hfr.1
      rw [Array.size_append, size_loadCells, Array.size_replicate] at h2
      omega
    rw [he]
    congr 1
    apply Array.ext
    · rw [Array.size_extract, hsz, hexsz]
      omega
    · intro i h1 h2
      rw [Array.getElem_extract, ← getD_eq_getElem ex' h2, ← hh i h2]
      rw [getD_eq_getElem]
  · rw [if_neg hk]

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
        cases hb : enumFast S.toArray r st ex with
        | none => rw [hb] at h; exact absurd h (by simp)
        | some ex' =>
            rw [hb] at h
            simp only [Except.ok.injEq, Prod.mk.injEq] at h
            obtain ⟨rfl, rfl⟩ := h
            rw [enumFast_eq hS hrep hn hlen] at hb
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

/-! ### Boolean entry points: one decidable check yields a theorem about the model -/

/-- A sound report, specialised to a known sequence function. -/
theorem ReportSound.of_seq {rep : Report} (h : ReportSound rep) {steps : ℕ}
    (seq : Fin (steps + 2) → Bool) (hseq : rep.seq = List.ofFn seq) :
    ∀ r ∈ rep.results, ResultSound steps seq r := by
  obtain ⟨steps', seq', hs', hall⟩ := h
  have he : List.ofFn seq' = List.ofFn seq := hs'.symm.trans hseq
  have hlen : steps' = steps := by
    have := congrArg List.length he
    simp only [List.length_ofFn] at this
    omega
  subst hlen
  have : seq' = seq := List.ofFn_injective he
  subst this
  exact hall

/-- Boolean test: the bytes are accepted for the sequence `S` and some
reported result for `(q, t)` has its probability interval inside `[lo, hi]`. -/
def certifiesInterval (b : ByteArray) (S : List Bool) (q : ℚ) (t : ℕ) (lo hi : ℚ) : Bool :=
  match verifyBytes b with
  | .ok rep =>
      rep.seq == S &&
        rep.results.any (fun r => r.q == q && r.t == t && decide (lo ≤ r.pLo) && decide (r.pHi ≤ hi))
  | .error _ => false

/-- **One Boolean check gives a probability enclosure for the Lean model.**
If `certifiesInterval` evaluates to `true` (for instance by `decide +kernel`
on a small certificate, or by running the compiled checker), then the exact
equilibrium probability `P(k ≥ t)` of the HP model of `seq` at Boltzmann
factor `q` lies in `[lo, hi]`. -/
theorem certifiesInterval_sound {steps : ℕ} (seq : Fin (steps + 2) → Bool) (b : ByteArray)
    (q : ℚ) (t : ℕ) (lo hi : ℚ) (h : certifiesInterval b (List.ofFn seq) q t lo hi = true) :
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq t) /
      Executable.massQ (weightQ steps q seq) univ ∈ Set.Icc lo hi := by
  unfold certifiesInterval at h
  cases hv : verifyBytes b with
  | error e => rw [hv] at h; exact absurd h (by simp)
  | ok rep =>
      rw [hv] at h
      simp only [Bool.and_eq_true, beq_iff_eq, List.any_eq_true, decide_eq_true_eq] at h
      obtain ⟨hseq, r, hr, ⟨⟨⟨hq, ht⟩, hlo⟩, hhi⟩⟩ := h
      obtain ⟨_, _, _, _, hp⟩ := (verifyBytes_sound b rep hv).of_seq seq hseq r hr
      rw [hq, ht] at hp
      exact ⟨hlo.trans hp.1, hp.2.trans hhi⟩

/-- Boolean test: the bytes are accepted for the sequence `S` without any
bound entry and some reported result for `(q, t)` has probability `p`. -/
def certifiesExact (b : ByteArray) (S : List Bool) (q : ℚ) (t : ℕ) (p : ℚ) : Bool :=
  match verifyBytes b with
  | .ok rep =>
      rep.entries == 0 && rep.seq == S &&
        rep.results.any (fun r => r.q == q && r.t == t && r.pLo == p)
  | .error _ => false

/-- **One Boolean check gives the exact probability of the Lean model.** -/
theorem certifiesExact_sound {steps : ℕ} (seq : Fin (steps + 2) → Bool) (b : ByteArray)
    (q : ℚ) (t : ℕ) (p : ℚ) (h : certifiesExact b (List.ofFn seq) q t p = true) :
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq t) /
      Executable.massQ (weightQ steps q seq) univ = p := by
  unfold certifiesExact at h
  cases hv : verifyBytes b with
  | error e => rw [hv] at h; exact absurd h (by simp)
  | ok rep =>
      rw [hv] at h
      simp only [Bool.and_eq_true, beq_iff_eq, List.any_eq_true] at h
      obtain ⟨⟨h0, hseq⟩, r, hr, ⟨⟨hq, ht⟩, hp⟩⟩ := h
      obtain ⟨steps', seq', hs', hall⟩ := verifyBytes_exact b rep hv h0
      have he : List.ofFn seq' = List.ofFn seq := hs'.symm.trans hseq
      have hlen : steps' = steps := by
        have := congrArg List.length he
        simp only [List.length_ofFn] at this
        omega
      subst hlen
      have : seq' = seq := List.ofFn_injective he
      subst this
      obtain ⟨_, _, _, _, _, _, _, hpl, _⟩ := hall r hr
      rw [hq, ht] at hpl
      rw [← hpl, hp]

end Kernel
end FoldBounds
