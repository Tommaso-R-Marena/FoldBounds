module

public import FoldBounds.FrontierNext

/-!
# Monotone refinement of the finite-memory relaxation

Enlarging the memory `k` of the relaxation can only shrink the set of
relaxed walks and tighten every per-state contact window.  We prove this
as a simulation between the two transition systems:

* `Rel st' st` — `st'` (memory `k'`) and `st` (memory `k`) agree on the
  endpoint, direction, cross score and index, the history of `st` is a
  prefix of the history of `st'`, and `st'` has a larger certified lower
  and a smaller certified upper unplaced–unplaced contact count;
* `step_rel` — every step allowed from `st'` is allowed from `st`, and the
  successors are again related, whenever `k ≤ k'`;
* `spec_mono` — the relaxed path-sums are ordered for every nonnegative
  state function that is monotone along `Rel`;
* `memUppers_mono` — **refinement theorem**: for `k ≤ k'`, all three
  envelope upper bounds (partition function, event, complement) computed
  by the version-2 checker with memory `k'` are at most those with memory
  `k`.  More memory never widens a bound.
-/
@[expose] public section

namespace FoldBounds
namespace Refinement

open LatticeModel Frontier Transfer Geometry MemTransfer FrontierCheck FrontierNext Finset

/-- Simulation relation between a memory-`k'` state and a memory-`k` state. -/
def Rel (st' st : MState) : Prop :=
  st'.cur = st.cur ∧ st'.dir = st.dir ∧ st'.cross = st.cross ∧ st'.idx = st.idx ∧
    st.hist <+: st'.hist ∧ st.ffLo ≤ st'.ffLo ∧ st'.ffHi ≤ st.ffHi

theorem Rel_refl (st : MState) : Rel st st :=
  ⟨rfl, rfl, rfl, rfl, List.prefix_refl _, le_refl _, le_refl _⟩

theorem lookupH_prefix {h h' : List (Site × ℕ)} (hp : h <+: h') {i : ℕ} {x : Site × ℕ}
    (hx : lookupH h i = some x) : lookupH h' i = some x := by
  obtain ⟨t, rfl⟩ := hp
  unfold lookupH at hx ⊢
  rw [List.find?_append, hx]
  rfl

theorem loAt_mono (S : List Bool) (m : ℕ) {h h' : List (Site × ℕ)} (hp : h <+: h') (j : ℕ)
    (s : Site) : loAt S m h j s ≤ loAt S m h' j s := by
  refine Finset.sum_le_sum (fun i _ => ?_)
  cases hl : lookupH h i with
  | none => simp
  | some x =>
      rw [lookupH_prefix hp hl]

theorem lo_add_out_anti (S : List Bool) (m : ℕ) {h h' : List (Site × ℕ)} (hp : h <+: h')
    (j : ℕ) (s : Site) :
    loAt S m h' j s + outAt S m h' j ≤ loAt S m h j s + outAt S m h j := by
  unfold loAt outAt
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  refine Finset.sum_le_sum (fun i _ => ?_)
  cases hl : lookupH h i with
  | none =>
      cases hl' : lookupH h' i with
      | none => simp
      | some x =>
          simp only [Option.map_some, Option.getD_some, Option.isNone_some, Bool.and_false,
            Option.map_none, Option.getD_none, Option.isNone_none, Bool.and_true]
          split_ifs <;> simp_all
  | some x =>
      rw [lookupH_prefix hp hl]

theorem upAt_anti (S : List Bool) (m : ℕ) {h h' : List (Site × ℕ)} (hp : h <+: h') (j : ℕ)
    (s : Site) : upAt S m h' j s ≤ upAt S m h j s := by
  unfold upAt
  exact min_le_min (le_refl _) (lo_add_out_anti S m hp j s)

theorem take_cons_prefix {α : Type*} (a : α) {h h' : List α} (hp : h <+: h') {k k' : ℕ}
    (hk : k ≤ k') : (a :: h).take k <+: (a :: h').take k' :=
  (List.take_prefix_take_left (l := a :: h) hk).trans (List.IsPrefix.take (List.cons_prefix_cons.2 ⟨rfl, hp⟩) k')

theorem step_rel (S : List Bool) (Wp : List Site) (m : ℕ) {k k' : ℕ} (hk : k ≤ k')
    {st' st : MState} (hr : Rel st' st) (e : Fin 4) (ha : memAllowed Wp st' e = true) :
    memAllowed Wp st e = true ∧ Rel (memStep S Wp m k' st' e) (memStep S Wp m k st e) := by
  obtain ⟨hcur, hdir, hcross, hidx, hhist, hlo, hhi⟩ := hr
  refine ⟨?_, ?_⟩
  · simp only [memAllowed, Bool.and_eq_true, decide_eq_true_eq, Bool.not_eq_true',
      List.any_eq_false] at ha ⊢
    obtain ⟨⟨h1, h2⟩, h3⟩ := ha
    rw [hdir, hcur] at *
    refine ⟨⟨h1, h2⟩, fun x hx => h3 x (hhist.subset hx)⟩
  · simp only [Rel, memStep, hcur, hcross, hidx, true_and]
    refine ⟨take_cons_prefix _ hhist hk, ?_, ?_⟩
    · have := loAt_mono S m hhist st.idx (advance st.cur e); omega
    · have := upAt_anti S m hhist st.idx (advance st.cur e); omega

theorem spec_mono (S : List Bool) (Wp : List Site) (m : ℕ) {k k' : ℕ} (hk : k ≤ k')
    (H : MState → ℚ) (hH0 : ∀ st, 0 ≤ H st) (hH : ∀ st' st, Rel st' st → H st' ≤ H st) :
    ∀ r st' st, Rel st' st →
      spec (memNext S Wp m k') H r st' ≤ spec (memNext S Wp m k) H r st := by
  intro r
  induction r with
  | zero => intro st' st hr; exact hH _ _ hr
  | succ r ih =>
      intro st' st hr
      simp only [spec, memNext, List.map_map]
      rw [sum_map_filter_eq, sum_map_filter_eq, sum_allDirs, sum_allDirs]
      refine Finset.sum_le_sum (fun e _ => ?_)
      by_cases ha : memAllowed Wp st' e = true
      · obtain ⟨ha2, hr2⟩ := step_rel S Wp m hk hr e ha
        simp only [ha, ha2, if_true, Function.comp]
        exact ih _ _ hr2
      · simp only [ha, Bool.false_eq_true, if_false]
        split_ifs
        · exact spec_nonneg _ _ hH0 _ _
        · exact le_refl _

/-! ### Monotonicity of the geometric cap -/

theorem lmax_mono {l l' : List ℕ} (h : ∀ x ∈ l, x ∈ l') : lmax l ≤ lmax l' :=
  lmax_le (fun x hx => le_lmax (h x hx))

theorem xSpan_mono {K K' : List Site} (h : ∀ a ∈ K, a ∈ K') : xSpan K ≤ xSpan K' := by
  apply lmax_mono
  intro x hx
  simp only [List.mem_flatMap, List.mem_map] at hx ⊢
  obtain ⟨a, ha, b, hb, rfl⟩ := hx
  exact ⟨a, h a ha, b, h b hb, rfl⟩

theorem ySpan_mono {K K' : List Site} (h : ∀ a ∈ K, a ∈ K') : ySpan K ≤ ySpan K' := by
  apply lmax_mono
  intro x hx
  simp only [List.mem_flatMap, List.mem_map] at hx ⊢
  obtain ⟨a, ha, b, hb, rfl⟩ := hx
  exact ⟨a, h a ha, b, h b hb, rfl⟩

theorem le_foldr_min (f : ℕ → ℕ) (X init : ℕ) (l : List ℕ) (hi : X ≤ init)
    (hl : ∀ c ∈ l, X ≤ f c) : X ≤ l.foldr (fun c acc => min acc (f c)) init := by
  induction l with
  | nil => exact hi
  | cons a l ih =>
      simp only [List.foldr_cons]
      exact le_min (ih (fun c hc => hl c (by simp [hc]))) (hl a (by simp))

theorem foldr_min_le_init (f : ℕ → ℕ) (init : ℕ) (l : List ℕ) :
    l.foldr (fun c acc => min acc (f c)) init ≤ init := by
  induction l with
  | nil => exact le_refl _
  | cons a l ih => exact (min_le_left _ _).trans ih

theorem foldr_min_le_mem (f : ℕ → ℕ) (init : ℕ) (l : List ℕ) {c : ℕ} (hc : c ∈ l) :
    l.foldr (fun c acc => min acc (f c)) init ≤ f c := by
  induction l with
  | nil => simp at hc
  | cons a l ih =>
      simp only [List.foldr_cons]
      rcases List.mem_cons.1 hc with rfl | hc
      · exact min_le_right _ _
      · exact (min_le_left _ _).trans (ih hc)

theorem minPerim_mono (N : ℕ) {w h w' h' : ℕ} (hw : w ≤ w') (hh : h ≤ h') :
    minPerim N w h ≤ minPerim N w' h' := by
  unfold minPerim
  apply le_foldr_min _ _ _ _ (foldr_min_le_init _ _ _)
  intro c hc
  simp only [List.mem_filter, decide_eq_true_eq] at hc
  have hc' : c ∈ (List.range (N + 1)).filter (fun c => w ≤ c) := by
    simp only [List.mem_filter, decide_eq_true_eq]; exact ⟨hc.1, hw.trans hc.2⟩
  refine (foldr_min_le_mem _ _ _ hc').trans ?_
  exact Nat.add_le_add_left (max_le_max hh (le_refl _)) _

theorem geomCap_anti (N : ℕ) {K K' : List Site} (h : ∀ a ∈ K, a ∈ K') :
    geomCap N K' ≤ geomCap N K := by
  unfold geomCap
  have := minPerim_mono N (xSpan_mono h) (ySpan_mono h)
  omega

theorem knownSites_mono (p : List (Fin 4)) {st' st : MState} (hr : Rel st' st) :
    ∀ a ∈ knownSites p st, a ∈ knownSites p st' := by
  obtain ⟨hcur, _, _, _, hhist, _, _⟩ := hr
  intro a ha
  simp only [knownSites, List.mem_append, List.mem_cons, List.mem_map] at ha ⊢
  rcases ha with ha | rfl | ⟨x, hx, rfl⟩
  · exact Or.inl ha
  · exact Or.inr (Or.inl hcur.symm)
  · exact Or.inr (Or.inr ⟨x, hhist.subset hx, rfl⟩)

theorem mExp_anti (steps : ℕ) (S : List Bool) (p : List (Fin 4)) {st' st : MState}
    (hr : Rel st' st) : mExp steps S p st' ≤ mExp steps S p st := by
  have hg := geomCap_anti (steps + 2) (knownSites_mono p hr)
  obtain ⟨_, _, hcross, _, _, _, hhi⟩ := hr
  unfold mExp
  rw [hcross]
  exact min_le_min (min_le_min (by omega) (le_refl _)) hg

theorem mLow_mono (S : List Bool) (p : List (Fin 4)) {st' st : MState} (hr : Rel st' st) :
    mLow S p st ≤ mLow S p st' := by
  obtain ⟨_, _, hcross, _, _, hlo, _⟩ := hr
  unfold mLow; omega

/-- **Refinement theorem.** For `k ≤ k'`, every envelope upper bound of
the memory-`k'` leaf is at most the corresponding memory-`k` bound. -/
theorem memUppers_mono (steps : ℕ) (S : List Bool) {q : ℚ} (hq : 1 ≤ q) (thr : ℕ)
    (p : List (Fin 4)) {k k' : ℕ} (hk : k ≤ k') :
    (memUppers k' steps S q thr p).1 ≤ (memUppers k steps S q thr p).1 ∧
    (memUppers k' steps S q thr p).2.1 ≤ (memUppers k steps S q thr p).2.1 ∧
    (memUppers k' steps S q thr p).2.2 ≤ (memUppers k steps S q thr p).2.2 := by
  have hq0 : (0 : ℚ) ≤ q := by linarith
  have hpos : ∀ k : ℕ, (0 : ℚ) ≤ q ^ k := fun k => pow_nonneg hq0 k
  have tbl : ∀ H : MState → ℚ, (∀ st, 0 ≤ H st) → (∀ st' st, Rel st' st → H st' ≤ H st) →
      wsum H (memTable steps S k' p) ≤ wsum H (memTable steps S k p) := by
    intro H hH0 hH
    simp only [memTable, wsum_dpM, wsum_cons, wsum_nil, one_mul, add_zero]
    exact spec_mono S _ _ hk H hH0 hH _ _ _ (Rel_refl _)
  have mZ : ∀ st' st, Rel st' st → hZ steps S q p st' ≤ hZ steps S q p st :=
    fun st' st hr => qpow_mono hq (mExp_anti steps S p hr)
  have mA : ∀ st' st, Rel st' st → hA steps S q thr p st' ≤ hA steps S q thr p st := by
    intro st' st hr
    have := mExp_anti steps S p hr
    unfold hA
    split_ifs with h1 h2
    · exact qpow_mono hq this
    · omega
    · exact hpos _
    · exact le_refl _
  have mN : ∀ st' st, Rel st' st → hN steps S q thr p st' ≤ hN steps S q thr p st := by
    intro st' st hr
    have h1 := mExp_anti steps S p hr
    have h2 := mLow_mono S p hr
    unfold hN
    split_ifs with a b
    · exact qpow_mono hq h1
    · omega
    · exact hpos _
    · exact le_refl _
  have hZ0 : ∀ st, 0 ≤ hZ steps S q p st := fun _ => hpos _
  have hA0 : ∀ st, 0 ≤ hA steps S q thr p st := fun _ => by unfold hA; split_ifs <;> simp [hpos]
  have hN0 : ∀ st, 0 ≤ hN steps S q thr p st := fun _ => by unfold hN; split_ifs <;> simp [hpos]
  simp only [memUppers]
  exact ⟨min_le_min (le_refl _) (tbl _ hZ0 mZ), min_le_min (le_refl _) (tbl _ hA0 mA),
    min_le_min (le_refl _) (tbl _ hN0 mN)⟩

/-- **Never worse than the version-1 counting bounds.** The memory-leaf
envelope uppers are at most the version-1 basic (`tr = false`) and
degree-capped (`tr = true`) counting × cap bounds, for every memory. -/
theorem memUppers_le_basic (k steps : ℕ) (S : List Bool) {q : ℚ} (hq : 1 ≤ q) (thr : ℕ)
    (p : List (Fin 4)) (tr : Bool) :
    (memUppers k steps S q thr p).1 ≤ basicZ tr steps S q p ∧
    (memUppers k steps S q thr p).2.1 ≤ basicA tr steps S q thr p ∧
    (memUppers k steps S q thr p).2.2 ≤ basicN tr steps S q thr p := by
  have hq0 : (0 : ℚ) ≤ q := by linarith
  have hk : kcN steps S p ≤ kcap tr S p := by
    have h1 : kcN steps S p ≤ kcap true S p := by unfold kcN; exact min_le_left _ _
    have h2 := futureCap_le_futurePairs S (p.length + 2)
    cases tr
    · simp only [kcap, kmax, Bool.false_eq_true, if_false] at h1 ⊢
      simp only [if_true] at h1; omega
    · exact h1
  have hZ : capZ steps S q p ≤ basicZ tr steps S q p := by
    unfold capZ basicZ
    exact mul_le_mul_of_nonneg_left (qpow_mono hq hk) (by positivity)
  have hZ0 : 0 ≤ capZ steps S q p := by unfold capZ; positivity
  simp only [memUppers]
  refine ⟨(min_le_left _ _).trans hZ, (min_le_left _ _).trans ?_, (min_le_left _ _).trans ?_⟩
  · unfold capA basicA
    split_ifs with h1 h2
    · exact hZ
    · omega
    · exact hZ0.trans hZ
    · exact le_refl _
  · unfold capN basicN
    split_ifs
    · exact hZ
    · exact le_refl _

end Refinement
end FoldBounds
