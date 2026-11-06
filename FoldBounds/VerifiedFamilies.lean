module

public import FoldBounds.FrontierFamily

/-!
# Partially directed constructive families

If the endpoint of a valid prefix `p` is the *unique* site of the prefix
walk with maximal coordinate along an axis direction `m` (`pdOK m p`),
then every *partially directed* continuation — a word over
`{m, up, down}` (the main direction and its two perpendiculars) with no
immediate `up`/`down` reversal — is a valid self-avoiding completion of
`p`.  The number of such words of length `r` is given by an exact linear
recurrence (`pdTriple`), growing like `(1 + √2)^r ≈ 2.414^r` and always at
least `2^r` (the diagonal family of `FrontierFamily` has `2^r` members).

Main results:

* `trace_pd` — partially directed words trace self-avoiding paths that
  stay strictly ahead of (or beside, in the starting column) their start;
* `pdFamily_subset_ext` — the family is contained in the completion set;
* `card_pdFamily` — its cardinality is exactly `pdCount r` (computed by
  the recurrence, never by enumeration);
* `mem_pdFamily_test` — membership implies the executable test `pdTest`,
  used by the checker to keep witness lower bounds disjoint from the family.
-/
@[expose] public section

namespace FoldBounds
namespace Families

open LatticeModel Frontier Finset

/-! ### Orientations -/

/-- Perpendicular "up" direction for main direction `m`. -/
def oUp (m : Fin 4) : Fin 4 := if m.val < 2 then 2 else 0

/-- Perpendicular "down" direction for main direction `m`. -/
def oDown (m : Fin 4) : Fin 4 := if m.val < 2 then 3 else 1

/-- Coordinate along the main direction. -/
def oF (m : Fin 4) (s : Site) : ℤ :=
  if m.val = 0 then s.1 else if m.val = 1 then -s.1 else if m.val = 2 then s.2 else -s.2

/-- Perpendicular coordinate. -/
def oG (m : Fin 4) (s : Site) : ℤ := if m.val < 2 then s.2 else s.1

theorem oF_main (m : Fin 4) (s : Site) : oF m (advance s m) = oF m s + 1 := by
  fin_cases m <;> simp [oF, advance] <;> ring

theorem oG_main (m : Fin 4) (s : Site) : oG m (advance s m) = oG m s := by
  fin_cases m <;> simp [oG, advance]

theorem oF_up (m : Fin 4) (s : Site) : oF m (advance s (oUp m)) = oF m s := by
  fin_cases m <;> simp [oF, oUp, advance]

theorem oG_up (m : Fin 4) (s : Site) : oG m (advance s (oUp m)) = oG m s + 1 := by
  fin_cases m <;> simp [oG, oUp, advance]

theorem oF_down (m : Fin 4) (s : Site) : oF m (advance s (oDown m)) = oF m s := by
  fin_cases m <;> simp [oF, oDown, advance]

theorem oG_down (m : Fin 4) (s : Site) : oG m (advance s (oDown m)) = oG m s - 1 := by
  fin_cases m <;> simp [oG, oDown, advance]

theorem main_ne_up (m : Fin 4) : m ≠ oUp m := by fin_cases m <;> decide
theorem main_ne_down (m : Fin 4) : m ≠ oDown m := by fin_cases m <;> decide
theorem up_ne_down (m : Fin 4) : oUp m ≠ oDown m := by fin_cases m <;> decide

/-! ### Partially directed words -/

/-- Partially directed words of length `r`; context `c`: `0` free, `1` last
step was `up`, `2` last step was `down`. -/
def pdWords (m : Fin 4) : ℕ → ℕ → List (List (Fin 4))
  | 0, _ => [[]]
  | r + 1, c =>
      (pdWords m r 0).map (m :: ·) ++
      (if c ≠ 2 then (pdWords m r 1).map (oUp m :: ·) else []) ++
      (if c ≠ 1 then (pdWords m r 2).map (oDown m :: ·) else [])

/-- Counts `(free, after-up, after-down)` by the linear recurrence. -/
def pdTriple : ℕ → ℕ × ℕ × ℕ
  | 0 => (1, 1, 1)
  | r + 1 =>
      let t := pdTriple r
      (t.1 + t.2.1 + t.2.2, t.1 + t.2.1, t.1 + t.2.2)

/-- Number of partially directed continuations of length `r`. -/
def pdCount (r : ℕ) : ℕ := (pdTriple r).1

/-- Component of the triple for context `c`. -/
def tripleAt (t : ℕ × ℕ × ℕ) (c : ℕ) : ℕ := if c = 1 then t.2.1 else if c = 2 then t.2.2 else t.1

theorem length_pdWords (m : Fin 4) : ∀ r c, c ≤ 2 →
    (pdWords m r c).length = tripleAt (pdTriple r) c := by
  intro r
  induction r with
  | zero => intro c _; simp [pdWords, pdTriple, tripleAt]
  | succ r ih =>
      intro c hc
      simp only [pdWords, List.length_append, List.length_map]
      rw [ih 0 (by omega)]
      have h1 := ih 1 (by omega)
      have h2 := ih 2 (by omega)
      simp only [pdTriple, tripleAt] at h1 h2 ⊢
      interval_cases c <;> simp_all

theorem mem_pdWords_length (m : Fin 4) : ∀ r c w, w ∈ pdWords m r c → w.length = r := by
  intro r
  induction r with
  | zero => intro c w hw; simp [pdWords] at hw; simp [hw]
  | succ r ih =>
      intro c w hw
      simp only [pdWords, List.mem_append, List.mem_map] at hw
      rcases hw with (⟨w', h, rfl⟩ | hw) | hw
      · simp [ih _ _ h]
      · split_ifs at hw
        · obtain ⟨w', h, rfl⟩ := List.mem_map.1 hw; simp [ih _ _ h]
        · simp at hw
      · split_ifs at hw
        · obtain ⟨w', h, rfl⟩ := List.mem_map.1 hw; simp [ih _ _ h]
        · simp at hw

theorem nodup_pdWords (m : Fin 4) : ∀ r c, (pdWords m r c).Nodup := by
  intro r
  induction r with
  | zero => intro c; simp [pdWords]
  | succ r ih =>
      intro c
      have inj : ∀ d : Fin 4, Function.Injective (fun w : List (Fin 4) => d :: w) :=
        fun d a b h => by simpa using h
      have hm := (ih 0).map (inj m)
      have hu : (if c ≠ 2 then (pdWords m r 1).map (oUp m :: ·) else []).Nodup := by
        split_ifs
        · exact (ih 1).map (inj _)
        · exact List.nodup_nil
      have hd : (if c ≠ 1 then (pdWords m r 2).map (oDown m :: ·) else []).Nodup := by
        split_ifs
        · exact (ih 2).map (inj _)
        · exact List.nodup_nil
      have head_of : ∀ (d : Fin 4) (L : List (List (Fin 4))) (P : Prop) [Decidable P] (x),
          x ∈ (if P then L.map (d :: ·) else []) → x.head? = some d := by
        intro d L P _ x hx
        split_ifs at hx
        · obtain ⟨w, _, rfl⟩ := List.mem_map.1 hx; rfl
        · simp at hx
      simp only [pdWords]
      refine List.nodup_append.2 ⟨List.nodup_append.2 ⟨hm, hu, ?_⟩, hd, ?_⟩
      · intro x hx y hy hxy
        subst hxy
        obtain ⟨w, _, rfl⟩ := List.mem_map.1 hx
        have := head_of _ _ _ _ hy
        simp at this
        exact main_ne_up m this
      · intro x hx y hy hxy
        subst hxy
        have h2 := head_of _ _ _ _ hy
        rcases List.mem_append.1 hx with hx | hx
        · obtain ⟨w, _, rfl⟩ := List.mem_map.1 hx
          simp at h2; exact main_ne_down m h2
        · have h1 := head_of _ _ _ _ hx
          rw [h1] at h2; simp at h2; exact up_ne_down m h2

/-- Admissible offsets in the starting column, by context. -/
def colOK (c : ℕ) (x : ℤ) : Prop := if c = 1 then 0 < x else if c = 2 then x < 0 else x ≠ 0

/-- **Self-avoidance of partially directed words.** -/
theorem trace_pd (m : Fin 4) : ∀ r c w, c ≤ 2 → w ∈ pdWords m r c → ∀ s : Site,
    (trace s w).Nodup ∧ ∀ t ∈ trace s w, oF m s < oF m t ∨
      (oF m t = oF m s ∧ colOK c (oG m t - oG m s)) := by
  intro r
  induction r with
  | zero => intro c w _ hw s; simp [pdWords] at hw; subst hw; simp [trace]
  | succ r ih =>
      intro c w hc hw s
      simp only [pdWords, List.mem_append, List.mem_map] at hw
      rcases hw with (⟨w', hw', rfl⟩ | hw) | hw
      · obtain ⟨hn, ht⟩ := ih 0 w' (by omega) hw' (advance s m)
        have hF := oF_main m s
        simp only [trace]
        refine ⟨List.nodup_cons.2 ⟨fun hm => ?_, hn⟩, ?_⟩
        · rcases ht _ hm with h | ⟨_, h⟩
          · omega
          · simp [colOK] at h
        · intro t ht'
          rcases List.mem_cons.1 ht' with rfl | ht'
          · left; omega
          · left; rcases ht t ht' with h | ⟨h, _⟩ <;> omega
      · split_ifs at hw with hc2
        · obtain ⟨w', hw', rfl⟩ := List.mem_map.1 hw
          obtain ⟨hn, ht⟩ := ih 1 w' (by omega) hw' (advance s (oUp m))
          have hF := oF_up m s
          have hG := oG_up m s
          simp only [trace]
          refine ⟨List.nodup_cons.2 ⟨fun hm => ?_, hn⟩, ?_⟩
          · rcases ht _ hm with h | ⟨_, h⟩
            · omega
            · simp [colOK] at h
          · intro t ht'
            rcases List.mem_cons.1 ht' with rfl | ht'
            · right; refine ⟨hF, ?_⟩; simp only [colOK]; split_ifs <;> omega
            · rcases ht t ht' with h | ⟨h1, h2⟩
              · left; omega
              · right; refine ⟨by omega, ?_⟩
                simp only [colOK, if_true] at h2
                simp only [colOK]; split_ifs <;> omega
        · simp at hw
      · split_ifs at hw with hc1
        · obtain ⟨w', hw', rfl⟩ := List.mem_map.1 hw
          obtain ⟨hn, ht⟩ := ih 2 w' (by omega) hw' (advance s (oDown m))
          have hF := oF_down m s
          have hG := oG_down m s
          simp only [trace]
          refine ⟨List.nodup_cons.2 ⟨fun hm => ?_, hn⟩, ?_⟩
          · rcases ht _ hm with h | ⟨_, h⟩
            · omega
            · simp [colOK] at h
          · intro t ht'
            rcases List.mem_cons.1 ht' with rfl | ht'
            · right; refine ⟨hF, ?_⟩; simp only [colOK]; split_ifs <;> omega
            · rcases ht t ht' with h | ⟨h1, h2⟩
              · left; omega
              · right; refine ⟨by omega, ?_⟩
                simp only [colOK, show (2 : ℕ) ≠ 1 by decide, if_false, if_true] at h2
                simp only [colOK]; split_ifs <;> omega
        · simp at hw

/-! ### The family at a prefix -/

/-- The endpoint of `p` is the unique prefix site maximising `oF m`. -/
def pdOK (m : Fin 4) (p : List (Fin 4)) : Bool :=
  (walkFrom p).all (fun s => s == endFrom first p || decide (oF m s < oF m (endFrom first p)))

/-- The family as a finite set of complete words. -/
def pdFamily (m : Fin 4) (p : List (Fin 4)) (r : ℕ) : Finset (List (Fin 4)) :=
  ((pdWords m r 0).map (p ++ ·)).toFinset

theorem card_pdFamily (m : Fin 4) (p : List (Fin 4)) (r : ℕ) :
    (pdFamily m p r).card = pdCount r := by
  rw [pdFamily, List.toFinset_card_of_nodup ((nodup_pdWords m r 0).map
    (List.append_right_injective p)), List.length_map, length_pdWords m r 0 (by omega)]
  rfl

theorem pdFamily_subset_ext {steps : ℕ} {m : Fin 4} {p : List (Fin 4)}
    (hn : (walkFrom p).Nodup) (hd : pdOK m p = true) (hle : p.length ≤ steps) :
    pdFamily m p (steps - p.length) ⊆ ext steps p := by
  intro l hl
  simp only [pdFamily, List.mem_toFinset, List.mem_map] at hl
  obtain ⟨w, hw, rfl⟩ := hl
  rw [mem_ext]
  refine ⟨by simp [mem_pdWords_length m _ _ _ hw]; omega, ?_, List.prefix_append _ _⟩
  rw [walkFrom_append]
  obtain ⟨htn, hlt⟩ := trace_pd m _ 0 w (by omega) hw (endFrom first p)
  refine List.nodup_append.2 ⟨hn, htn, ?_⟩
  intro x hx y hy hxy
  subst hxy
  have h1 := List.all_eq_true.1 hd x hx
  simp only [Bool.or_eq_true, beq_iff_eq, decide_eq_true_eq] at h1
  rcases hlt x hy with h2 | ⟨h2, h3⟩
  · rcases h1 with rfl | h1 <;> omega
  · rcases h1 with rfl | h1
    · simp [colOK] at h3
    · omega

/-- Executable membership test for partially directed words. -/
def pdTest (m : Fin 4) : ℕ → List (Fin 4) → Bool
  | _, [] => true
  | c, d :: w =>
      if d = m then pdTest m 0 w
      else if d = oUp m ∧ c ≠ 2 then pdTest m 1 w
      else if d = oDown m ∧ c ≠ 1 then pdTest m 2 w
      else false

theorem pdTest_of_mem (m : Fin 4) : ∀ r c w, w ∈ pdWords m r c → pdTest m c w = true := by
  intro r
  induction r with
  | zero => intro c w hw; simp [pdWords] at hw; subst hw; rfl
  | succ r ih =>
      intro c w hw
      simp only [pdWords, List.mem_append, List.mem_map] at hw
      rcases hw with (⟨w', hw', rfl⟩ | hw) | hw
      · simp [pdTest, ih _ _ hw']
      · split_ifs at hw with hc2
        · obtain ⟨w', hw', rfl⟩ := List.mem_map.1 hw
          simp [pdTest, (main_ne_up m).symm, hc2, ih _ _ hw']
        · simp at hw
      · split_ifs at hw with hc1
        · obtain ⟨w', hw', rfl⟩ := List.mem_map.1 hw
          simp [pdTest, (main_ne_down m).symm, (up_ne_down m).symm, hc1, ih _ _ hw']
        · simp at hw

theorem mem_pdFamily_test {m : Fin 4} {p w : List (Fin 4)} {r : ℕ}
    (h : p ++ w ∈ pdFamily m p r) : pdTest m 0 w = true := by
  simp only [pdFamily, List.mem_toFinset, List.mem_map] at h
  obtain ⟨w', hw', he⟩ := h
  have : w' = w := List.append_cancel_left he
  subst this
  exact pdTest_of_mem m _ _ _ hw'

/-- The partially directed family is never smaller than the diagonal one. -/
theorem two_pow_le_pdCount (r : ℕ) : 2 ^ r ≤ pdCount r := by
  have key : ∀ r, 2 ^ r ≤ (pdTriple r).1 ∧ 2 ^ r ≤ (pdTriple r).2.1 + (pdTriple r).2.2 - 2 ^ r
      ∧ 2 ^ r ≤ (pdTriple r).2.1 ∧ 2 ^ r ≤ (pdTriple r).2.2 := by
    intro r
    induction r with
    | zero => simp [pdTriple]
    | succ r ih =>
        simp only [pdTriple, pow_succ]
        omega
  exact (key r).1

end Families
end FoldBounds
