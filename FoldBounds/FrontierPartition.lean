module

public import FoldBounds.FrontierSoundness

/-!
# Accepted frontiers partition the conformation space; transfer ≤ basic

* `sum_ext_eq_sum_leaves`: for an accepted trie, the completion set of its
  root is the *disjoint* union of the completion sets of its leaves (all
  `exact`, `bound` and `dead` entries), in the strong form that every
  completion sum over the root equals the sum over leaves of the
  completion sums.
* `leaf_cover_unique`: consequently every valid complete conformation lies
  in the completion set of exactly one frontier leaf (and dead leaves cover
  nothing).
* `relaxed_le_basic`: the transfer-relaxation envelope sum never exceeds the
  basic `3^(n-m) · C` envelope, so enabling `transfer` can only tighten the
  bound (the `min` in `envelopeUppers` is a redundant safety net).
-/
@[expose] public section

namespace FoldBounds
namespace FrontierCheck

open LatticeModel Frontier Transfer Finset

/-- Frontier leaves (prefixes of `exact`, `bound`, `dead` entries). -/
def leaves : List (Fin 4) → FTree → List (List (Fin 4))
  | p, .node t0 t1 t2 t3 =>
      leaves (p ++ [0]) t0 ++ leaves (p ++ [1]) t1 ++ leaves (p ++ [2]) t2 ++
        leaves (p ++ [3]) t3
  | p, .exact => [p]
  | p, .bound _ _ => [p]
  | p, .dead => [p]

/-- **Accepted frontiers partition the completion set.** -/
theorem sum_ext_eq_sum_leaves {M : Type*} [AddCommMonoid M] (steps : ℕ) (S : List Bool)
    (q : ℚ) (thr : ℕ) (F : List (Fin 4) → M) :
    ∀ (t : FTree) (p : List (Fin 4)) (e : Enclosure),
      checkNode steps S q thr p t = some e →
      ∑ l ∈ ext steps p, F l = ((leaves p t).map (fun x => ∑ l ∈ ext steps x, F l)).sum := by
  intro t
  induction t with
  | exact => intro p e _; simp [leaves]
  | bound _ _ => intro p e _; simp [leaves]
  | dead => intro p e _; simp [leaves]
  | node t0 t1 t2 t3 ih0 ih1 ih2 ih3 =>
      intro p e h
      simp only [checkNode] at h
      split_ifs at h with hc
      cases h0 : checkNode steps S q thr (p ++ [0]) t0 <;>
      cases h1 : checkNode steps S q thr (p ++ [1]) t1 <;>
      cases h2 : checkNode steps S q thr (p ++ [2]) t2 <;>
      cases h3 : checkNode steps S q thr (p ++ [3]) t3 <;>
      simp only [h0, h1, h2, h3, reduceCtorEq] at h
      rw [sum_ext_split hc.1, Fin.sum_univ_four, ih0 _ _ h0, ih1 _ _ h1, ih2 _ _ h2,
        ih3 _ _ h3]
      simp only [leaves, List.map_append, List.sum_append]

/-- **Unique coverage.** In an accepted certificate every valid complete
conformation belongs to the completion set of exactly one frontier leaf. -/
theorem leaf_cover_unique (S : List Bool) (q : ℚ) (thr : ℕ) (t : FTree) (e : Enclosure)
    (h : checkNode (S.length - 2) S q thr [] t = some e)
    (l : List (Fin 4)) (hl : l ∈ validWords (S.length - 2)) :
    ((leaves [] t).map (fun x => if l ∈ ext (S.length - 2) x then 1 else 0)).sum = 1 := by
  have key := sum_ext_eq_sum_leaves (S.length - 2) S q thr
    (fun l' => if l' = l then (1 : ℕ) else 0) t [] e h
  rw [ext_nil, Finset.sum_ite_eq' _ l, if_pos hl] at key
  have hmap : (leaves [] t).map (fun x => ∑ l' ∈ ext (S.length - 2) x,
      (fun l' => if l' = l then (1 : ℕ) else 0) l') =
      (leaves [] t).map (fun x => if l ∈ ext (S.length - 2) x then 1 else 0) :=
    List.map_congr_left (fun x _ => Finset.sum_ite_eq' _ _ _)
  rw [hmap] at key
  exact key.symm

/-! ### Transfer never loses against the basic bound -/

theorem length_relaxedNext_le (S : List Bool) (Wp : List Site) (m : ℕ) (st : RState) :
    (relaxedNext S Wp m st).length ≤ 3 := by
  unfold relaxedNext allowed
  rw [List.length_map]
  have h : (allDirs.filter (fun e => decide (e ≠ opp st.2.1) &&
      !(Wp.contains (advance st.1 e)))) =
      (allDirs.filter (fun e => decide (e ≠ opp st.2.1))).filter
        (fun e => !(Wp.contains (advance st.1 e))) := by
    rw [List.filter_filter]
    congr 1
    funext e
    exact Bool.and_comm _ _
  rw [h]
  refine (List.length_filter_le _ _).trans ?_
  have key : ∀ d : Fin 4, (allDirs.filter (fun e => decide (e ≠ opp d))).length ≤ 3 := by
    decide
  exact key _

theorem spec_le_pow {σ : Type*} (next : σ → List σ) (G : σ → ℚ) (C : ℚ) (hC : 0 ≤ C)
    (hG : ∀ s, G s ≤ C) (hlen : ∀ s, (next s).length ≤ 3) :
    ∀ r s, spec next G r s ≤ 3 ^ r * C := by
  intro r
  induction r with
  | zero => intro s; simpa using hG s
  | succ r ih =>
      intro s
      rw [spec]
      calc ((next s).map (spec next G r)).sum
          ≤ ((next s).map (spec next G r)).length • (3 ^ r * C) :=
            List.sum_le_card_nsmul _ _ (by
              intro x hx
              obtain ⟨y, _, rfl⟩ := List.mem_map.1 hx
              exact ih y)
        _ = ((next s).length : ℚ) * (3 ^ r * C) := by simp [nsmul_eq_mul]
        _ ≤ 3 * (3 ^ r * C) := by
            apply mul_le_mul_of_nonneg_right _ (by positivity)
            exact_mod_cast hlen s
        _ = 3 ^ (r + 1) * C := by ring

/-- The relaxed transfer table never exceeds the basic counting envelope. -/
theorem relaxed_le_basic (steps : ℕ) (S : List Bool) (p : List (Fin 4)) (g : ℕ → ℚ)
    (C : ℚ) (hC : 0 ≤ C) (hg : ∀ k, g k ≤ C) :
    wsum (fun st => g st.2.2.1) (relaxedTable steps S p) ≤ 3 ^ (steps - p.length) * C := by
  rw [relaxedTable, wsum_dp, initTable, wsum_cons, wsum_nil, one_mul, add_zero]
  exact spec_le_pow _ _ C hC (fun st => hg _) (length_relaxedNext_le S _ _) _ _

/-- The degree-capped counting envelope is never larger than Python's
basic envelope `3^(n-m) q^{K_max}`. -/
theorem basicZ_degree_le_basic (steps : ℕ) (S : List Bool) {q : ℚ} (hq : 1 ≤ q)
    (p : List (Fin 4)) : basicZ true steps S q p ≤ basicZ false steps S q p := by
  unfold basicZ kcap kmax
  simp only [if_true, Bool.false_eq_true, if_false]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  exact qpow_mono hq (by have := futureCap_le_futurePairs S (p.length + 2); omega)

/-- **Hierarchy theorem.** The transfer-relaxation total-mass envelope is
never larger than the degree-capped counting envelope, hence never larger
than the basic `3^(n-m) q^{K_max}` envelope. -/
theorem transfer_total_le_basic (steps : ℕ) (S : List Bool) {q : ℚ} (hq : 1 ≤ q)
    (p : List (Fin 4)) :
    wsum (fun st => gZ true S q p st.2.2.1) (relaxedTable steps S p) ≤
      basicZ false steps S q p := by
  refine le_trans ?_ (basicZ_degree_le_basic steps S hq p)
  unfold basicZ
  exact relaxed_le_basic steps S p _ _ (pow_nonneg (by linarith) _)
    (fun k => qpow_mono hq (envExp_le_kcap _ _ _ _))

end FrontierCheck
end FoldBounds
