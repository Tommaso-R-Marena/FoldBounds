module

public import FoldBounds.Frontier

/-!
# Certified conformation families: non-enumerative lower bounds

If the last site `e` of a valid prefix `p` maximises a diagonal potential
`φ(x, y) = ±x ± y` over all prefix sites, then *every* continuation using
only the two steps that increase `φ` by one is self-avoiding: `φ` is
strictly increasing along the new sites, and every new site lies strictly
above every prefix site.  These `2^r` continuations (`r` remaining steps)
are distinct genuine completions of `p`, so

  `∑_{c ∈ C(p)} F(c) ≥ 2^r · min_{c ∈ C(p)} F(c)`

for any nonnegative `F`, without listing a single completion.  With
`F = q^k` and `k(c) ≥ known(p)` this gives the lower bound
`2^r q^{known(p)}` used by the checker.
-/
@[expose] public section

namespace FoldBounds
namespace Frontier

open LatticeModel Finset

/-- The `x` step of a diagonal family (`true`: `+x`, `false`: `-x`). -/
def dirX (a : Bool) : Fin 4 := if a then 0 else 1

/-- The `y` step of a diagonal family (`true`: `+y`, `false`: `-y`). -/
def dirY (b : Bool) : Fin 4 := if b then 2 else 3

/-- Diagonal potential `(±x) + (±y)`. -/
def pot (a b : Bool) (s : Site) : ℤ :=
  (if a then s.1 else -s.1) + (if b then s.2 else -s.2)

theorem pot_advance_dirX (a b : Bool) (s : Site) :
    pot a b (advance s (dirX a)) = pot a b s + 1 := by
  cases a <;> cases b <;> simp [pot, advance, dirX] <;> ring

theorem pot_advance_dirY (a b : Bool) (s : Site) :
    pot a b (advance s (dirY b)) = pot a b s + 1 := by
  cases a <;> cases b <;> simp [pot, advance, dirY] <;> ring

theorem dirX_ne_dirY (a b : Bool) : dirX a ≠ dirY b := by
  cases a <;> cases b <;> decide

/-- Along a diagonal-family word the potential strictly increases. -/
theorem trace_diag (a b : Bool) :
    ∀ (w : List (Fin 4)), (∀ d ∈ w, d = dirX a ∨ d = dirY b) → ∀ e : Site,
      (trace e w).Nodup ∧ ∀ s ∈ trace e w, pot a b e < pot a b s := by
  intro w
  induction w with
  | nil => intro _ e; simp [trace]
  | cons d w ih =>
      intro hw e
      have hd : pot a b (advance e d) = pot a b e + 1 := by
        rcases hw d (by simp) with h | h
        · rw [h, pot_advance_dirX]
        · rw [h, pot_advance_dirY]
      obtain ⟨hn, hlt⟩ := ih (fun x hx => hw x (by simp [hx])) (advance e d)
      simp only [trace]
      refine ⟨List.nodup_cons.2 ⟨fun hm => ?_, hn⟩, ?_⟩
      · have := hlt _ hm; omega
      · intro s hs
        rcases List.mem_cons.1 hs with rfl | hs
        · omega
        · have := hlt s hs; omega

/-- Locally checkable family condition: the prefix end maximises the
diagonal potential over the prefix walk. -/
def diagOK (a b : Bool) (p : List (Fin 4)) : Bool :=
  (walkFrom p).all (fun s => decide (pot a b s ≤ pot a b (endFrom first p)))

/-- Some diagonal family is available at `p`. -/
def famOK (p : List (Fin 4)) : Bool :=
  diagOK true true p || diagOK true false p || diagOK false true p || diagOK false false p

/-- A family completion is a valid self-avoiding walk. -/
theorem nodup_family {a b : Bool} {p : List (Fin 4)} (hn : (walkFrom p).Nodup)
    (hd : diagOK a b p = true) (w : List (Fin 4)) (hw : ∀ d ∈ w, d = dirX a ∨ d = dirY b) :
    (walkFrom (p ++ w)).Nodup := by
  rw [walkFrom_append]
  obtain ⟨htn, hlt⟩ := trace_diag a b w hw (endFrom first p)
  refine List.nodup_append.2 ⟨hn, htn, ?_⟩
  intro x hx y hy hxy
  subst hxy
  have h1 := of_decide_eq_true (List.all_eq_true.1 hd x hx)
  have h2 := hlt x hy
  omega

/-- The family word selected by a Boolean pattern. -/
def famWord (a b : Bool) {r : ℕ} (f : Fin r → Bool) : List (Fin 4) :=
  List.ofFn (fun i => if f i then dirX a else dirY b)

theorem famWord_injective (a b : Bool) (r : ℕ) :
    Function.Injective (famWord a b (r := r)) := by
  intro f g h
  funext i
  have := congrArg (fun l => l[i.1]?) h
  simp only [famWord, List.getElem?_ofFn] at this
  have hne := dirX_ne_dirY a b
  simp only [i.2, dite_true, Option.some.injEq, Fin.eta] at this
  cases hf : f i <;> cases hg : g i
  · rfl
  · rw [hf, hg] at this; exact absurd this.symm (by simpa using hne)
  · rw [hf, hg] at this; exact absurd this (by simpa using hne)
  · rfl

/-- The family as a finite set of complete words. -/
def family (a b : Bool) (p : List (Fin 4)) (r : ℕ) : Finset (List (Fin 4)) :=
  (Finset.univ : Finset (Fin r → Bool)).image (fun f => p ++ famWord a b f)

theorem card_family (a b : Bool) (p : List (Fin 4)) (r : ℕ) :
    (family a b p r).card = 2 ^ r := by
  rw [family, Finset.card_image_of_injective]
  · simp
  · intro f g h
    exact famWord_injective a b r (List.append_cancel_left h)

theorem family_subset_ext {steps : ℕ} {a b : Bool} {p : List (Fin 4)}
    (hn : (walkFrom p).Nodup) (hd : diagOK a b p = true) (hle : p.length ≤ steps) :
    family a b p (steps - p.length) ⊆ ext steps p := by
  intro l hl
  obtain ⟨f, _, rfl⟩ := Finset.mem_image.1 hl
  rw [mem_ext]
  refine ⟨by simp [famWord]; omega, nodup_family hn hd _ ?_, List.prefix_append _ _⟩
  intro d hdm
  simp only [famWord, List.mem_ofFn] at hdm
  obtain ⟨i, rfl⟩ := hdm
  by_cases h : f i <;> simp [h]

/-- **Family lower bound.** If some diagonal family is available at a
valid prefix `p`, then for every nonnegative `F` bounded below by `c` on the
completion set, `2^(steps-|p|) · c ≤ ∑_{l ∈ C(p)} F l`. -/
theorem family_lower {steps : ℕ} {p : List (Fin 4)} (hn : (walkFrom p).Nodup)
    (hle : p.length ≤ steps) (hf : famOK p = true) (F : List (Fin 4) → ℚ)
    (hF : ∀ l, 0 ≤ F l) (c : ℚ) (hc : ∀ l ∈ ext steps p, c ≤ F l) :
    2 ^ (steps - p.length) * c ≤ ∑ l ∈ ext steps p, F l := by
  have key : ∀ a b, diagOK a b p = true →
      2 ^ (steps - p.length) * c ≤ ∑ l ∈ ext steps p, F l := by
    intro a b hd
    have hsub := family_subset_ext (steps := steps) hn hd hle
    calc (2 : ℚ) ^ (steps - p.length) * c
        = ∑ _l ∈ family a b p (steps - p.length), c := by
          rw [Finset.sum_const, card_family, nsmul_eq_mul]; push_cast; ring
      _ ≤ ∑ l ∈ family a b p (steps - p.length), F l :=
          Finset.sum_le_sum (fun l hl => hc l (hsub hl))
      _ ≤ ∑ l ∈ ext steps p, F l :=
          Finset.sum_le_sum_of_subset_of_nonneg hsub (fun l _ _ => hF l)
  unfold famOK at hf
  simp only [Bool.or_eq_true] at hf
  rcases hf with ((h | h) | h) | h
  · exact key _ _ h
  · exact key _ _ h
  · exact key _ _ h
  · exact key _ _ h

end Frontier
end FoldBounds
