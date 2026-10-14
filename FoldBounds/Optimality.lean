module

public import FoldBounds.Probability

/-!
# FoldBounds.Optimality — the block-certificate enclosure of `P(A)` is data-optimal

`prob_mem_of_blocks` encloses `P(A)` using only the *certificate data*
`(S, E|_S, blk, ℓ, u)`: exact energies on the explored set and block energy intervals on
the rest.  Here we show both endpoints are **attained** by energy landscapes `E'` that are
consistent with exactly the same data.  Consequently no sound procedure whose output
depends only on that data can return a strictly smaller interval: tightening requires new
information (more enumeration, finer blocks, or sharper block intervals) — which is
precisely what the refinement theorems guarantee will help.
-/

@[expose] public section

open Finset Real

namespace FoldBounds

variable {C : Type*} [Fintype C] [DecidableEq C] {ι : Type*}

/-- Landscape that agrees with `E` on `S` and puts every unexplored conformation of `A`
at energy `eA (blk c)` and every unexplored conformation of `Aᶜ` at `eB (blk c)`. -/
noncomputable def extremalLandscape (E : C → ℝ) (blk : C → ι) (eA eB : ι → ℝ)
    (S A : Finset C) : C → ℝ :=
  fun c => if c ∈ S then E c else if c ∈ A then eA (blk c) else eB (blk c)

omit [Fintype C] in
private theorem ZOn_extremal_inter (β : ℝ) (E : C → ℝ) (blk : C → ι) (eA eB : ι → ℝ)
    (S A B : Finset C) :
    ZOn β (extremalLandscape E blk eA eB S A) (B ∩ S) = ZOn β E (B ∩ S) :=
  Finset.sum_congr rfl fun c hc => by
    simp [weight, extremalLandscape, (Finset.mem_inter.mp hc).2]

omit [Fintype C] in
private theorem ZOn_extremal_A (β : ℝ) (E : C → ℝ) (blk : C → ι) (eA eB : ι → ℝ)
    (S A : Finset C) :
    ZOn β (extremalLandscape E blk eA eB S A) (A \ S) = blockBound β blk eA (A \ S) :=
  Finset.sum_congr rfl fun c hc => by
    obtain ⟨h1, h2⟩ := Finset.mem_sdiff.mp hc
    simp [weight, extremalLandscape, h1, h2]

private theorem ZOn_extremal_Ac (β : ℝ) (E : C → ℝ) (blk : C → ι) (eA eB : ι → ℝ)
    (S A : Finset C) :
    ZOn β (extremalLandscape E blk eA eB S A) (Aᶜ \ S) = blockBound β blk eB (Aᶜ \ S) :=
  Finset.sum_congr rfl fun c hc => by
    obtain ⟨h1, h2⟩ := Finset.mem_sdiff.mp hc
    have h1' : c ∉ A := Finset.mem_compl.mp h1
    simp [weight, extremalLandscape, h1', h2]

private theorem prob_extremal (β : ℝ) (E : C → ℝ) (blk : C → ι) (eA eB : ι → ℝ)
    (S A : Finset C) :
    prob β (extremalLandscape E blk eA eB S A) A =
      (ZOn β E (A ∩ S) + blockBound β blk eA (A \ S)) /
        ((ZOn β E (A ∩ S) + blockBound β blk eA (A \ S)) +
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk eB (Aᶜ \ S))) := by
  have hsplit : ∀ B : Finset C, ZOn β (extremalLandscape E blk eA eB S A) B =
      ZOn β (extremalLandscape E blk eA eB S A) (B ∩ S) +
        ZOn β (extremalLandscape E blk eA eB S A) (B \ S) :=
    fun B => (Finset.sum_inter_add_sum_diff B S _).symm
  rw [prob_eq_ratio, hsplit A, hsplit Aᶜ, ZOn_extremal_inter, ZOn_extremal_inter, ZOn_extremal_A, ZOn_extremal_Ac]

omit [Fintype C] [DecidableEq C] in
private theorem extremal_consistent (E : C → ℝ) (blk : C → ι) (ℓ u eA eB : ι → ℝ)
    (S A : Finset C) [DecidableEq C]
    (hA : ∀ c, ℓ (blk c) ≤ eA (blk c) ∧ eA (blk c) ≤ u (blk c) ∨ c ∈ S)
    (hB : ∀ c, ℓ (blk c) ≤ eB (blk c) ∧ eB (blk c) ≤ u (blk c) ∨ c ∈ S) :
    (∀ c ∈ S, extremalLandscape E blk eA eB S A c = E c) ∧
      (∀ c ∉ S, ℓ (blk c) ≤ extremalLandscape E blk eA eB S A c ∧
        extremalLandscape E blk eA eB S A c ≤ u (blk c)) := by
  refine ⟨fun c hc => by simp [extremalLandscape, hc], fun c hc => ?_⟩
  simp only [extremalLandscape, hc, if_false]
  split_ifs
  · exact (hA c).resolve_right hc
  · exact (hB c).resolve_right hc

/-- **Data-optimality of the lower endpoint.** For any valid block certificate there is an
energy landscape `E'` with the *same* certificate data (equal to `E` on `S`, within
`[ℓ, u]` blockwise off `S`) whose probability `P_{E'}(A)` equals the certified lower
endpoint of `prob_mem_of_blocks`. -/
theorem prob_lower_endpoint_attained (β : ℝ) (E : C → ℝ) (blk : C → ι) (ℓ u : ι → ℝ)
    (S A : Finset C) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    ∃ E' : C → ℝ, (∀ c ∈ S, E' c = E c) ∧
      (∀ c ∉ S, ℓ (blk c) ≤ E' c ∧ E' c ≤ u (blk c)) ∧
      prob β E' A =
        (ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) /
          ((ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) +
            (ZOn β E (Aᶜ ∩ S) + blockBound β blk ℓ (Aᶜ \ S))) := by
  have hlu : ∀ c, c ∉ S → ℓ (blk c) ≤ u (blk c) := fun c hc => (hℓ c hc).trans (hu c hc)
  obtain ⟨h1, h2⟩ := extremal_consistent E blk ℓ u u ℓ S A
    (fun c => by by_cases hc : c ∈ S <;> simp [hc, hlu c])
    (fun c => by by_cases hc : c ∈ S <;> simp [hc, hlu c])
  exact ⟨_, h1, h2, prob_extremal β E blk u ℓ S A⟩

/-- **Data-optimality of the upper endpoint.** -/
theorem prob_upper_endpoint_attained (β : ℝ) (E : C → ℝ) (blk : C → ι) (ℓ u : ι → ℝ)
    (S A : Finset C) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    ∃ E' : C → ℝ, (∀ c ∈ S, E' c = E c) ∧
      (∀ c ∉ S, ℓ (blk c) ≤ E' c ∧ E' c ≤ u (blk c)) ∧
      prob β E' A =
        (ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) /
          ((ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) +
            (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S))) := by
  have hlu : ∀ c, c ∉ S → ℓ (blk c) ≤ u (blk c) := fun c hc => (hℓ c hc).trans (hu c hc)
  obtain ⟨h1, h2⟩ := extremal_consistent E blk ℓ u ℓ u S A
    (fun c => by by_cases hc : c ∈ S <;> simp [hc, hlu c])
    (fun c => by by_cases hc : c ∈ S <;> simp [hc, hlu c])
  exact ⟨_, h1, h2, prob_extremal β E blk ℓ u S A⟩

end FoldBounds
