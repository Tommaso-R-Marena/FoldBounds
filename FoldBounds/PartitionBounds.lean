module

public import FoldBounds.Basic

/-!
# FoldBounds.PartitionBounds — certified enclosures of `Z` from partial enumeration

A branch-and-bound style computation of `Z` splits the conformation space into

* an *explored* set `S : Finset C`, whose Boltzmann weights are summed exactly, and
* the *unexplored* remainder `Sᶜ`, which is covered by *blocks* (`blk : C → ι`), each
  carrying a certified energy interval `ℓ i ≤ E c ≤ u i` for all `c` in block `i`.

The key quantity is the block bound
`blockBound β blk e T = ∑_{c ∈ T} exp (-β * e (blk c)) = ∑_i #{c ∈ T | blk c = i} · exp (-β e i)`,
which only needs block cardinalities and one number per block.  For `0 ≤ β`:

* `blockBound β blk ℓ T` is an upper bound for `Z_T` (energy *lower* bounds),
* `blockBound β blk u T` is a lower bound for `Z_T` (energy *upper* bounds).

We also prove the refinement monotonicity statements (finer blocks / sharper energy
bounds / more enumeration never loosen the enclosure) and an explicit width bound
showing that the enclosure becomes tight when block energy gaps shrink.
-/

@[expose] public section

open Finset Real

namespace FoldBounds

variable {C : Type*}

section Pointwise

variable {β : ℝ} (E : C → ℝ)

/-- An energy lower bound gives a weight upper bound (needs `0 ≤ β`). -/
theorem weight_le_of_le_energy (hβ : 0 ≤ β) {e : ℝ} {c : C} (h : e ≤ E c) :
    weight β E c ≤ Real.exp (-β * e) := by
  unfold weight
  exact Real.exp_le_exp.mpr (by nlinarith)

/-- An energy upper bound gives a weight lower bound (needs `0 ≤ β`). -/
theorem le_weight_of_energy_le (hβ : 0 ≤ β) {e : ℝ} {c : C} (h : E c ≤ e) :
    Real.exp (-β * e) ≤ weight β E c := by
  unfold weight
  exact Real.exp_le_exp.mpr (by nlinarith)

end Pointwise

section Uniform

variable {β : ℝ} (E : C → ℝ)

/-- Uniform tail bound: if every conformation in `T` has energy `≥ Emin`, then
`Z_T ≤ |T| · exp (-β Emin)`. -/
theorem ZOn_le_card_mul_exp (hβ : 0 ≤ β) (T : Finset C) (Emin : ℝ)
    (h : ∀ c ∈ T, Emin ≤ E c) : ZOn β E T ≤ T.card * Real.exp (-β * Emin) := by
  have := Finset.sum_le_card_nsmul T (weight β E) _
    (fun c hc => weight_le_of_le_energy E hβ (h c hc))
  simpa [ZOn, nsmul_eq_mul] using this

/-- Uniform tail bound: if every conformation in `T` has energy `≤ Emax`, then
`|T| · exp (-β Emax) ≤ Z_T`. -/
theorem card_mul_exp_le_ZOn (hβ : 0 ≤ β) (T : Finset C) (Emax : ℝ)
    (h : ∀ c ∈ T, E c ≤ Emax) : T.card * Real.exp (-β * Emax) ≤ ZOn β E T := by
  have := Finset.card_nsmul_le_sum T (weight β E) _
    (fun c hc => le_weight_of_energy_le E hβ (h c hc))
  simpa [ZOn, nsmul_eq_mul] using this

end Uniform

section Block

variable {ι : Type*}

/-- Block bound `∑_{c ∈ T} exp (-β e(blk c))`.  With `e` = energy lower bounds per block it
is an upper bound on `Z_T`; with `e` = energy upper bounds it is a lower bound. -/
noncomputable def blockBound (β : ℝ) (blk : C → ι) (e : ι → ℝ) (T : Finset C) : ℝ :=
  ∑ c ∈ T, Real.exp (-β * e (blk c))

/-- The block bound only needs block cardinalities: the "computable form". -/
theorem blockBound_eq_sum_card [Fintype ι] [DecidableEq ι] (β : ℝ) (blk : C → ι)
    (e : ι → ℝ) (T : Finset C) :
    blockBound β blk e T = ∑ i, ((T.filter (fun c => blk c = i)).card : ℝ) *
      Real.exp (-β * e i) := by
  unfold blockBound
  rw [← Finset.sum_fiberwise T blk]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.sum_congr rfl (g := fun _ => Real.exp (-β * e i)), Finset.sum_const,
    nsmul_eq_mul]
  intro c hc
  rw [(Finset.mem_filter.mp hc).2]

theorem blockBound_nonneg (β : ℝ) (blk : C → ι) (e : ι → ℝ) (T : Finset C) :
    0 ≤ blockBound β blk e T :=
  Finset.sum_nonneg fun _ _ => (Real.exp_pos _).le

variable {β : ℝ} (E : C → ℝ) (blk : C → ι)

/-- **Block upper bound.** Certified per-block energy lower bounds give an upper bound
on the partial partition function. -/
theorem ZOn_le_blockBound (hβ : 0 ≤ β) (ℓ : ι → ℝ) (T : Finset C)
    (hℓ : ∀ c ∈ T, ℓ (blk c) ≤ E c) : ZOn β E T ≤ blockBound β blk ℓ T :=
  Finset.sum_le_sum fun c hc => weight_le_of_le_energy E hβ (hℓ c hc)

/-- **Block lower bound.** Certified per-block energy upper bounds give a lower bound
on the partial partition function. -/
theorem blockBound_le_ZOn (hβ : 0 ≤ β) (u : ι → ℝ) (T : Finset C)
    (hu : ∀ c ∈ T, E c ≤ u (blk c)) : blockBound β blk u T ≤ ZOn β E T :=
  Finset.sum_le_sum fun c hc => le_weight_of_energy_le E hβ (hu c hc)

/-- If the energy is exactly constant on blocks (`E = e ∘ blk` on `T`), the block bound is
exact, for every real `β`.  This is the density-of-states / symmetry reduction formula
`Z_T = ∑_i #(block i) · exp (-β e_i)` (combine with `blockBound_eq_sum_card`). -/
theorem blockBound_eq_ZOn_of_exact (β : ℝ) (e : ι → ℝ) (T : Finset C)
    (he : ∀ c ∈ T, E c = e (blk c)) : blockBound β blk e T = ZOn β E T :=
  Finset.sum_congr rfl fun c hc => by rw [weight, he c hc]

/-- **Theorem 3 (generic form).** Any certified enclosure `L ≤ Z_{Sᶜ} ≤ U` of the unexplored
contribution yields a certified enclosure of `Z`. -/
theorem Z_mem_Icc_of_tail [Fintype C] [DecidableEq C] (β : ℝ) (S : Finset C)
    {L U : ℝ} (hL : L ≤ ZOn β E Sᶜ) (hU : ZOn β E Sᶜ ≤ U) :
    Z β E ∈ Set.Icc (ZOn β E S + L) (ZOn β E S + U) := by
  rw [← ZOn_add_ZOn_compl β E S]
  exact ⟨by linarith, by linarith⟩

/-- **Theorem 3 (block-certificate form).** Exact enumeration of `S` plus per-block energy
intervals on the unexplored conformations gives a certified enclosure of `Z`. -/
theorem Z_mem_Icc_of_blocks [Fintype C] [DecidableEq C] (hβ : 0 ≤ β) (S : Finset C)
    (ℓ u : ι → ℝ) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    Z β E ∈ Set.Icc (ZOn β E S + blockBound β blk u Sᶜ)
      (ZOn β E S + blockBound β blk ℓ Sᶜ) :=
  Z_mem_Icc_of_tail E β S
    (blockBound_le_ZOn E blk hβ u Sᶜ fun c hc => hu c (Finset.mem_compl.mp hc))
    (ZOn_le_blockBound E blk hβ ℓ Sᶜ fun c hc => hℓ c (Finset.mem_compl.mp hc))

/-- **Theorem 3 (uniform form).** With a single global energy interval
`Emin ≤ E ≤ Emax` on the unexplored conformations:
`Z_S + |Sᶜ| e^{-β Emax} ≤ Z ≤ Z_S + |Sᶜ| e^{-β Emin}`. -/
theorem Z_mem_Icc_of_uniform [Fintype C] [DecidableEq C] (hβ : 0 ≤ β) (S : Finset C)
    (Emin Emax : ℝ) (hmin : ∀ c ∉ S, Emin ≤ E c) (hmax : ∀ c ∉ S, E c ≤ Emax) :
    Z β E ∈ Set.Icc (ZOn β E S + Sᶜ.card * Real.exp (-β * Emax))
      (ZOn β E S + Sᶜ.card * Real.exp (-β * Emin)) :=
  Z_mem_Icc_of_tail E β S
    (card_mul_exp_le_ZOn E hβ Sᶜ Emax fun c hc => hmax c (Finset.mem_compl.mp hc))
    (ZOn_le_card_mul_exp E hβ Sᶜ Emin fun c hc => hmin c (Finset.mem_compl.mp hc))

/-! ### Refinement monotonicity (Theorem 5, partition-function part) -/

/-- Refining blocks / sharpening energy lower bounds (pointwise larger lower bound for
every conformation) can only decrease the upper bound. -/
theorem blockBound_antitone_refine {ι' : Type*} (hβ : 0 ≤ β) (blk' : C → ι')
    (ℓ : ι → ℝ) (ℓ' : ι' → ℝ) (T : Finset C) (h : ∀ c ∈ T, ℓ (blk c) ≤ ℓ' (blk' c)) :
    blockBound β blk' ℓ' T ≤ blockBound β blk ℓ T :=
  Finset.sum_le_sum fun c hc => Real.exp_le_exp.mpr (by nlinarith [h c hc])

/-- Refining a partition: `blk = proj ∘ blk'` and each child bound is at least its parent's.
Then the refined upper bound is no larger. -/
theorem blockBound_refine_partition {ι' : Type*} (hβ : 0 ≤ β) (blk' : C → ι')
    (proj : ι' → ι) (hproj : ∀ c, proj (blk' c) = blk c) (ℓ : ι → ℝ) (ℓ' : ι' → ℝ)
    (hℓ : ∀ j, ℓ (proj j) ≤ ℓ' j) (T : Finset C) :
    blockBound β blk' ℓ' T ≤ blockBound β blk ℓ T :=
  blockBound_antitone_refine blk hβ blk' ℓ ℓ' T fun c _ => hproj c ▸ hℓ (blk' c)

/-- Sharpening energy upper bounds can only increase the lower bound. -/
theorem blockBound_monotone_refine {ι' : Type*} (hβ : 0 ≤ β) (blk' : C → ι')
    (u : ι → ℝ) (u' : ι' → ℝ) (T : Finset C) (h : ∀ c ∈ T, u' (blk' c) ≤ u (blk c)) :
    blockBound β blk u T ≤ blockBound β blk' u' T :=
  Finset.sum_le_sum fun c hc => Real.exp_le_exp.mpr (by nlinarith [h c hc])


private theorem compl_sdiff_compl_eq [Fintype C] [DecidableEq C] (S S' : Finset C) :
    Sᶜ \ S'ᶜ = S' \ S := by
  ext c; simp [and_comm]

/-- Enumerating more conformations (`S ⊆ S'`) with the same valid block lower bounds never
increases the certified upper bound on `Z`. -/
theorem upper_antitone_enumeration [Fintype C] [DecidableEq C] (hβ : 0 ≤ β)
    (ℓ : ι → ℝ) {S S' : Finset C} (hSS' : S ⊆ S') (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) :
    ZOn β E S' + blockBound β blk ℓ S'ᶜ ≤ ZOn β E S + blockBound β blk ℓ Sᶜ := by
  have h1 : ZOn β E (S' \ S) + ZOn β E S = ZOn β E S' := Finset.sum_sdiff hSS'
  have hc : S'ᶜ ⊆ Sᶜ := Finset.compl_subset_compl.mpr hSS'
  have h2 : blockBound β blk ℓ (Sᶜ \ S'ᶜ) + blockBound β blk ℓ S'ᶜ =
      blockBound β blk ℓ Sᶜ := Finset.sum_sdiff hc
  rw [compl_sdiff_compl_eq] at h2
  have h3 : ZOn β E (S' \ S) ≤ blockBound β blk ℓ (S' \ S) :=
    ZOn_le_blockBound E blk hβ ℓ _ fun c hc => hℓ c (Finset.mem_sdiff.mp hc).2
  linarith

/-- Enumerating more conformations never decreases the certified lower bound on `Z`. -/
theorem lower_monotone_enumeration [Fintype C] [DecidableEq C] (hβ : 0 ≤ β)
    (u : ι → ℝ) {S S' : Finset C} (hSS' : S ⊆ S') (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    ZOn β E S + blockBound β blk u Sᶜ ≤ ZOn β E S' + blockBound β blk u S'ᶜ := by
  have h1 : ZOn β E (S' \ S) + ZOn β E S = ZOn β E S' := Finset.sum_sdiff hSS'
  have hc : S'ᶜ ⊆ Sᶜ := Finset.compl_subset_compl.mpr hSS'
  have h2 : blockBound β blk u (Sᶜ \ S'ᶜ) + blockBound β blk u S'ᶜ =
      blockBound β blk u Sᶜ := Finset.sum_sdiff hc
  rw [compl_sdiff_compl_eq] at h2
  have h3 : blockBound β blk u (S' \ S) ≤ ZOn β E (S' \ S) :=
    blockBound_le_ZOn E blk hβ u _ fun c hc => hu c (Finset.mem_sdiff.mp hc).2
  linarith

/-! ### Width of the block enclosure (convergence certificate) -/

/-- **Width bound.** For any real `β` and any block data `ℓ, u`:
`blockBound ℓ - blockBound u ≤ β · ∑_{c ∈ T} exp (-β ℓ(blk c)) · (u(blk c) - ℓ(blk c))`.
Hence shrinking block energy gaps forces the `Z` enclosure to collapse. -/
theorem blockBound_sub_le (β : ℝ) (ℓ u : ι → ℝ) (T : Finset C) :
    blockBound β blk ℓ T - blockBound β blk u T ≤
      β * ∑ c ∈ T, Real.exp (-β * ℓ (blk c)) * (u (blk c) - ℓ (blk c)) := by
  unfold blockBound
  rw [← Finset.sum_sub_distrib, Finset.mul_sum]
  refine Finset.sum_le_sum fun c _ => ?_
  have hsplit : Real.exp (-β * u (blk c)) =
      Real.exp (-β * ℓ (blk c)) * Real.exp (-(β * (u (blk c) - ℓ (blk c)))) := by
    rw [← Real.exp_add]; ring_nf
  have h1 := Real.add_one_le_exp (-(β * (u (blk c) - ℓ (blk c))))
  have hpos := Real.exp_pos (-β * ℓ (blk c))
  rw [hsplit]
  nlinarith

/-- **Relative width bound.** If every block has energy gap at most `δ` (and `0 ≤ β`),
the width of the block enclosure is at most `β δ` times the upper bound. -/
theorem blockBound_sub_le_mul (hβ : 0 ≤ β) (ℓ u : ι → ℝ) (T : Finset C) (δ : ℝ)
    (hgap : ∀ c ∈ T, u (blk c) - ℓ (blk c) ≤ δ) :
    blockBound β blk ℓ T - blockBound β blk u T ≤ β * δ * blockBound β blk ℓ T := by
  refine (blockBound_sub_le blk β ℓ u T).trans ?_
  rw [mul_assoc]
  refine mul_le_mul_of_nonneg_left ?_ hβ
  unfold blockBound
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun c hc => ?_
  have := Real.exp_pos (-β * ℓ (blk c))
  nlinarith [hgap c hc]

end Block

end FoldBounds
