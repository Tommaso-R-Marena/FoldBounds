# FoldBounds — exact theorem statements

Auto-extracted from the Lean sources (declaration text up to `:=`; section variables such as `{C : Type*} [Fintype C] [DecidableEq C] [Nonempty C]` declared with `variable` in the file are implicit arguments of each statement — see the file headers). Every declaration below builds without `sorry`; the headline theorems are audited in `FoldBounds/AxiomAudit.lean` and depend only on `propext`, `Classical.choice`, `Quot.sound`.

## `FoldBounds/Basic.lean`

Section variables: `variable {C : Type*}`; `variable [Fintype C]`; `variable (β : ℝ) (E : C → ℝ)`

### `FoldBounds.weight`
Boltzmann weight `exp (-β E(c))` of a conformation `c`.
```lean
noncomputable def weight (β : ℝ) (E : C → ℝ) (c : C) : ℝ
```

### `FoldBounds.ZOn`
Partial partition function `Z_A = ∑_{c ∈ A} exp (-β E(c))`.
```lean
noncomputable def ZOn (β : ℝ) (E : C → ℝ) (A : Finset C) : ℝ
```

### `FoldBounds.weight_pos`
```lean
theorem weight_pos (β : ℝ) (E : C → ℝ) (c : C) : 0 < weight β E c
```

### `FoldBounds.ZOn_empty`
```lean
@[simp] theorem ZOn_empty (β : ℝ) (E : C → ℝ) : ZOn β E ∅ = 0
```

### `FoldBounds.ZOn_union_of_disjoint`
```lean
theorem ZOn_union_of_disjoint [DecidableEq C] (β : ℝ) (E : C → ℝ) {A B : Finset C}
    (h : Disjoint A B) : ZOn β E (A ∪ B) = ZOn β E A + ZOn β E B
```

### `FoldBounds.ZOn_nonneg`
```lean
theorem ZOn_nonneg (β : ℝ) (E : C → ℝ) (A : Finset C) : 0 ≤ ZOn β E A
```

### `FoldBounds.ZOn_pos`
```lean
theorem ZOn_pos (β : ℝ) (E : C → ℝ) {A : Finset C} (hA : A.Nonempty) : 0 < ZOn β E A
```

### `FoldBounds.ZOn_mono`
```lean
theorem ZOn_mono (β : ℝ) (E : C → ℝ) {A B : Finset C} (h : A ⊆ B) :
    ZOn β E A ≤ ZOn β E B
```

### `FoldBounds.Z`
Partition function `Z = ∑_c exp (-β E(c))`.
```lean
noncomputable def Z (β : ℝ) (E : C → ℝ) : ℝ
```

### `FoldBounds.prob`
Equilibrium (Boltzmann) probability `P(A) = Z_A / Z`.
```lean
noncomputable def prob (β : ℝ) (E : C → ℝ) (A : Finset C) : ℝ
```

### `FoldBounds.ZOn_univ`
```lean
@[simp] theorem ZOn_univ : ZOn β E Finset.univ = Z β E
```

### `FoldBounds.Z_pos`
**Theorem 1.** The partition function of a nonempty finite conformation space is positive (for every real `β`).
```lean
theorem Z_pos [Nonempty C] : 0 < Z β E
```

### `FoldBounds.ZOn_le_Z`
```lean
theorem ZOn_le_Z (A : Finset C) : ZOn β E A ≤ Z β E
```

### `FoldBounds.ZOn_add_ZOn_compl`
`Z = Z_A + Z_{Aᶜ}`: the basic decomposition behind every bound in this library.
```lean
theorem ZOn_add_ZOn_compl [DecidableEq C] (A : Finset C) :
    ZOn β E A + ZOn β E Aᶜ = Z β E
```

### `FoldBounds.prob_nonneg`
```lean
theorem prob_nonneg (A : Finset C) : 0 ≤ prob β E A
```

### `FoldBounds.prob_le_one`
```lean
theorem prob_le_one (A : Finset C) : prob β E A ≤ 1
```

### `FoldBounds.prob_mem_Icc`
**Theorem 2.** `0 ≤ P(A) ≤ 1`.
```lean
theorem prob_mem_Icc (A : Finset C) : prob β E A ∈ Set.Icc (0 : ℝ) 1
```

### `FoldBounds.prob_univ`
```lean
theorem prob_univ [Nonempty C] : prob β E Finset.univ = 1
```

### `FoldBounds.prob_empty`
```lean
@[simp] theorem prob_empty : prob β E ∅ = 0
```

### `FoldBounds.prob_pos`
```lean
theorem prob_pos [Nonempty C] {A : Finset C} (hA : A.Nonempty) : 0 < prob β E A
```

### `FoldBounds.prob_union_of_disjoint`
```lean
theorem prob_union_of_disjoint [DecidableEq C] {A B : Finset C} (h : Disjoint A B) :
    prob β E (A ∪ B) = prob β E A + prob β E B
```

### `FoldBounds.prob_compl`
```lean
theorem prob_compl [DecidableEq C] [Nonempty C] (A : Finset C) :
    prob β E Aᶜ = 1 - prob β E A
```

### `FoldBounds.prob_eq_ratio`
`P(A) = Z_A / (Z_A + Z_{Aᶜ})`: the numerator also appears in the denominator.  This dependency is what the coupled probability bounds exploit.
```lean
theorem prob_eq_ratio [DecidableEq C] (A : Finset C) :
    prob β E A = ZOn β E A / (ZOn β E A + ZOn β E Aᶜ)
```

## `FoldBounds/PartitionBounds.lean`

Section variables: `variable {C : Type*}`; `variable {β : ℝ} (E : C → ℝ)`; `variable {β : ℝ} (E : C → ℝ)`; `variable {ι : Type*}`; `variable {β : ℝ} (E : C → ℝ) (blk : C → ι)`

### `FoldBounds.weight_le_of_le_energy`
An energy lower bound gives a weight upper bound (needs `0 ≤ β`).
```lean
theorem weight_le_of_le_energy (hβ : 0 ≤ β) {e : ℝ} {c : C} (h : e ≤ E c) :
    weight β E c ≤ Real.exp (-β * e)
```

### `FoldBounds.le_weight_of_energy_le`
An energy upper bound gives a weight lower bound (needs `0 ≤ β`).
```lean
theorem le_weight_of_energy_le (hβ : 0 ≤ β) {e : ℝ} {c : C} (h : E c ≤ e) :
    Real.exp (-β * e) ≤ weight β E c
```

### `FoldBounds.ZOn_le_card_mul_exp`
Uniform tail bound: if every conformation in `T` has energy `≥ Emin`, then `Z_T ≤ |T| · exp (-β Emin)`.
```lean
theorem ZOn_le_card_mul_exp (hβ : 0 ≤ β) (T : Finset C) (Emin : ℝ)
    (h : ∀ c ∈ T, Emin ≤ E c) : ZOn β E T ≤ T.card * Real.exp (-β * Emin)
```

### `FoldBounds.card_mul_exp_le_ZOn`
Uniform tail bound: if every conformation in `T` has energy `≤ Emax`, then `|T| · exp (-β Emax) ≤ Z_T`.
```lean
theorem card_mul_exp_le_ZOn (hβ : 0 ≤ β) (T : Finset C) (Emax : ℝ)
    (h : ∀ c ∈ T, E c ≤ Emax) : T.card * Real.exp (-β * Emax) ≤ ZOn β E T
```

### `FoldBounds.blockBound`
Block bound `∑_{c ∈ T} exp (-β e(blk c))`.  With `e` = energy lower bounds per block it is an upper bound on `Z_T`; with `e` = energy upper bounds it is a lower bound.
```lean
noncomputable def blockBound (β : ℝ) (blk : C → ι) (e : ι → ℝ) (T : Finset C) : ℝ
```

### `FoldBounds.blockBound_eq_sum_card`
The block bound only needs block cardinalities: the "computable form".
```lean
theorem blockBound_eq_sum_card [Fintype ι] [DecidableEq ι] (β : ℝ) (blk : C → ι)
    (e : ι → ℝ) (T : Finset C) :
    blockBound β blk e T = ∑ i, ((T.filter (fun c => blk c = i)).card : ℝ) *
      Real.exp (-β * e i)
```

### `FoldBounds.blockBound_nonneg`
```lean
theorem blockBound_nonneg (β : ℝ) (blk : C → ι) (e : ι → ℝ) (T : Finset C) :
    0 ≤ blockBound β blk e T
```

### `FoldBounds.ZOn_le_blockBound`
**Block upper bound.** Certified per-block energy lower bounds give an upper bound on the partial partition function.
```lean
theorem ZOn_le_blockBound (hβ : 0 ≤ β) (ℓ : ι → ℝ) (T : Finset C)
    (hℓ : ∀ c ∈ T, ℓ (blk c) ≤ E c) : ZOn β E T ≤ blockBound β blk ℓ T
```

### `FoldBounds.blockBound_le_ZOn`
**Block lower bound.** Certified per-block energy upper bounds give a lower bound on the partial partition function.
```lean
theorem blockBound_le_ZOn (hβ : 0 ≤ β) (u : ι → ℝ) (T : Finset C)
    (hu : ∀ c ∈ T, E c ≤ u (blk c)) : blockBound β blk u T ≤ ZOn β E T
```

### `FoldBounds.blockBound_eq_ZOn_of_exact`
If the energy is exactly constant on blocks (`E = e ∘ blk` on `T`), the block bound is exact, for every real `β`.  This is the density-of-states / symmetry reduction formula `Z_T = ∑_i #(block i) · exp (-β e_i)` (combine with `blockBound_eq_sum_card`).
```lean
theorem blockBound_eq_ZOn_of_exact (β : ℝ) (e : ι → ℝ) (T : Finset C)
    (he : ∀ c ∈ T, E c = e (blk c)) : blockBound β blk e T = ZOn β E T
```

### `FoldBounds.Z_mem_Icc_of_tail`
**Theorem 3 (generic form).** Any certified enclosure `L ≤ Z_{Sᶜ} ≤ U` of the unexplored contribution yields a certified enclosure of `Z`.
```lean
theorem Z_mem_Icc_of_tail [Fintype C] [DecidableEq C] (β : ℝ) (S : Finset C)
    {L U : ℝ} (hL : L ≤ ZOn β E Sᶜ) (hU : ZOn β E Sᶜ ≤ U) :
    Z β E ∈ Set.Icc (ZOn β E S + L) (ZOn β E S + U)
```

### `FoldBounds.Z_mem_Icc_of_blocks`
**Theorem 3 (block-certificate form).** Exact enumeration of `S` plus per-block energy intervals on the unexplored conformations gives a certified enclosure of `Z`.
```lean
theorem Z_mem_Icc_of_blocks [Fintype C] [DecidableEq C] (hβ : 0 ≤ β) (S : Finset C)
    (ℓ u : ι → ℝ) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    Z β E ∈ Set.Icc (ZOn β E S + blockBound β blk u Sᶜ)
      (ZOn β E S + blockBound β blk ℓ Sᶜ)
```

### `FoldBounds.Z_mem_Icc_of_uniform`
**Theorem 3 (uniform form).** With a single global energy interval `Emin ≤ E ≤ Emax` on the unexplored conformations: `Z_S + |Sᶜ| e^{-β Emax} ≤ Z ≤ Z_S + |Sᶜ| e^{-β Emin}`.
```lean
theorem Z_mem_Icc_of_uniform [Fintype C] [DecidableEq C] (hβ : 0 ≤ β) (S : Finset C)
    (Emin Emax : ℝ) (hmin : ∀ c ∉ S, Emin ≤ E c) (hmax : ∀ c ∉ S, E c ≤ Emax) :
    Z β E ∈ Set.Icc (ZOn β E S + Sᶜ.card * Real.exp (-β * Emax))
      (ZOn β E S + Sᶜ.card * Real.exp (-β * Emin))
```

### `FoldBounds.blockBound_antitone_refine`
Refining blocks / sharpening energy lower bounds (pointwise larger lower bound for every conformation) can only decrease the upper bound.
```lean
theorem blockBound_antitone_refine {ι' : Type*} (hβ : 0 ≤ β) (blk' : C → ι')
    (ℓ : ι → ℝ) (ℓ' : ι' → ℝ) (T : Finset C) (h : ∀ c ∈ T, ℓ (blk c) ≤ ℓ' (blk' c)) :
    blockBound β blk' ℓ' T ≤ blockBound β blk ℓ T
```

### `FoldBounds.blockBound_refine_partition`
Refining a partition: `blk = proj ∘ blk'` and each child bound is at least its parent's. Then the refined upper bound is no larger.
```lean
theorem blockBound_refine_partition {ι' : Type*} (hβ : 0 ≤ β) (blk' : C → ι')
    (proj : ι' → ι) (hproj : ∀ c, proj (blk' c) = blk c) (ℓ : ι → ℝ) (ℓ' : ι' → ℝ)
    (hℓ : ∀ j, ℓ (proj j) ≤ ℓ' j) (T : Finset C) :
    blockBound β blk' ℓ' T ≤ blockBound β blk ℓ T
```

### `FoldBounds.blockBound_monotone_refine`
Sharpening energy upper bounds can only increase the lower bound.
```lean
theorem blockBound_monotone_refine {ι' : Type*} (hβ : 0 ≤ β) (blk' : C → ι')
    (u : ι → ℝ) (u' : ι' → ℝ) (T : Finset C) (h : ∀ c ∈ T, u' (blk' c) ≤ u (blk c)) :
    blockBound β blk u T ≤ blockBound β blk' u' T
```

### `FoldBounds.upper_antitone_enumeration`
Enumerating more conformations (`S ⊆ S'`) with the same valid block lower bounds never increases the certified upper bound on `Z`.
```lean
theorem upper_antitone_enumeration [Fintype C] [DecidableEq C] (hβ : 0 ≤ β)
    (ℓ : ι → ℝ) {S S' : Finset C} (hSS' : S ⊆ S') (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) :
    ZOn β E S' + blockBound β blk ℓ S'ᶜ ≤ ZOn β E S + blockBound β blk ℓ Sᶜ
```

### `FoldBounds.lower_monotone_enumeration`
Enumerating more conformations never decreases the certified lower bound on `Z`.
```lean
theorem lower_monotone_enumeration [Fintype C] [DecidableEq C] (hβ : 0 ≤ β)
    (u : ι → ℝ) {S S' : Finset C} (hSS' : S ⊆ S') (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    ZOn β E S + blockBound β blk u Sᶜ ≤ ZOn β E S' + blockBound β blk u S'ᶜ
```

### `FoldBounds.blockBound_sub_le`
**Width bound.** For any real `β` and any block data `ℓ, u`: `blockBound ℓ - blockBound u ≤ β · ∑_{c ∈ T} exp (-β ℓ(blk c)) · (u(blk c) - ℓ(blk c))`. Hence shrinking block energy gaps forces the `Z` enclosure to collapse.
```lean
theorem blockBound_sub_le (β : ℝ) (ℓ u : ι → ℝ) (T : Finset C) :
    blockBound β blk ℓ T - blockBound β blk u T ≤
      β * ∑ c ∈ T, Real.exp (-β * ℓ (blk c)) * (u (blk c) - ℓ (blk c))
```

### `FoldBounds.blockBound_sub_le_mul`
**Relative width bound.** If every block has energy gap at most `δ` (and `0 ≤ β`), the width of the block enclosure is at most `β δ` times the upper bound.
```lean
theorem blockBound_sub_le_mul (hβ : 0 ≤ β) (ℓ u : ι → ℝ) (T : Finset C) (δ : ℝ)
    (hgap : ∀ c ∈ T, u (blk c) - ℓ (blk c) ≤ δ) :
    blockBound β blk ℓ T - blockBound β blk u T ≤ β * δ * blockBound β blk ℓ T
```

## `FoldBounds/Probability.lean`

Section variables: `variable {C : Type*} [Fintype C] [DecidableEq C] [Nonempty C] (β : ℝ) (E : C → ℝ)`; `variable {C : Type*} [Fintype C] [DecidableEq C] [Nonempty C] {ι : Type*}`

### `FoldBounds.ratio_le_ratio`
`a / (a + b)` is monotone in `a` and antitone in `b` on the nonnegative quadrant.
```lean
theorem ratio_le_ratio {a b a' b' : ℝ} (ha : 0 ≤ a) (hb' : 0 ≤ b') (haa : a ≤ a')
    (hbb : b' ≤ b) (h1 : 0 < a + b) (h2 : 0 < a' + b') :
    a / (a + b) ≤ a' / (a' + b')
```

### `FoldBounds.ratio_mem_coupled`
**Coupled enclosure (abstract form).** For nonnegative interval data and any `a, b` in the box with `a + b > 0`: `aL/(aL+bU) ≤ a/(a+b) ≤ aU/(aU+bL)`.
```lean
theorem ratio_mem_coupled {a b aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : a ∈ Set.Icc aL aU) (hb : b ∈ Set.Icc bL bU) (hab : 0 < a + b) :
    a / (a + b) ∈ Set.Icc (aL / (aL + bU)) (aU / (aU + bL))
```

### `FoldBounds.coupled_bounds_mem_unit`
The coupled enclosure always lies in `[0, 1]`.
```lean
theorem coupled_bounds_mem_unit {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (haU : 0 ≤ aU)
    (hbL : 0 ≤ bL) (hbU : 0 ≤ bU) :
    0 ≤ aL / (aL + bU) ∧ aU / (aU + bL) ≤ 1
```

### `FoldBounds.naive_lower_le_coupled_lower`
The coupled lower bound is at least the naive lower bound `aL / (aU + bU)`.
```lean
theorem naive_lower_le_coupled_lower {aL aU bU : ℝ} (haL : 0 ≤ aL) (hbU : 0 ≤ bU)
    (haLU : aL ≤ aU) : aL / (aU + bU) ≤ aL / (aL + bU)
```

### `FoldBounds.coupled_upper_le_naive_upper`
The coupled upper bound is at most the naive upper bound `aU / (aL + bL)`.
```lean
theorem coupled_upper_le_naive_upper {aL aU bL : ℝ} (haU : 0 ≤ aU) (haLU : aL ≤ aU)
    (hpos : 0 < aL + bL) : aU / (aU + bL) ≤ aU / (aL + bL)
```

### `FoldBounds.coupled_upper_lt_naive_upper`
Strict improvement: if the numerator interval has positive width and `aU > 0`, the coupled upper bound is strictly smaller than the naive one.
```lean
theorem coupled_upper_lt_naive_upper {aL aU bL : ℝ} (haU : 0 < aU) (haLU : aL < aU)
    (hpos : 0 < aL + bL) : aU / (aU + bL) < aU / (aL + bL)
```

### `FoldBounds.coupled_lower_isLeast`
**Optimality of the coupled lower bound**: it is the exact minimum of `a/(a+b)` over the box `[aL, aU] × [bL, bU]`.
```lean
theorem coupled_lower_isLeast {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (haLU : aL ≤ aU) (hbLU : bL ≤ bU) (hpos : 0 < aL + bL) :
    IsLeast {p | ∃ a ∈ Set.Icc aL aU, ∃ b ∈ Set.Icc bL bU, p = a / (a + b)}
      (aL / (aL + bU))
```

### `FoldBounds.coupled_upper_isGreatest`
**Optimality of the coupled upper bound**: it is the exact maximum of `a/(a+b)` over the box `[aL, aU] × [bL, bU]`.
```lean
theorem coupled_upper_isGreatest {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (haLU : aL ≤ aU) (hbLU : bL ≤ bU) (hpos : 0 < aL + bL) :
    IsGreatest {p | ∃ a ∈ Set.Icc aL aU, ∃ b ∈ Set.Icc bL bU, p = a / (a + b)}
      (aU / (aU + bL))
```

### `FoldBounds.coupled_width_le`
**Width of the coupled enclosure.** The probability enclosure is never wider than the relative width of the induced `Z` enclosure: `aU/(aU+bL) - aL/(aL+bU) ≤ ((aU - aL) + (bU - bL)) / (aL + bL)`. Combined with `blockBound_sub_le_mul` this turns block energy gaps into a guaranteed probability resolution (a stopping criterion for branch-and-bound).
```lean
theorem coupled_width_le {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL) (ha : aL ≤ aU)
    (hb : bL ≤ bU) (hpos : 0 < aL + bL) :
    aU / (aU + bL) - aL / (aL + bU) ≤ ((aU - aL) + (bU - bL)) / (aL + bL)
```

### `FoldBounds.coupled_enclosure_mono`
**Tightening intervals tightens the enclosure.** If `[aL', aU'] ⊆ [aL, aU]` and `[bL', bU'] ⊆ [bL, bU]` (nonnegative data, `aL + bL > 0`) then the coupled enclosure computed from the primed intervals is contained in the original one.
```lean
theorem coupled_enclosure_mono {aL aU bL bU aL' aU' bL' bU' : ℝ} (haL : 0 ≤ aL)
    (hbL : 0 ≤ bL) (hpos : 0 < aL + bL) (h1 : aL ≤ aL') (h2 : aU' ≤ aU) (h3 : bL ≤ bL')
    (h4 : bU' ≤ bU) (ha' : aL' ≤ aU') (hb' : bL' ≤ bU') :
    Set.Icc (aL' / (aL' + bU')) (aU' / (aU' + bL')) ⊆
      Set.Icc (aL / (aL + bU)) (aU / (aU + bL))
```

### `FoldBounds.mem_Icc_max_min`
Intersecting two valid enclosures of the same quantity is again a valid enclosure, contained in both: refinement by intersection can never invalidate a bound.
```lean
theorem mem_Icc_max_min {x l₁ u₁ l₂ u₂ : ℝ} (h₁ : x ∈ Set.Icc l₁ u₁)
    (h₂ : x ∈ Set.Icc l₂ u₂) :
    x ∈ Set.Icc (max l₁ l₂) (min u₁ u₂) ∧
      Set.Icc (max l₁ l₂) (min u₁ u₂) ⊆ Set.Icc l₁ u₁ ∩ Set.Icc l₂ u₂
```

### `FoldBounds.prob_mem_naive`
**Theorem 4 (naive enclosure).** Separate certified bounds on `Z_A` and on `Z`.
```lean
theorem prob_mem_naive (A : Finset C) {aL aU Zlo Zhi : ℝ} (haL : 0 ≤ aL) (hZlo : 0 < Zlo)
    (ha : ZOn β E A ∈ Set.Icc aL aU) (hZ : Z β E ∈ Set.Icc Zlo Zhi) :
    prob β E A ∈ Set.Icc (aL / Zhi) (aU / Zlo)
```

### `FoldBounds.prob_mem_coupled`
**Theorem 4 (coupled enclosure).** Certified bounds `aL ≤ Z_A ≤ aU`, `bL ≤ Z_{Aᶜ} ≤ bU` with `0 ≤ aL, 0 ≤ bL` give `aL/(aL+bU) ≤ P(A) ≤ aU/(aU+bL)`.
```lean
theorem prob_mem_coupled (A : Finset C) {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : ZOn β E A ∈ Set.Icc aL aU) (hb : ZOn β E Aᶜ ∈ Set.Icc bL bU) :
    prob β E A ∈ Set.Icc (aL / (aL + bU)) (aU / (aU + bL))
```

### `FoldBounds.prob_enclosure_refine`
**Theorem 5 (robust refinement).** Given two independently certified pairs of intervals for `(Z_A, Z_{Aᶜ})` (e.g. before and after a refinement step), the coupled enclosure built from their intersections is valid and contained in the coupled enclosure of the first pair.  In particular refinement never invalidates the enclosure.
```lean
theorem prob_enclosure_refine (A : Finset C) {aL aU bL bU aL' aU' bL' bU' : ℝ}
    (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : ZOn β E A ∈ Set.Icc aL aU) (hb : ZOn β E Aᶜ ∈ Set.Icc bL bU)
    (ha' : ZOn β E A ∈ Set.Icc aL' aU') (hb' : ZOn β E Aᶜ ∈ Set.Icc bL' bU') :
    prob β E A ∈ Set.Icc (max aL aL' / (max aL aL' + min bU bU'))
        (min aU aU' / (min aU aU' + max bL bL')) ∧
      Set.Icc (max aL aL' / (max aL aL' + min bU bU'))
          (min aU aU' / (min aU aU' + max bL bL')) ⊆
        Set.Icc (aL / (aL + bU)) (aU / (aU + bL))
```

### `FoldBounds.prob_mem_of_blocks`
**Theorem 4 (from partial enumeration).** Explored set `S`, block energy certificates `ℓ ≤ E ≤ u` on the unexplored conformations.  The four interval endpoints for `Z_A = Z_{A∩S} + Z_{A\S}` and `Z_{Aᶜ}` are computed from exact sums over explored conformations plus block bounds, and the coupled enclosure follows.
```lean
theorem prob_mem_of_blocks {ι : Type*} (hβ : 0 ≤ β) (blk : C → ι) (ℓ u : ι → ℝ)
    (S A : Finset C) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    prob β E A ∈ Set.Icc
      ((ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) /
        ((ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) +
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk ℓ (Aᶜ \ S))))
      ((ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) /
        ((ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) +
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S))))
```

### `FoldBounds.blocks_lower_sum_pos`
The certified lower end `aL + bL` of the block enclosure of `Z` is positive.
```lean
theorem blocks_lower_sum_pos (β : ℝ) (E : C → ℝ) (blk : C → ι) (u : ι → ℝ)
    (S A : Finset C) :
    0 < (ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) +
      (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S))
```

### `FoldBounds.prob_width_of_blocks`
**Convergence certificate for branch-and-bound.** If every unexplored block has energy gap at most `δ`, the width of the certified probability enclosure of `prob_mem_of_blocks` is at most `β δ · (unexplored upper mass) / (certified lower Z)`. Hence refining blocks until `δ → 0` (or enumerating everything) makes the enclosure collapse to the exact value.
```lean
theorem prob_width_of_blocks {β : ℝ} (hβ : 0 ≤ β) (E : C → ℝ) (blk : C → ι)
    (ℓ u : ι → ℝ) (S A : Finset C) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c)
    (hu : ∀ c ∉ S, E c ≤ u (blk c)) (δ : ℝ) (hgap : ∀ c ∉ S, u (blk c) - ℓ (blk c) ≤ δ) :
    (ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) /
        ((ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) +
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S))) -
      (ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) /
        ((ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) +
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk ℓ (Aᶜ \ S))) ≤
      β * δ * (blockBound β blk ℓ (A \ S) + blockBound β blk ℓ (Aᶜ \ S)) /
        ((ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) +
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S)))
```

## `FoldBounds/FreeEnergy.lean`

Section variables: `variable {C : Type*} [Fintype C]`; `variable {β : ℝ} (E : C → ℝ)`

### `FoldBounds.freeEnergy`
Helmholtz free energy `F = -(1/β) log Z`.
```lean
noncomputable def freeEnergy (β : ℝ) (E : C → ℝ) : ℝ
```

### `FoldBounds.freeEnergyOn`
Free energy of a macrostate `A`: `F_A = -(1/β) log Z_A`.
```lean
noncomputable def freeEnergyOn (β : ℝ) (E : C → ℝ) (A : Finset C) : ℝ
```

### `FoldBounds.neg_log_div_mem_Icc`
Monotone transfer of an enclosure of a positive quantity through `x ↦ -(log x)/β`.
```lean
theorem neg_log_div_mem_Icc (hβ : 0 < β) {x lo hi : ℝ} (hlo : 0 < lo)
    (hx : x ∈ Set.Icc lo hi) :
    -Real.log x / β ∈ Set.Icc (-Real.log hi / β) (-Real.log lo / β)
```

### `FoldBounds.freeEnergy_mem_Icc`
**Theorem 6.** A certified enclosure `0 < Zlo ≤ Z ≤ Zhi` gives the certified enclosure `-(log Zhi)/β ≤ F ≤ -(log Zlo)/β` for `0 < β`.
```lean
theorem freeEnergy_mem_Icc (hβ : 0 < β) {Zlo Zhi : ℝ} (hZlo : 0 < Zlo)
    (hZ : Z β E ∈ Set.Icc Zlo Zhi) :
    freeEnergy β E ∈ Set.Icc (-Real.log Zhi / β) (-Real.log Zlo / β)
```

### `FoldBounds.freeEnergy_le_energy`
The free energy never exceeds the energy of any conformation (`F ≤ E_min`).
```lean
theorem freeEnergy_le_energy [Nonempty C] (hβ : 0 < β) (c : C) :
    freeEnergy β E ≤ E c
```

### `FoldBounds.le_freeEnergy_of_energy_ge`
Lower bound by the ground-state energy and the entropy of the full space: `Emin - (log |C|)/β ≤ F` whenever `Emin ≤ E c` for all `c`.
```lean
theorem le_freeEnergy_of_energy_ge [Nonempty C] (hβ : 0 < β) (Emin : ℝ)
    (h : ∀ c, Emin ≤ E c) :
    Emin - Real.log (Fintype.card C) / β ≤ freeEnergy β E
```

### `FoldBounds.freeEnergy_le_freeEnergyOn`
Variational (enumeration) bound: for any nonempty explored set `S`, `F ≤ F_S`.
```lean
theorem freeEnergy_le_freeEnergyOn (hβ : 0 < β) {S : Finset C} (hS : S.Nonempty) :
    freeEnergy β E ≤ freeEnergyOn β E S
```

### `FoldBounds.prob_eq_exp_freeEnergy`
Free-energy representation of probabilities: `P(A) = exp (-β (F_A - F))`.
```lean
theorem prob_eq_exp_freeEnergy [Nonempty C] (hβ : 0 < β) {A : Finset C}
    (hA : A.Nonempty) :
    prob β E A = Real.exp (-β * (freeEnergyOn β E A - freeEnergy β E))
```

### `FoldBounds.freeEnergy_diff_mem_Icc`
**Certified free-energy differences.** A certified probability enclosure `0 < pL ≤ P(A) ≤ pU` yields `-(log pU)/β ≤ F_A - F ≤ -(log pL)/β`.
```lean
theorem freeEnergy_diff_mem_Icc [Nonempty C] (hβ : 0 < β) {A : Finset C}
    (hA : A.Nonempty) {pL pU : ℝ} (hpL : 0 < pL) (hp : prob β E A ∈ Set.Icc pL pU) :
    freeEnergyOn β E A - freeEnergy β E ∈
      Set.Icc (-Real.log pU / β) (-Real.log pL / β)
```

### `FoldBounds.certifiedZLower_pos`
The certified lower end of the block enclosure of `Z` is strictly positive.
```lean
theorem certifiedZLower_pos [DecidableEq C] [Nonempty C] {ι : Type*} (blk : C → ι)
    (u : ι → ℝ) (S : Finset C) :
    0 < ZOn β E S + blockBound β blk u Sᶜ
```

### `FoldBounds.freeEnergy_mem_Icc_of_blocks`
**Theorem 6 (from partial enumeration).** Exact enumeration of `S` and block energy certificates on `Sᶜ` give a certified free-energy enclosure, for `0 < β`.
```lean
theorem freeEnergy_mem_Icc_of_blocks [DecidableEq C] [Nonempty C] {ι : Type*}
    (hβ : 0 < β) (blk : C → ι) (ℓ u : ι → ℝ) (S : Finset C)
    (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    freeEnergy β E ∈
      Set.Icc (-Real.log (ZOn β E S + blockBound β blk ℓ Sᶜ) / β)
        (-Real.log (ZOn β E S + blockBound β blk u Sᶜ) / β)
```

## `FoldBounds/Symmetry.lean`

Section variables: `variable {C : Type*} [Fintype C] {D : Type*} [Fintype D] [DecidableEq D]`; `variable (β : ℝ) (E : C → ℝ)`; `variable {G : Type*} [Group G] [MulAction G C]`; `variable [Fintype (orbitRel.Quotient G C)] [DecidableEq (orbitRel.Quotient G C)]`

### `FoldBounds.multiplicity`
Multiplicity of a reduced state `d`: the number of conformations mapped to `d`.
```lean
def multiplicity (f : C → D) (d : D) : ℕ
```

### `FoldBounds.Z_eq_sum_multiplicity`
**Counting formula.** If `E = e ∘ f` then `Z = ∑_d m_d · exp (-β e_d)` (any real `β`).
```lean
theorem Z_eq_sum_multiplicity (f : C → D) (e : D → ℝ) (hE : ∀ c, E c = e (f c)) :
    Z β E = ∑ d, (multiplicity f d : ℝ) * Real.exp (-β * e d)
```

### `FoldBounds.Z_mem_Icc_of_reduced`
**Reduced certified enclosure.** Explored reduced states `R ⊆ D` are evaluated exactly; every unexplored reduced state carries an energy interval `ℓ d ≤ e d ≤ u d`.  Then `∑_{d∈R} m_d e^{-β e_d} + ∑_{d∉R} m_d e^{-β u_d} ≤ Z ≤ ∑_{d∈R} m_d e^{-β e_d} + ∑_{d∉R} m_d e^{-β ℓ_d}`.
```lean
theorem Z_mem_Icc_of_reduced (hβ : 0 ≤ β) (f : C → D) (e : D → ℝ)
    (hE : ∀ c, E c = e (f c)) (R : Finset D) (ℓ u : D → ℝ)
    (hℓ : ∀ d ∉ R, ℓ d ≤ e d) (hu : ∀ d ∉ R, e d ≤ u d) :
    Z β E ∈ Set.Icc
      (∑ d ∈ R, (multiplicity f d : ℝ) * Real.exp (-β * e d) +
        ∑ d ∈ Rᶜ, (multiplicity f d : ℝ) * Real.exp (-β * u d))
      (∑ d ∈ R, (multiplicity f d : ℝ) * Real.exp (-β * e d) +
        ∑ d ∈ Rᶜ, (multiplicity f d : ℝ) * Real.exp (-β * ℓ d))
```

### `FoldBounds.energy_eq_of_mem_orbit`
Energy invariance under a group action makes the energy constant on orbits.
```lean
theorem energy_eq_of_mem_orbit (hinv : ∀ (g : G) c, E (g • c) = E c) {c c' : C}
    (h : c ∈ orbit G c') : E c = E c'
```

### `FoldBounds.energy_eq_energy_out`
The energy of a conformation equals the energy of its orbit's chosen representative.
```lean
theorem energy_eq_energy_out (hinv : ∀ (g : G) c, E (g • c) = E c) (c : C) :
    E c = E (Quotient.out (Quotient.mk (orbitRel G C) c))
```

### `FoldBounds.Z_eq_sum_orbits`
**Symmetry reduction.** For an energy-preserving action of `G` on `C`, `Z = ∑_{orbits ω} |ω| · exp (-β E(rep ω))`, where `|ω|` is the number of conformations in the orbit and `rep ω` any chosen representative.
```lean
theorem Z_eq_sum_orbits (hinv : ∀ (g : G) c, E (g • c) = E c) :
    Z β E = ∑ ω : orbitRel.Quotient G C,
      (multiplicity (Quotient.mk (orbitRel G C)) ω : ℝ) * Real.exp (-β * E (Quotient.out ω))
```

### `FoldBounds.multiplicity_orbit_eq`
The multiplicity of an orbit in the reduction formula is its size.
```lean
theorem multiplicity_orbit_eq (ω : orbitRel.Quotient G C) :
    multiplicity (Quotient.mk (orbitRel G C)) ω =
      (Finset.univ.filter fun c => c ∈ orbit G (Quotient.out ω)).card
```

## `FoldBounds/Compositional.lean`

Section variables: `variable {C₁ C₂ : Type*} [Fintype C₁] [Fintype C₂]`; `variable (β : ℝ) (E₁ : C₁ → ℝ) (E₂ : C₂ → ℝ) (J : C₁ → C₂ → ℝ)`

### `FoldBounds.coupledEnergy`
Two-domain energy `E(x, y) = E₁ x + E₂ y + J x y`.
```lean
def coupledEnergy (E₁ : C₁ → ℝ) (E₂ : C₂ → ℝ) (J : C₁ → C₂ → ℝ) : C₁ × C₂ → ℝ
```

### `FoldBounds.Z_coupledEnergy`
```lean
theorem Z_coupledEnergy :
    Z β (coupledEnergy E₁ E₂ J) =
      ∑ x, ∑ y, weight β E₁ x * weight β E₂ y * Real.exp (-β * J x y)
```

### `FoldBounds.Z_coupledEnergy_zero`
**Exact factorisation** for non-interacting domains (any real `β`).
```lean
theorem Z_coupledEnergy_zero :
    Z β (coupledEnergy E₁ E₂ fun _ _ => 0) = Z β E₁ * Z β E₂
```

### `FoldBounds.Z_coupledEnergy_mem_Icc`
**Compositional partition-function bound.**
```lean
theorem Z_coupledEnergy_mem_Icc {β : ℝ} (hβ : 0 ≤ β) {Jlo Jhi : ℝ}
    (hJ : ∀ x y, J x y ∈ Set.Icc Jlo Jhi) :
    Z β (coupledEnergy E₁ E₂ J) ∈
      Set.Icc (Z β E₁ * Z β E₂ * Real.exp (-β * Jhi))
        (Z β E₁ * Z β E₂ * Real.exp (-β * Jlo))
```

### `FoldBounds.freeEnergy_coupledEnergy_mem_Icc`
**Compositional free-energy bound**: `F₁ + F₂ + Jlo ≤ F ≤ F₁ + F₂ + Jhi` for `0 < β`.
```lean
theorem freeEnergy_coupledEnergy_mem_Icc [Nonempty C₁] [Nonempty C₂] {β : ℝ} (hβ : 0 < β)
    {Jlo Jhi : ℝ} (hJ : ∀ x y, J x y ∈ Set.Icc Jlo Jhi) :
    freeEnergy β (coupledEnergy E₁ E₂ J) ∈
      Set.Icc (freeEnergy β E₁ + freeEnergy β E₂ + Jlo)
        (freeEnergy β E₁ + freeEnergy β E₂ + Jhi)
```

## `FoldBounds/Variational.lean`

Section variables: `variable {C : Type*} [Fintype C]`

### `FoldBounds.refAverage`
Reference-ensemble average `⟨V⟩₀ = ∑ w₀ V / Z₀`.
```lean
noncomputable def refAverage (β : ℝ) (E₀ V : C → ℝ) : ℝ
```

### `FoldBounds.Z_ge_gibbs_bogoliubov`
**Gibbs–Bogoliubov lower bound on `Z`** (any real `β`).
```lean
theorem Z_ge_gibbs_bogoliubov [Nonempty C] (β : ℝ) (E₀ V : C → ℝ) :
    Z β E₀ * Real.exp (-β * refAverage β E₀ V) ≤ Z β (fun c => E₀ c + V c)
```

### `FoldBounds.freeEnergy_le_gibbs_bogoliubov`
**Gibbs–Bogoliubov upper bound on `F`**: `F(E₀ + V) ≤ F(E₀) + ⟨V⟩₀` for `0 < β`.
```lean
theorem freeEnergy_le_gibbs_bogoliubov [Nonempty C] {β : ℝ} (hβ : 0 < β) (E₀ V : C → ℝ) :
    freeEnergy β (fun c => E₀ c + V c) ≤ freeEnergy β E₀ + refAverage β E₀ V
```

### `FoldBounds.refAverage_le`
The reference average never exceeds a uniform upper bound on the perturbation, so the Gibbs–Bogoliubov bound is at least as strong as the worst-case bound.
```lean
theorem refAverage_le [Nonempty C] (β : ℝ) (E₀ V : C → ℝ) {Vhi : ℝ} (h : ∀ c, V c ≤ Vhi) :
    refAverage β E₀ V ≤ Vhi
```

### `FoldBounds.Z_coupledEnergy_ge_gibbs`
**Compositional Gibbs–Bogoliubov bound.** With the decoupled two-domain landscape as reference, `Z ≥ Z₁ Z₂ · exp (-β ⟨J⟩₁₂)` where `⟨J⟩₁₂` is the average coupling in the product ensemble.  By `refAverage_le` this is never weaker than `Z ≥ Z₁ Z₂ e^{-β Jhi}`.
```lean
theorem Z_coupledEnergy_ge_gibbs {C₁ C₂ : Type*} [Fintype C₁] [Fintype C₂] [Nonempty C₁]
    [Nonempty C₂] (β : ℝ) (E₁ : C₁ → ℝ) (E₂ : C₂ → ℝ) (J : C₁ → C₂ → ℝ) :
    Z β E₁ * Z β E₂ *
        Real.exp (-β * refAverage β (coupledEnergy E₁ E₂ fun _ _ => 0) (fun p => J p.1 p.2)) ≤
      Z β (coupledEnergy E₁ E₂ J)
```

## `FoldBounds/Optimality.lean`

Section variables: `variable {C : Type*} [Fintype C] [DecidableEq C] {ι : Type*}`

### `FoldBounds.extremalLandscape`
Landscape that agrees with `E` on `S` and puts every unexplored conformation of `A` at energy `eA (blk c)` and every unexplored conformation of `Aᶜ` at `eB (blk c)`.
```lean
noncomputable def extremalLandscape (E : C → ℝ) (blk : C → ι) (eA eB : ι → ℝ)
    (S A : Finset C) : C → ℝ
```

### `FoldBounds.prob_lower_endpoint_attained`
**Data-optimality of the lower endpoint.** For any valid block certificate there is an energy landscape `E'` with the *same* certificate data (equal to `E` on `S`, within `[ℓ, u]` blockwise off `S`) whose probability `P_{E'}(A)` equals the certified lower endpoint of `prob_mem_of_blocks`.
```lean
theorem prob_lower_endpoint_attained (β : ℝ) (E : C → ℝ) (blk : C → ι) (ℓ u : ι → ℝ)
    (S A : Finset C) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    ∃ E' : C → ℝ, (∀ c ∈ S, E' c = E c) ∧
      (∀ c ∉ S, ℓ (blk c) ≤ E' c ∧ E' c ≤ u (blk c)) ∧
      prob β E' A =
        (ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) /
          ((ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) +
            (ZOn β E (Aᶜ ∩ S) + blockBound β blk ℓ (Aᶜ \ S)))
```

### `FoldBounds.prob_upper_endpoint_attained`
**Data-optimality of the upper endpoint.**
```lean
theorem prob_upper_endpoint_attained (β : ℝ) (E : C → ℝ) (blk : C → ι) (ℓ u : ι → ℝ)
    (S A : Finset C) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    ∃ E' : C → ℝ, (∀ c ∈ S, E' c = E c) ∧
      (∀ c ∉ S, ℓ (blk c) ≤ E' c ∧ E' c ≤ u (blk c)) ∧
      prob β E' A =
        (ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) /
          ((ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) +
            (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S)))
```

## `FoldBounds/Certificate.lean`

Section variables: `variable {C : Type*} [Fintype C] (β : ℝ) (E : C → ℝ)`

### `FoldBounds.coupledLowerQ`
Coupled lower bound computed in exact rational arithmetic.
```lean
def coupledLowerQ (aL bU : ℚ) : ℚ
```

### `FoldBounds.coupledUpperQ`
Coupled upper bound computed in exact rational arithmetic.
```lean
def coupledUpperQ (aU bL : ℚ) : ℚ
```

### `FoldBounds.ZOn_mem_Icc_of_rat`
Rational weight enclosures sum (exactly, in `ℚ`) to an enclosure of `Z_T`.
```lean
theorem ZOn_mem_Icc_of_rat (T : Finset C) (wlo whi : C → ℚ)
    (hw : ∀ c ∈ T, ((wlo c : ℚ) : ℝ) ≤ weight β E c ∧ weight β E c ≤ ((whi c : ℚ) : ℝ)) :
    ZOn β E T ∈ Set.Icc (((∑ c ∈ T, wlo c : ℚ)) : ℝ) (((∑ c ∈ T, whi c : ℚ)) : ℝ)
```

### `FoldBounds.prob_mem_Icc_of_rat`
**Soundness of the rational probability certificate.** If rational numbers `aL ≤ Z_A ≤ aU`, `bL ≤ Z_{Aᶜ} ≤ bU` are certified (as real inequalities), then the rationally computed `coupledLowerQ aL bU`, `coupledUpperQ aU bL` enclose `P(A)`.
```lean
theorem prob_mem_Icc_of_rat [DecidableEq C] [Nonempty C] (A : Finset C)
    {aL aU bL bU : ℚ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : ZOn β E A ∈ Set.Icc (aL : ℝ) aU) (hb : ZOn β E Aᶜ ∈ Set.Icc (bL : ℝ) bU) :
    prob β E A ∈ Set.Icc ((coupledLowerQ aL bU : ℚ) : ℝ) ((coupledUpperQ aU bL : ℚ) : ℝ)
```

### `FoldBounds.toyEnergy`
Toy energy landscape.
```lean
noncomputable def toyEnergy : Fin 3 → ℝ
```

### `FoldBounds.exp_neg_one_mem`
```lean
theorem exp_neg_one_mem : Real.exp (-1) ∈ Set.Icc (0.3678 : ℝ) 0.3679
```

### `FoldBounds.exp_neg_two_mem`
```lean
theorem exp_neg_two_mem : Real.exp (-2) ∈ Set.Icc (0.1352 : ℝ) 0.1354
```

### `FoldBounds.toy_prob_mem`
**Certified numerical result**: `0.6652 ≤ P({0}) ≤ 0.6654` for the toy landscape, obtained from `prob_mem_Icc_of_rat` with exact rational arithmetic.
```lean
theorem toy_prob_mem : prob 1 toyEnergy {0} ∈ Set.Icc (0.6652 : ℝ) 0.6654
```
