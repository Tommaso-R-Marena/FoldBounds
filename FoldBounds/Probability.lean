module

public import FoldBounds.PartitionBounds

/-!
# FoldBounds.Probability — certified enclosures of `P(A)`

Write `a = Z_A` and `b = Z_{Aᶜ}`; then `P(A) = a / (a + b)` (`prob_eq_ratio`).
Given certified intervals `a ∈ [aL, aU]`, `b ∈ [bL, bU]` we prove:

* **naive enclosure** `aL / Zhi ≤ P(A) ≤ aU / Zlo` from separate bounds on `Z_A` and `Z`;
* **coupled enclosure** `aL / (aL + bU) ≤ P(A) ≤ aU / (aU + bL)`, which uses that the
  numerator also occurs in the denominator;
* the coupled enclosure is never worse than the naive one built from the same data
  (`Zlo = aL + bL`, `Zhi = aU + bU`), is strictly better on the upper side whenever the
  numerator interval has positive width, and always lies inside `[0, 1]`;
* the coupled enclosure is **optimal**: its endpoints are the exact minimum and maximum
  of `a / (a + b)` over the box, so no method that only uses these two intervals can do
  better;
* refinement statements: shrinking intervals shrinks the enclosure, and intersecting
  valid intervals never invalidates it.
-/

@[expose] public section

open Finset Real

namespace FoldBounds

/-! ### The ratio `a / (a + b)` -/

/-- `a / (a + b)` is monotone in `a` and antitone in `b` on the nonnegative quadrant. -/
theorem ratio_le_ratio {a b a' b' : ℝ} (ha : 0 ≤ a) (hb' : 0 ≤ b') (haa : a ≤ a')
    (hbb : b' ≤ b) (h1 : 0 < a + b) (h2 : 0 < a' + b') :
    a / (a + b) ≤ a' / (a' + b') := by
  rw [div_le_div_iff₀ h1 h2]
  nlinarith [mul_le_mul haa hbb hb' (ha.trans haa)]

/-- **Coupled enclosure (abstract form).** For nonnegative interval data and any `a, b`
in the box with `a + b > 0`: `aL/(aL+bU) ≤ a/(a+b) ≤ aU/(aU+bL)`. -/
theorem ratio_mem_coupled {a b aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : a ∈ Set.Icc aL aU) (hb : b ∈ Set.Icc bL bU) (hab : 0 < a + b) :
    a / (a + b) ∈ Set.Icc (aL / (aL + bU)) (aU / (aU + bL)) := by
  obtain ⟨ha1, ha2⟩ := ha
  obtain ⟨hb1, hb2⟩ := hb
  constructor
  · rcases (add_nonneg haL (hbL.trans (hb1.trans hb2))).eq_or_lt with h | h
    · rw [← h, div_zero]; exact div_nonneg (haL.trans ha1) hab.le
    · exact ratio_le_ratio haL (hbL.trans hb1) ha1 hb2 h hab
  · rcases (haL.trans ha1).eq_or_lt with h | h
    · rw [← h, zero_div]; exact div_nonneg (h ▸ ha2) (add_nonneg (h ▸ ha2) hbL)
    · exact ratio_le_ratio (haL.trans ha1) hbL ha2 hb1 hab (by linarith)

/-- The coupled enclosure always lies in `[0, 1]`. -/
theorem coupled_bounds_mem_unit {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (haU : 0 ≤ aU)
    (hbL : 0 ≤ bL) (hbU : 0 ≤ bU) :
    0 ≤ aL / (aL + bU) ∧ aU / (aU + bL) ≤ 1 :=
  ⟨div_nonneg haL (add_nonneg haL hbU),
    div_le_one_of_le₀ (le_add_of_nonneg_right hbL) (add_nonneg haU hbL)⟩

/-- The coupled lower bound is at least the naive lower bound `aL / (aU + bU)`. -/
theorem naive_lower_le_coupled_lower {aL aU bU : ℝ} (haL : 0 ≤ aL) (hbU : 0 ≤ bU)
    (haLU : aL ≤ aU) : aL / (aU + bU) ≤ aL / (aL + bU) := by
  rcases (add_nonneg haL hbU).eq_or_lt with h | h
  · have : aL = 0 := by linarith
    simp [this]
  · exact div_le_div_of_nonneg_left haL h (by linarith)

/-- The coupled upper bound is at most the naive upper bound `aU / (aL + bL)`. -/
theorem coupled_upper_le_naive_upper {aL aU bL : ℝ} (haU : 0 ≤ aU) (haLU : aL ≤ aU)
    (hpos : 0 < aL + bL) : aU / (aU + bL) ≤ aU / (aL + bL) :=
  div_le_div_of_nonneg_left haU hpos (by linarith)

/-- Strict improvement: if the numerator interval has positive width and `aU > 0`, the
coupled upper bound is strictly smaller than the naive one. -/
theorem coupled_upper_lt_naive_upper {aL aU bL : ℝ} (haU : 0 < aU) (haLU : aL < aU)
    (hpos : 0 < aL + bL) : aU / (aU + bL) < aU / (aL + bL) :=
  div_lt_div_of_pos_left haU hpos (by linarith)

/-- **Optimality of the coupled lower bound**: it is the exact minimum of `a/(a+b)` over
the box `[aL, aU] × [bL, bU]`. -/
theorem coupled_lower_isLeast {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (haLU : aL ≤ aU) (hbLU : bL ≤ bU) (hpos : 0 < aL + bL) :
    IsLeast {p | ∃ a ∈ Set.Icc aL aU, ∃ b ∈ Set.Icc bL bU, p = a / (a + b)}
      (aL / (aL + bU)) := by
  refine ⟨⟨aL, ⟨le_rfl, haLU⟩, bU, ⟨hbLU, le_rfl⟩, rfl⟩, ?_⟩
  rintro p ⟨a, ha, b, hb, rfl⟩
  exact (ratio_mem_coupled haL hbL ha hb (by linarith [ha.1, hb.1])).1

/-- **Optimality of the coupled upper bound**: it is the exact maximum of `a/(a+b)` over
the box `[aL, aU] × [bL, bU]`. -/
theorem coupled_upper_isGreatest {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (haLU : aL ≤ aU) (hbLU : bL ≤ bU) (hpos : 0 < aL + bL) :
    IsGreatest {p | ∃ a ∈ Set.Icc aL aU, ∃ b ∈ Set.Icc bL bU, p = a / (a + b)}
      (aU / (aU + bL)) := by
  refine ⟨⟨aU, ⟨haLU, le_rfl⟩, bL, ⟨le_rfl, hbLU⟩, rfl⟩, ?_⟩
  rintro p ⟨a, ha, b, hb, rfl⟩
  exact (ratio_mem_coupled haL hbL ha hb (by linarith [ha.1, hb.1])).2

/-- **Width of the coupled enclosure.** The probability enclosure is never wider than
the relative width of the induced `Z` enclosure:
`aU/(aU+bL) - aL/(aL+bU) ≤ ((aU - aL) + (bU - bL)) / (aL + bL)`.
Combined with `blockBound_sub_le_mul` this turns block energy gaps into a guaranteed
probability resolution (a stopping criterion for branch-and-bound). -/
theorem coupled_width_le {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL) (ha : aL ≤ aU)
    (hb : bL ≤ bU) (hpos : 0 < aL + bL) :
    aU / (aU + bL) - aL / (aL + bU) ≤ ((aU - aL) + (bU - bL)) / (aL + bL) := by
  have h1 : 0 < aU + bL := by linarith
  have h2 : 0 < aL + bU := by linarith
  rw [div_sub_div _ _ h1.ne' h2.ne', div_le_div_iff₀ (mul_pos h1 h2) hpos]
  have hq : 0 ≤ aU - aL := by linarith
  have hs : 0 ≤ bU - bL := by linarith
  nlinarith [mul_nonneg (mul_nonneg hq hs) (add_nonneg haL hbL),
    mul_nonneg (mul_nonneg hq hq) (add_nonneg haL hbL),
    mul_nonneg (mul_nonneg hs hs) (add_nonneg haL hbL), mul_nonneg (mul_nonneg hq hs) hq,
    mul_nonneg (mul_nonneg hq hs) hs, mul_nonneg (mul_nonneg haL hbL) hq,
    mul_nonneg (mul_nonneg haL hbL) hs, mul_nonneg (mul_nonneg haL haL) hq,
    mul_nonneg (mul_nonneg hbL hbL) hs, mul_nonneg (mul_nonneg haL haL) hs,
    mul_nonneg (mul_nonneg hbL hbL) hq]

/-! ### Refinement (Theorem 5, probability part) -/

/-- **Tightening intervals tightens the enclosure.** If `[aL', aU'] ⊆ [aL, aU]` and
`[bL', bU'] ⊆ [bL, bU]` (nonnegative data, `aL + bL > 0`) then the coupled enclosure
computed from the primed intervals is contained in the original one. -/
theorem coupled_enclosure_mono {aL aU bL bU aL' aU' bL' bU' : ℝ} (haL : 0 ≤ aL)
    (hbL : 0 ≤ bL) (hpos : 0 < aL + bL) (h1 : aL ≤ aL') (h2 : aU' ≤ aU) (h3 : bL ≤ bL')
    (h4 : bU' ≤ bU) (ha' : aL' ≤ aU') (hb' : bL' ≤ bU') :
    Set.Icc (aL' / (aL' + bU')) (aU' / (aU' + bL')) ⊆
      Set.Icc (aL / (aL + bU)) (aU / (aU + bL)) := by
  apply Set.Icc_subset_Icc
  · exact ratio_le_ratio haL (by linarith) h1 h4 (by linarith) (by linarith)
  · exact ratio_le_ratio (by linarith) hbL h2 h3 (by linarith) (by linarith)

/-- Intersecting two valid enclosures of the same quantity is again a valid enclosure,
contained in both: refinement by intersection can never invalidate a bound. -/
theorem mem_Icc_max_min {x l₁ u₁ l₂ u₂ : ℝ} (h₁ : x ∈ Set.Icc l₁ u₁)
    (h₂ : x ∈ Set.Icc l₂ u₂) :
    x ∈ Set.Icc (max l₁ l₂) (min u₁ u₂) ∧
      Set.Icc (max l₁ l₂) (min u₁ u₂) ⊆ Set.Icc l₁ u₁ ∩ Set.Icc l₂ u₂ :=
  ⟨⟨max_le h₁.1 h₂.1, le_min h₁.2 h₂.2⟩, fun _ hy =>
    ⟨⟨le_trans (le_max_left _ _) hy.1, le_trans hy.2 (min_le_left _ _)⟩,
     ⟨le_trans (le_max_right _ _) hy.1, le_trans hy.2 (min_le_right _ _)⟩⟩⟩

variable {C : Type*} [Fintype C] [DecidableEq C] [Nonempty C] (β : ℝ) (E : C → ℝ)

omit [DecidableEq C] in
/-- **Theorem 4 (naive enclosure).** Separate certified bounds on `Z_A` and on `Z`. -/
theorem prob_mem_naive (A : Finset C) {aL aU Zlo Zhi : ℝ} (haL : 0 ≤ aL) (hZlo : 0 < Zlo)
    (ha : ZOn β E A ∈ Set.Icc aL aU) (hZ : Z β E ∈ Set.Icc Zlo Zhi) :
    prob β E A ∈ Set.Icc (aL / Zhi) (aU / Zlo) :=
  ⟨div_le_div₀ (ZOn_nonneg β E A) ha.1 (Z_pos β E) hZ.2,
   div_le_div₀ (haL.trans (ha.1.trans ha.2)) ha.2 hZlo hZ.1⟩

/-- **Theorem 4 (coupled enclosure).** Certified bounds `aL ≤ Z_A ≤ aU`,
`bL ≤ Z_{Aᶜ} ≤ bU` with `0 ≤ aL, 0 ≤ bL` give `aL/(aL+bU) ≤ P(A) ≤ aU/(aU+bL)`. -/
theorem prob_mem_coupled (A : Finset C) {aL aU bL bU : ℝ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : ZOn β E A ∈ Set.Icc aL aU) (hb : ZOn β E Aᶜ ∈ Set.Icc bL bU) :
    prob β E A ∈ Set.Icc (aL / (aL + bU)) (aU / (aU + bL)) := by
  rw [prob_eq_ratio]
  exact ratio_mem_coupled haL hbL ha hb (by rw [ZOn_add_ZOn_compl]; exact Z_pos β E)

/-- **Theorem 5 (robust refinement).** Given two independently certified pairs of
intervals for `(Z_A, Z_{Aᶜ})` (e.g. before and after a refinement step), the coupled
enclosure built from their intersections is valid and contained in the coupled enclosure
of the first pair.  In particular refinement never invalidates the enclosure. -/
theorem prob_enclosure_refine (A : Finset C) {aL aU bL bU aL' aU' bL' bU' : ℝ}
    (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : ZOn β E A ∈ Set.Icc aL aU) (hb : ZOn β E Aᶜ ∈ Set.Icc bL bU)
    (ha' : ZOn β E A ∈ Set.Icc aL' aU') (hb' : ZOn β E Aᶜ ∈ Set.Icc bL' bU') :
    prob β E A ∈ Set.Icc (max aL aL' / (max aL aL' + min bU bU'))
        (min aU aU' / (min aU aU' + max bL bL')) ∧
      Set.Icc (max aL aL' / (max aL aL' + min bU bU'))
          (min aU aU' / (min aU aU' + max bL bL')) ⊆
        Set.Icc (aL / (aL + bU)) (aU / (aU + bL)) := by
  obtain ⟨hA, -⟩ := mem_Icc_max_min ha ha'
  obtain ⟨hB, -⟩ := mem_Icc_max_min hb hb'
  refine ⟨prob_mem_coupled β E A (le_max_of_le_left haL) (le_max_of_le_left hbL) hA hB, ?_⟩
  by_cases hpos : 0 < aL + bL
  · exact coupled_enclosure_mono haL hbL hpos (le_max_left _ _) (min_le_left _ _)
      (le_max_left _ _) (min_le_left _ _) (hA.1.trans hA.2) (hB.1.trans hB.2)
  · -- degenerate data `aL = bL = 0`: compare endpoints directly
    have h0 : aL = 0 ∧ bL = 0 := by constructor <;> linarith
    obtain ⟨rfl, rfl⟩ := h0
    have hbU' : 0 ≤ min bU bU' := (ZOn_nonneg β E Aᶜ).trans hB.2
    have haU : 0 ≤ aU := (ZOn_nonneg β E A).trans ha.2
    apply Set.Icc_subset_Icc
    · rw [zero_div]
      exact div_nonneg (le_max_left _ _) (add_nonneg (le_max_left _ _) hbU')
    · have hmU : 0 ≤ min aU aU' := (ZOn_nonneg β E A).trans hA.2
      have hmL : 0 ≤ max (0:ℝ) bL' := le_max_left _ _
      rcases hmU.eq_or_lt with h | h
      · rw [← h, zero_div, add_zero]; exact div_nonneg haU haU
      · rw [add_zero, div_self (ne_of_gt (h.trans_le (min_le_left _ _)))]
        exact div_le_one_of_le₀ (le_add_of_nonneg_right hmL) (by linarith)

/-- **Theorem 4 (from partial enumeration).** Explored set `S`, block energy certificates
`ℓ ≤ E ≤ u` on the unexplored conformations.  The four interval endpoints for
`Z_A = Z_{A∩S} + Z_{A\S}` and `Z_{Aᶜ}` are computed from exact sums over explored
conformations plus block bounds, and the coupled enclosure follows. -/
theorem prob_mem_of_blocks {ι : Type*} (hβ : 0 ≤ β) (blk : C → ι) (ℓ u : ι → ℝ)
    (S A : Finset C) (hℓ : ∀ c ∉ S, ℓ (blk c) ≤ E c) (hu : ∀ c ∉ S, E c ≤ u (blk c)) :
    prob β E A ∈ Set.Icc
      ((ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) /
        ((ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) +
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk ℓ (Aᶜ \ S))))
      ((ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) /
        ((ZOn β E (A ∩ S) + blockBound β blk ℓ (A \ S)) +
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S)))) := by
  have key : ∀ B : Finset C, ZOn β E B ∈ Set.Icc
      (ZOn β E (B ∩ S) + blockBound β blk u (B \ S))
      (ZOn β E (B ∩ S) + blockBound β blk ℓ (B \ S)) := by
    intro B
    have hsplit : ZOn β E (B ∩ S) + ZOn β E (B \ S) = ZOn β E B :=
      Finset.sum_inter_add_sum_diff B S _
    have h1 := blockBound_le_ZOn E blk hβ u (B \ S)
      fun c hc => hu c (Finset.mem_sdiff.mp hc).2
    have h2 := ZOn_le_blockBound E blk hβ ℓ (B \ S)
      fun c hc => hℓ c (Finset.mem_sdiff.mp hc).2
    exact ⟨by linarith, by linarith⟩
  exact prob_mem_coupled β E A
    (add_nonneg (ZOn_nonneg β E _) (blockBound_nonneg β blk u _))
    (add_nonneg (ZOn_nonneg β E _) (blockBound_nonneg β blk u _)) (key A) (key Aᶜ)

end FoldBounds

namespace FoldBounds

variable {C : Type*} [Fintype C] [DecidableEq C] [Nonempty C] {ι : Type*}

/-- The certified lower end `aL + bL` of the block enclosure of `Z` is positive. -/
theorem blocks_lower_sum_pos (β : ℝ) (E : C → ℝ) (blk : C → ι) (u : ι → ℝ)
    (S A : Finset C) :
    0 < (ZOn β E (A ∩ S) + blockBound β blk u (A \ S)) +
      (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S)) := by
  classical
  set g : C → ℝ := fun c => if c ∈ S then weight β E c else Real.exp (-β * u (blk c))
  have hsplit : ∀ B : Finset C,
      ZOn β E (B ∩ S) + blockBound β blk u (B \ S) = ∑ c ∈ B, g c := by
    intro B
    rw [Finset.sum_ite, Finset.filter_mem_eq_inter, Finset.inter_comm,
      Finset.filter_not, ZOn, blockBound]
    congr 1
    refine Finset.sum_congr ?_ fun _ _ => rfl
    ext c
    simp only [Finset.mem_filter, Finset.mem_sdiff]
    tauto
  rw [hsplit, hsplit, Finset.sum_add_sum_compl]
  refine Finset.sum_pos (fun c _ => ?_) Finset.univ_nonempty
  simp only [g]
  split_ifs
  · exact weight_pos β E c
  · exact Real.exp_pos _

/-- **Convergence certificate for branch-and-bound.** If every unexplored block has energy
gap at most `δ`, the width of the certified probability enclosure of
`prob_mem_of_blocks` is at most `β δ · (unexplored upper mass) / (certified lower Z)`.
Hence refining blocks until `δ → 0` (or enumerating everything) makes the enclosure
collapse to the exact value. -/
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
          (ZOn β E (Aᶜ ∩ S) + blockBound β blk u (Aᶜ \ S))) := by
  have hle : ∀ B : Finset C, blockBound β blk u (B \ S) ≤ blockBound β blk ℓ (B \ S) :=
    fun B => blockBound_monotone_refine blk hβ blk u ℓ _ fun c hc =>
      ((hℓ c (Finset.mem_sdiff.mp hc).2).trans (hu c (Finset.mem_sdiff.mp hc).2))
  have hw : ∀ B : Finset C, blockBound β blk ℓ (B \ S) - blockBound β blk u (B \ S) ≤
      β * δ * blockBound β blk ℓ (B \ S) :=
    fun B => blockBound_sub_le_mul blk hβ ℓ u _ δ fun c hc => hgap c (Finset.mem_sdiff.mp hc).2
  have hpos := blocks_lower_sum_pos β E blk u S A
  refine (coupled_width_le (add_nonneg (ZOn_nonneg β E _) (blockBound_nonneg β blk u _))
    (add_nonneg (ZOn_nonneg β E _) (blockBound_nonneg β blk u _))
    (by linarith [hle A]) (by linarith [hle Aᶜ]) hpos).trans ?_
  apply div_le_div_of_nonneg_right _ hpos.le
  linarith [hw A, hw Aᶜ]

end FoldBounds
