module

public import Mathlib

/-!
# FoldBounds.Basic — Boltzmann partition functions on finite conformation spaces

Setting: `C` is a finite type of conformations, `E : C → ℝ` an energy function and
`β : ℝ` an inverse temperature.  We define

* `weight β E c = exp (-β * E c)`  (Boltzmann weight),
* `Z β E = ∑ c, weight β E c`       (partition function),
* `ZOn β E A = ∑ c ∈ A, weight β E c` (partial partition function of `A : Finset C`),
* `prob β E A = ZOn β E A / Z β E`  (equilibrium probability of `A`).

None of the results in this file needs `0 ≤ β`: Boltzmann weights are positive for
every real `β`.  Sign assumptions on `β` only enter when energy bounds are turned into
weight bounds (see `FoldBounds.PartitionBounds`).
-/

@[expose] public section

open Finset Real

namespace FoldBounds

variable {C : Type*}

/-- Boltzmann weight `exp (-β E(c))` of a conformation `c`. -/
noncomputable def weight (β : ℝ) (E : C → ℝ) (c : C) : ℝ := Real.exp (-β * E c)

/-- Partial partition function `Z_A = ∑_{c ∈ A} exp (-β E(c))`. -/
noncomputable def ZOn (β : ℝ) (E : C → ℝ) (A : Finset C) : ℝ := ∑ c ∈ A, weight β E c

theorem weight_pos (β : ℝ) (E : C → ℝ) (c : C) : 0 < weight β E c := Real.exp_pos _

@[simp] theorem ZOn_empty (β : ℝ) (E : C → ℝ) : ZOn β E ∅ = 0 := by simp [ZOn]

theorem ZOn_union_of_disjoint [DecidableEq C] (β : ℝ) (E : C → ℝ) {A B : Finset C}
    (h : Disjoint A B) : ZOn β E (A ∪ B) = ZOn β E A + ZOn β E B :=
  Finset.sum_union h

theorem ZOn_nonneg (β : ℝ) (E : C → ℝ) (A : Finset C) : 0 ≤ ZOn β E A :=
  Finset.sum_nonneg fun c _ => (weight_pos β E c).le

theorem ZOn_pos (β : ℝ) (E : C → ℝ) {A : Finset C} (hA : A.Nonempty) : 0 < ZOn β E A :=
  Finset.sum_pos (fun c _ => weight_pos β E c) hA

theorem ZOn_mono (β : ℝ) (E : C → ℝ) {A B : Finset C} (h : A ⊆ B) :
    ZOn β E A ≤ ZOn β E B :=
  Finset.sum_le_sum_of_subset_of_nonneg h fun c _ _ => (weight_pos β E c).le

variable [Fintype C]

/-- Partition function `Z = ∑_c exp (-β E(c))`. -/
noncomputable def Z (β : ℝ) (E : C → ℝ) : ℝ := ∑ c, weight β E c

/-- Equilibrium (Boltzmann) probability `P(A) = Z_A / Z`. -/
noncomputable def prob (β : ℝ) (E : C → ℝ) (A : Finset C) : ℝ := ZOn β E A / Z β E

variable (β : ℝ) (E : C → ℝ)

@[simp] theorem ZOn_univ : ZOn β E Finset.univ = Z β E := rfl

/-- **Theorem 1.** The partition function of a nonempty finite conformation space is
positive (for every real `β`). -/
theorem Z_pos [Nonempty C] : 0 < Z β E :=
  ZOn_pos β E Finset.univ_nonempty

theorem ZOn_le_Z (A : Finset C) : ZOn β E A ≤ Z β E :=
  ZOn_mono β E (Finset.subset_univ A)

/-- `Z = Z_A + Z_{Aᶜ}`: the basic decomposition behind every bound in this library. -/
theorem ZOn_add_ZOn_compl [DecidableEq C] (A : Finset C) :
    ZOn β E A + ZOn β E Aᶜ = Z β E :=
  Finset.sum_add_sum_compl A _

theorem prob_nonneg (A : Finset C) : 0 ≤ prob β E A :=
  div_nonneg (ZOn_nonneg β E A) (ZOn_nonneg β E Finset.univ)

theorem prob_le_one (A : Finset C) : prob β E A ≤ 1 :=
  div_le_one_of_le₀ (ZOn_le_Z β E A) (ZOn_nonneg β E Finset.univ)

/-- **Theorem 2.** `0 ≤ P(A) ≤ 1`. -/
theorem prob_mem_Icc (A : Finset C) : prob β E A ∈ Set.Icc (0 : ℝ) 1 :=
  ⟨prob_nonneg β E A, prob_le_one β E A⟩

theorem prob_univ [Nonempty C] : prob β E Finset.univ = 1 :=
  div_self (Z_pos β E).ne'

@[simp] theorem prob_empty : prob β E ∅ = 0 := by simp [prob]

theorem prob_pos [Nonempty C] {A : Finset C} (hA : A.Nonempty) : 0 < prob β E A :=
  div_pos (ZOn_pos β E hA) (Z_pos β E)

theorem prob_union_of_disjoint [DecidableEq C] {A B : Finset C} (h : Disjoint A B) :
    prob β E (A ∪ B) = prob β E A + prob β E B := by
  simp only [prob, ZOn_union_of_disjoint β E h, add_div]

theorem prob_compl [DecidableEq C] [Nonempty C] (A : Finset C) :
    prob β E Aᶜ = 1 - prob β E A := by
  have h := ZOn_add_ZOn_compl β E A
  have hZ := (Z_pos β E).ne'
  simp only [prob]
  field_simp
  linarith

/-- `P(A) = Z_A / (Z_A + Z_{Aᶜ})`: the numerator also appears in the denominator.  This
dependency is what the coupled probability bounds exploit. -/
theorem prob_eq_ratio [DecidableEq C] (A : Finset C) :
    prob β E A = ZOn β E A / (ZOn β E A + ZOn β E Aᶜ) := by
  rw [ZOn_add_ZOn_compl, prob]

end FoldBounds
