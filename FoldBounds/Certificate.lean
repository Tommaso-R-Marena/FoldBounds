module

public import FoldBounds.FreeEnergy

/-!
# FoldBounds.Certificate — from real-number theorems to exact rational certificates

The theorems in the other files are statements about real numbers.  A computer
implementation cannot evaluate `exp` exactly, so a *certificate* consists of

1. rational enclosures `wlo c ≤ exp (-β E c) ≤ whi c` of individual Boltzmann weights
   (these must be justified by a verified `exp`-enclosure — here we use Mathlib's
   `Real.exp_one_gt_d9`/`Real.exp_one_lt_d9`; in a production checker this would be a
   verified interval-arithmetic library), and
2. **exact rational arithmetic** (`ℚ`, no rounding) for sums and quotients.

The soundness lemmas below show that the rational computation produces exactly the real
quantities appearing in the abstract theorems (`Rat.cast` commutes with `+`, `/`, `∑`),
so no floating-point error analysis is needed.  Floating-point (IEEE) implementations are
**not** covered by these theorems; they would additionally require directed rounding and
a separate proof.

The file ends with a fully checked numerical instance.
-/

@[expose] public section

open Finset Real

namespace FoldBounds

/-- Coupled lower bound computed in exact rational arithmetic. -/
def coupledLowerQ (aL bU : ℚ) : ℚ := aL / (aL + bU)

/-- Coupled upper bound computed in exact rational arithmetic. -/
def coupledUpperQ (aU bL : ℚ) : ℚ := aU / (aU + bL)

variable {C : Type*} [Fintype C] (β : ℝ) (E : C → ℝ)

omit [Fintype C] in
/-- Rational weight enclosures sum (exactly, in `ℚ`) to an enclosure of `Z_T`. -/
theorem ZOn_mem_Icc_of_rat (T : Finset C) (wlo whi : C → ℚ)
    (hw : ∀ c ∈ T, ((wlo c : ℚ) : ℝ) ≤ weight β E c ∧ weight β E c ≤ ((whi c : ℚ) : ℝ)) :
    ZOn β E T ∈ Set.Icc (((∑ c ∈ T, wlo c : ℚ)) : ℝ) (((∑ c ∈ T, whi c : ℚ)) : ℝ) := by
  push_cast
  exact ⟨Finset.sum_le_sum fun c hc => (hw c hc).1, Finset.sum_le_sum fun c hc => (hw c hc).2⟩

/-- **Soundness of the rational probability certificate.** If rational numbers
`aL ≤ Z_A ≤ aU`, `bL ≤ Z_{Aᶜ} ≤ bU` are certified (as real inequalities), then the
rationally computed `coupledLowerQ aL bU`, `coupledUpperQ aU bL` enclose `P(A)`. -/
theorem prob_mem_Icc_of_rat [DecidableEq C] [Nonempty C] (A : Finset C)
    {aL aU bL bU : ℚ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : ZOn β E A ∈ Set.Icc (aL : ℝ) aU) (hb : ZOn β E Aᶜ ∈ Set.Icc (bL : ℝ) bU) :
    prob β E A ∈ Set.Icc ((coupledLowerQ aL bU : ℚ) : ℝ) ((coupledUpperQ aU bL : ℚ) : ℝ) := by
  simp only [coupledLowerQ, coupledUpperQ]
  push_cast
  exact prob_mem_coupled β E A (by exact_mod_cast haL) (by exact_mod_cast hbL) ha hb

/-! ### A fully checked instance

Three conformations with energies `0, 1, 2` (in units of `k_B T`, so `β = 1`); `A` is the
ground state `{0}`.  True value: `P(A) = 1/(1 + e⁻¹ + e⁻²) ≈ 0.66524`. -/

/-- Toy energy landscape. -/
noncomputable def toyEnergy : Fin 3 → ℝ := ![0, 1, 2]

theorem exp_neg_one_mem : Real.exp (-1) ∈ Set.Icc (0.3678 : ℝ) 0.3679 := by
  have h1 := Real.exp_one_gt_d9
  have h2 := Real.exp_one_lt_d9
  have hpos := Real.exp_pos 1
  rw [Real.exp_neg]
  constructor
  · rw [le_inv_comm₀ (by norm_num) hpos]; norm_num at *; linarith
  · rw [inv_le_comm₀ hpos (by norm_num)]; norm_num at *; linarith

theorem exp_neg_two_mem : Real.exp (-2) ∈ Set.Icc (0.1352 : ℝ) 0.1354 := by
  have h := exp_neg_one_mem
  have e : Real.exp (-2) = Real.exp (-1) * Real.exp (-1) := by
    rw [← Real.exp_add]; norm_num
  rw [e]
  constructor <;> nlinarith [h.1, h.2]

/-- **Certified numerical result**: `0.6652 ≤ P({0}) ≤ 0.6654` for the toy landscape,
obtained from `prob_mem_Icc_of_rat` with exact rational arithmetic. -/
theorem toy_prob_mem : prob 1 toyEnergy {0} ∈ Set.Icc (0.6652 : ℝ) 0.6654 := by
  have hc : ({0}ᶜ : Finset (Fin 3)) = {1, 2} := by decide
  have hA : ZOn 1 toyEnergy {0} ∈ Set.Icc ((1 : ℚ) : ℝ) ((1 : ℚ) : ℝ) := by
    simp [ZOn, weight, toyEnergy]
  have hB : ZOn 1 toyEnergy {0}ᶜ ∈ Set.Icc (((0.5030 : ℚ)) : ℝ) ((0.5033 : ℚ) : ℝ) := by
    rw [hc]
    simp only [ZOn, weight, toyEnergy]
    rw [Finset.sum_pair (by decide)]
    simp only [Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
      Matrix.tail_cons, Matrix.cons_val_zero]
    have h1 := exp_neg_one_mem
    have h2 := exp_neg_two_mem
    norm_num at h1 h2 ⊢
    constructor <;> linarith [h1.1, h1.2, h2.1, h2.2]
  have key := prob_mem_Icc_of_rat 1 toyEnergy {0} (by norm_num) (by norm_num) hA hB
  have hlo : (0.6652 : ℝ) ≤ ((coupledLowerQ 1 0.5033 : ℚ) : ℝ) := by
    rw [show (0.6652 : ℝ) = ((0.6652 : ℚ) : ℝ) by norm_num, Rat.cast_le]
    norm_num [coupledLowerQ]
  have hhi : ((coupledUpperQ 1 0.5030 : ℚ) : ℝ) ≤ 0.6654 := by
    rw [show (0.6654 : ℝ) = ((0.6654 : ℚ) : ℝ) by norm_num, Rat.cast_le]
    norm_num [coupledUpperQ]
  exact ⟨hlo.trans key.1, key.2.trans hhi⟩

end FoldBounds
