module

public import FoldBounds.Probability

/-!
# FoldBounds.CoupledRatio — sharp coupled linear-fractional enclosures (CEC §1)

Regions `k ∈ s` (a `Finset ι`) carry mass bounds `L k ≤ α k ≤ U k` on the A-mass `α k` and
ratio bounds `r k * α k ≤ β k ≤ R k * α k` on the B-mass `β k`.  The target quantity is
`Q(α, β) = (∑ k ∈ s, β k) / (∑ k ∈ s, α k)`.  Following `docs/coupled/CEC_THEORY.md` §1:

* **Theorems 1.1 / 1.2** (`forall_le_Q_iff`, `forall_Q_le_iff`): the exact lower and upper
  tests `0 ≤ gLo t` and `gHi t ≤ 0`.  Under the §1.1 data conditions (`WellFormed`) no
  non-emptiness hypothesis is needed; the soundness directions need no hypothesis at all.
* **Theorem 1.3** (partial: `isLeast_lamLo`, `isGreatest_lamHi`, `lamLo_attained_at_vertex`,
  `gLo_lamLo_eq_zero`, …): the extremes `lamLo = sInf QSet`, `lamHi = sSup QSet` exist and are
  attained when some `U k > 0`; they are attained at vertices `α k ∈ {L k, U k}` and are the
  roots of the tests.  The sorted-threshold description is not formalised.
* **Theorem 1.4** (`sepLo_le_lamLo`, `lamHi_le_sepHi`, `minRatio_le_lamLo`,
  `lamHi_le_maxRatio`): dominance over the separate enclosure and the extreme regional ratios.
* **Theorem 1.5** (`lamLo_eq_sepLo_iff`, `lamHi_eq_sepHi_iff`): when coupling does not help.
* **Theorem 1.6** (`lamLo_lamHi_of_uniform`): uniform-ratio exactness.
* **Theorem 1.7** (`QSet_mono`, `lamLo_mono`, `lamHi_anti`, `gLo_mono`, `QSet_refine`,
  `lamLo_refine`, `lamHi_refine`): monotone refinement and refinement by splitting.
* **§1.4, T1′** (`le_Q_of_lowerDualCert`, `Q_le_of_upperDualCert`): soundness of LP dual
  certificates over a polytope of A-masses.
* An **executable rational checker** (`checkLower`, `checkUpper`, `checkEnclosure`) with
  soundness and completeness theorems about the actual functions (`checkLower_sound`,
  `checkLower_iff`, `checkEnclosure_sound`, …).
-/

@[expose] public section

open Finset

namespace FoldBounds
namespace CoupledRatio

variable {ι : Type*}

/-! ### Data, admissible set, test functions -/

section Defs

variable (s : Finset ι) (L U r R : ι → ℝ)

/-- The §1.1 data conditions on every region: `0 ≤ L k ≤ U k` and `0 ≤ r k ≤ R k`. -/
def WellFormed : Prop :=
  ∀ k ∈ s, 0 ≤ L k ∧ L k ≤ U k ∧ 0 ≤ r k ∧ r k ≤ R k

/-- The §1.1 admissible set: nonnegative `(α, β)` with `L ≤ α ≤ U`, `r α ≤ β ≤ R α` on every
region and positive total A-mass. -/
def IsAdm (α β : ι → ℝ) : Prop :=
  (∀ k ∈ s, 0 ≤ α k ∧ 0 ≤ β k ∧ L k ≤ α k ∧ α k ≤ U k ∧ r k * α k ≤ β k ∧ β k ≤ R k * α k) ∧
    0 < ∑ k ∈ s, α k

/-- The target ratio `Q(α, β) = Σ β / Σ α`. -/
noncomputable def Q (α β : ι → ℝ) : ℝ := (∑ k ∈ s, β k) / ∑ k ∈ s, α k

/-- All values of `Q` on the admissible set. -/
def QSet : Set ℝ := {q | ∃ α β, IsAdm s L U r R α β ∧ Q s α β = q}

/-- Lower test function `g₋(t) = Σ_k min((r_k - t) L_k, (r_k - t) U_k)`. -/
noncomputable def gLo (t : ℝ) : ℝ := ∑ k ∈ s, min ((r k - t) * L k) ((r k - t) * U k)

/-- Upper test function `g₊(t) = Σ_k max((R_k - t) L_k, (R_k - t) U_k)`. -/
noncomputable def gHi (t : ℝ) : ℝ := ∑ k ∈ s, max ((R k - t) * L k) ((R k - t) * U k)

/-- The sharp lower endpoint `λ₋ = inf Q` over the admissible set. -/
noncomputable def lamLo : ℝ := sInf (QSet s L U r R)

/-- The sharp upper endpoint `λ₊ = sup Q` over the admissible set. -/
noncomputable def lamHi : ℝ := sSup (QSet s L U r R)

/-- Separate lower endpoint `S₋ = Σ r L / Σ U` (bounding `Z_B` and `Z_A` separately). -/
noncomputable def sepLo : ℝ := (∑ k ∈ s, r k * L k) / ∑ k ∈ s, U k

/-- Separate upper endpoint `S₊ = Σ R U / Σ L`. -/
noncomputable def sepHi : ℝ := (∑ k ∈ s, R k * U k) / ∑ k ∈ s, L k

end Defs

variable {s : Finset ι} {L U r R : ι → ℝ}

/-! ### Elementary one-region facts -/

/-- A linear function on `[a, b]` is at least the smaller endpoint value. -/
theorem min_mul_le_mul {c a b x : ℝ} (ha : a ≤ x) (hb : x ≤ b) :
    min (c * a) (c * b) ≤ c * x := by
  rcases le_total 0 c with hc | hc
  · exact (min_le_left _ _).trans (mul_le_mul_of_nonneg_left ha hc)
  · exact (min_le_right _ _).trans (mul_le_mul_of_nonpos_left hb hc)

/-- A linear function on `[a, b]` is at most the larger endpoint value. -/
theorem mul_le_max_mul {c a b x : ℝ} (ha : a ≤ x) (hb : x ≤ b) :
    c * x ≤ max (c * a) (c * b) := by
  rcases le_total 0 c with hc | hc
  · exact (mul_le_mul_of_nonneg_left hb hc).trans (le_max_right _ _)
  · exact (mul_le_mul_of_nonpos_left ha hc).trans (le_max_left _ _)

theorem min_mul_eq_ite {c a b : ℝ} (hab : a ≤ b) :
    min (c * a) (c * b) = c * (if 0 ≤ c then a else b) := by
  split_ifs with hc
  · exact min_eq_left (mul_le_mul_of_nonneg_left hab hc)
  · exact min_eq_right (mul_le_mul_of_nonpos_left hab (le_of_lt (not_le.mp hc)))

theorem max_mul_eq_ite {c a b : ℝ} (hab : a ≤ b) :
    max (c * a) (c * b) = c * (if 0 ≤ c then b else a) := by
  split_ifs with hc
  · exact max_eq_right (mul_le_mul_of_nonneg_left hab hc)
  · exact max_eq_left (mul_le_mul_of_nonpos_left hab (le_of_lt (not_le.mp hc)))

/-- The second form of `g₋` in the spec: `Σ (r_k - t) · (L_k if r_k ≥ t else U_k)`. -/
theorem gLo_eq_ite (hLU : ∀ k ∈ s, L k ≤ U k) (t : ℝ) :
    gLo s L U r t = ∑ k ∈ s, (r k - t) * (if t ≤ r k then L k else U k) := by
  refine sum_congr rfl fun k hk => ?_
  rw [min_mul_eq_ite (hLU k hk)]
  simp only [sub_nonneg]

/-- The second form of `g₊` in the spec: `Σ (R_k - t) · (U_k if R_k ≥ t else L_k)`. -/
theorem gHi_eq_ite (hLU : ∀ k ∈ s, L k ≤ U k) (t : ℝ) :
    gHi s L U R t = ∑ k ∈ s, (R k - t) * (if t ≤ R k then U k else L k) := by
  refine sum_congr rfl fun k hk => ?_
  rw [max_mul_eq_ite (hLU k hk)]
  simp only [sub_nonneg]

/-- `Σ (c_k - t) x_k = Σ c_k x_k - t Σ x_k`. -/
theorem sum_sub_mul_eq (c x : ι → ℝ) (t : ℝ) :
    ∑ k ∈ s, (c k - t) * x k = ∑ k ∈ s, c k * x k - t * ∑ k ∈ s, x k := by
  rw [mul_sum, ← sum_sub_distrib]
  exact sum_congr rfl fun k _ => by ring

/-! ### Theorems 1.1 and 1.2: the exact tests -/

/-- **Theorem 1.1 (⇐), soundness of the lower test.** No hypothesis on the data is needed. -/
theorem le_Q_of_gLo_nonneg {t : ℝ} (hg : 0 ≤ gLo s L U r t) {α β : ι → ℝ}
    (h : IsAdm s L U r R α β) : t ≤ Q s α β := by
  obtain ⟨hk, hpos⟩ := h
  rw [Q, le_div_iff₀ hpos]
  have key : gLo s L U r t ≤ ∑ k ∈ s, (r k - t) * α k := by
    refine sum_le_sum fun k hks => ?_
    obtain ⟨_, _, hL, hU, _, _⟩ := hk k hks
    exact min_mul_le_mul hL hU
  have hrb : ∑ k ∈ s, r k * α k ≤ ∑ k ∈ s, β k :=
    sum_le_sum fun k hks => (hk k hks).2.2.2.2.1
  rw [sum_sub_mul_eq] at key
  linarith

/-- **Theorem 1.2 (⇐), soundness of the upper test.** No hypothesis on the data is needed. -/
theorem Q_le_of_gHi_nonpos {t : ℝ} (hg : gHi s L U R t ≤ 0) {α β : ι → ℝ}
    (h : IsAdm s L U r R α β) : Q s α β ≤ t := by
  obtain ⟨hk, hpos⟩ := h
  rw [Q, div_le_iff₀ hpos]
  have key : ∑ k ∈ s, (R k - t) * α k ≤ gHi s L U R t := by
    refine sum_le_sum fun k hks => ?_
    obtain ⟨_, _, hL, hU, _, _⟩ := hk k hks
    exact mul_le_max_mul hL hU
  have hRb : ∑ k ∈ s, β k ≤ ∑ k ∈ s, R k * α k :=
    sum_le_sum fun k hks => (hk k hks).2.2.2.2.2
  rw [sum_sub_mul_eq] at key
  linarith

/-- A positive sum has a positive term. -/
theorem exists_pos_of_sum_pos {x : ι → ℝ} (h : 0 < ∑ k ∈ s, x k) : ∃ k ∈ s, 0 < x k := by
  have h' : ∑ k ∈ s, (0 : ℝ) < ∑ k ∈ s, x k := by simpa using h
  exact exists_lt_of_sum_lt h'

/-- If a nonnegative family has zero sum, `Σ (c_k - t) x_k = 0`. -/
theorem sum_sub_mul_eq_zero_of_sum_eq_zero {c x : ι → ℝ} {t : ℝ} (hx : ∀ k ∈ s, 0 ≤ x k)
    (h0 : ¬ 0 < ∑ k ∈ s, x k) : ∑ k ∈ s, (c k - t) * x k = 0 := by
  have hsum : ∑ k ∈ s, x k = 0 := le_antisymm (not_lt.mp h0) (sum_nonneg hx)
  rw [sum_eq_zero_iff_of_nonneg hx] at hsum
  exact sum_eq_zero fun k hk => by rw [hsum k hk, mul_zero]

/-- **Theorem 1.1 (⇒), completeness of the lower test.** If every admissible configuration has
`t ≤ Q`, then `0 ≤ g₋(t)`.  Only the §1.1 data conditions are used: when the admissible set is
empty (every `U k = 0`) the test function vanishes identically. -/
theorem gLo_nonneg_of_forall_le_Q (hwf : WellFormed s L U r R) {t : ℝ}
    (h : ∀ α β, IsAdm s L U r R α β → t ≤ Q s α β) : 0 ≤ gLo s L U r t := by
  by_contra hneg
  push_neg at hneg
  set α : ι → ℝ := fun k => if t ≤ r k then L k else U k with hαdef
  have hg : gLo s L U r t = ∑ k ∈ s, (r k - t) * α k :=
    gLo_eq_ite (fun k hk => (hwf k hk).2.1) t
  have hαnn : ∀ k ∈ s, 0 ≤ α k := fun k hk => by
    obtain ⟨h0, h1, _, _⟩ := hwf k hk
    simp only [hαdef]; split_ifs <;> linarith
  have hpos : 0 < ∑ k ∈ s, α k := by
    by_contra h0
    have := sum_sub_mul_eq_zero_of_sum_eq_zero (c := r) (t := t) hαnn h0
    linarith
  have hadm : IsAdm s L U r R α (fun k => r k * α k) := by
    refine ⟨fun k hk => ?_, hpos⟩
    obtain ⟨h0, h1, h2, h3⟩ := hwf k hk
    have ha := hαnn k hk
    refine ⟨ha, mul_nonneg h2 ha, ?_, ?_, le_rfl, mul_le_mul_of_nonneg_right h3 ha⟩ <;>
      simp only [hαdef] <;> split_ifs <;> linarith
  have hle := h α _ hadm
  rw [Q, le_div_iff₀ hpos] at hle
  rw [hg, sum_sub_mul_eq] at hneg
  linarith

/-- **Theorem 1.2 (⇒), completeness of the upper test.** -/
theorem gHi_nonpos_of_forall_Q_le (hwf : WellFormed s L U r R) {t : ℝ}
    (h : ∀ α β, IsAdm s L U r R α β → Q s α β ≤ t) : gHi s L U R t ≤ 0 := by
  by_contra hneg
  push_neg at hneg
  set α : ι → ℝ := fun k => if t ≤ R k then U k else L k with hαdef
  have hg : gHi s L U R t = ∑ k ∈ s, (R k - t) * α k :=
    gHi_eq_ite (fun k hk => (hwf k hk).2.1) t
  have hαnn : ∀ k ∈ s, 0 ≤ α k := fun k hk => by
    obtain ⟨h0, h1, _, _⟩ := hwf k hk
    simp only [hαdef]; split_ifs <;> linarith
  have hpos : 0 < ∑ k ∈ s, α k := by
    by_contra h0
    have := sum_sub_mul_eq_zero_of_sum_eq_zero (c := R) (t := t) hαnn h0
    linarith
  have hadm : IsAdm s L U r R α (fun k => R k * α k) := by
    refine ⟨fun k hk => ?_, hpos⟩
    obtain ⟨h0, h1, h2, h3⟩ := hwf k hk
    have ha := hαnn k hk
    refine ⟨ha, mul_nonneg (h2.trans h3) ha, ?_, ?_, mul_le_mul_of_nonneg_right h3 ha, le_rfl⟩ <;>
      simp only [hαdef] <;> split_ifs <;> linarith
  have hle := h α _ hadm
  rw [Q, div_le_iff₀ hpos] at hle
  rw [hg, sum_sub_mul_eq] at hneg
  linarith

/-- **Theorem 1.1 (exact lower test).** Under the §1.1 data conditions, for every real `t`:
`(∀ (α, β) admissible, t ≤ Q(α, β)) ↔ 0 ≤ g₋(t)`.  The spec's hypothesis `Adm ≠ ∅` is not
needed. -/
theorem forall_le_Q_iff (hwf : WellFormed s L U r R) (t : ℝ) :
    (∀ α β, IsAdm s L U r R α β → t ≤ Q s α β) ↔ 0 ≤ gLo s L U r t :=
  ⟨gLo_nonneg_of_forall_le_Q hwf, fun hg _ _ h => le_Q_of_gLo_nonneg hg h⟩

/-- **Theorem 1.2 (exact upper test).** Under the §1.1 data conditions, for every real `t`:
`(∀ (α, β) admissible, Q(α, β) ≤ t) ↔ g₊(t) ≤ 0`.  The spec's hypothesis `Adm ≠ ∅` is not
needed. -/
theorem forall_Q_le_iff (hwf : WellFormed s L U r R) (t : ℝ) :
    (∀ α β, IsAdm s L U r R α β → Q s α β ≤ t) ↔ gHi s L U R t ≤ 0 :=
  ⟨gHi_nonpos_of_forall_Q_le hwf, fun hg _ _ h => Q_le_of_gHi_nonpos hg h⟩

/-- Theorem 1.1 in the spec's literal form, with its (superfluous) hypothesis `Adm ≠ ∅`. -/
theorem forall_le_Q_iff_of_nonempty (hwf : WellFormed s L U r R)
    (_hne : ∃ α β, IsAdm s L U r R α β) (t : ℝ) :
    (∀ α β, IsAdm s L U r R α β → t ≤ Q s α β) ↔ 0 ≤ gLo s L U r t :=
  forall_le_Q_iff hwf t

/-- Theorem 1.2 in the spec's literal form, with its (superfluous) hypothesis `Adm ≠ ∅`. -/
theorem forall_Q_le_iff_of_nonempty (hwf : WellFormed s L U r R)
    (_hne : ∃ α β, IsAdm s L U r R α β) (t : ℝ) :
    (∀ α β, IsAdm s L U r R α β → Q s α β ≤ t) ↔ gHi s L U R t ≤ 0 :=
  forall_Q_le_iff hwf t

/-- Precise meaning of the non-emptiness hypothesis: under the §1.1 data conditions the
admissible set is nonempty iff some region has `U k > 0`. -/
theorem adm_nonempty_iff (hwf : WellFormed s L U r R) :
    (∃ α β, IsAdm s L U r R α β) ↔ ∃ k ∈ s, 0 < U k := by
  constructor
  · rintro ⟨α, β, hk, hpos⟩
    obtain ⟨k, hks, hαk⟩ := exists_pos_of_sum_pos hpos
    exact ⟨k, hks, hαk.trans_le (hk k hks).2.2.2.1⟩
  · rintro ⟨k, hks, hUk⟩
    refine ⟨U, fun j => r j * U j, fun j hj => ?_, ?_⟩
    · obtain ⟨h0, h1, h2, h3⟩ := hwf j hj
      have hU0 : 0 ≤ U j := h0.trans h1
      exact ⟨hU0, mul_nonneg h2 hU0, h1, le_rfl, le_rfl, mul_le_mul_of_nonneg_right h3 hU0⟩
    · exact sum_pos' (fun j hj => (hwf j hj).1.trans (hwf j hj).2.1) ⟨k, hks, hUk⟩

/-! ### Theorem 1.3 (partial): the extremes exist and are attained at vertices -/

section Vertices

variable [DecidableEq ι]

/-- The vertex of the mass box that takes the value `U` on `S` and `L` off `S`. -/
def vertexMass (L U : ι → ℝ) (S : Finset ι) (k : ι) : ℝ := if k ∈ S then U k else L k

theorem vertexMass_mem (L U : ι → ℝ) (S : Finset ι) (k : ι) :
    vertexMass L U S k = L k ∨ vertexMass L U S k = U k := by
  unfold vertexMass; split_ifs <;> simp

theorem vertexMass_nonneg (hwf : WellFormed s L U r R) (S : Finset ι) :
    ∀ k ∈ s, 0 ≤ vertexMass L U S k := fun k hk => by
  obtain ⟨h0, h1, _, _⟩ := hwf k hk
  unfold vertexMass; split_ifs <;> linarith

theorem vertexMass_isAdm_lo (hwf : WellFormed s L U r R) (S : Finset ι)
    (hpos : 0 < ∑ k ∈ s, vertexMass L U S k) :
    IsAdm s L U r R (vertexMass L U S) (fun k => r k * vertexMass L U S k) := by
  refine ⟨fun k hk => ?_, hpos⟩
  obtain ⟨h0, h1, h2, h3⟩ := hwf k hk
  have ha := vertexMass_nonneg hwf S k hk
  refine ⟨ha, mul_nonneg h2 ha, ?_, ?_, le_rfl, mul_le_mul_of_nonneg_right h3 ha⟩ <;>
    unfold vertexMass <;> split_ifs <;> linarith

theorem vertexMass_isAdm_hi (hwf : WellFormed s L U r R) (S : Finset ι)
    (hpos : 0 < ∑ k ∈ s, vertexMass L U S k) :
    IsAdm s L U r R (vertexMass L U S) (fun k => R k * vertexMass L U S k) := by
  refine ⟨fun k hk => ?_, hpos⟩
  obtain ⟨h0, h1, h2, h3⟩ := hwf k hk
  have ha := vertexMass_nonneg hwf S k hk
  refine ⟨ha, mul_nonneg (h2.trans h3) ha, ?_, ?_, mul_le_mul_of_nonneg_right h3 ha, le_rfl⟩ <;>
    unfold vertexMass <;> split_ifs <;> linarith

/-- Dinkelbach step: a value below the ratio of every positive-mass vertex passes the lower
test.  (The test is evaluated at the vertex `U` on `{r < m}`, `L` elsewhere.) -/
theorem gLo_nonneg_of_vertex_bound (hwf : WellFormed s L U r R) {m : ℝ}
    (hv : ∀ S ⊆ s, 0 < ∑ k ∈ s, vertexMass L U S k →
      m ≤ Q s (vertexMass L U S) (fun k => r k * vertexMass L U S k)) :
    0 ≤ gLo s L U r m := by
  set S0 := s.filter (fun k => r k < m) with hS0
  have hterm : gLo s L U r m = ∑ k ∈ s, (r k - m) * vertexMass L U S0 k := by
    refine sum_congr rfl fun k hk => ?_
    rw [min_mul_eq_ite (hwf k hk).2.1]
    have hmem : k ∈ S0 ↔ r k < m := by simp [hS0, hk]
    unfold vertexMass
    by_cases hrk : r k < m
    · rw [if_neg (by linarith), if_pos (hmem.mpr hrk)]
    · rw [if_pos (by linarith), if_neg (fun h => hrk (hmem.mp h))]
  rw [hterm]
  by_cases hpos : 0 < ∑ k ∈ s, vertexMass L U S0 k
  · have hle := hv S0 (filter_subset _ _) hpos
    rw [Q, le_div_iff₀ hpos] at hle
    rw [sum_sub_mul_eq]
    linarith
  · rw [sum_sub_mul_eq_zero_of_sum_eq_zero (vertexMass_nonneg hwf S0) hpos]

/-- Dinkelbach step for the upper test. -/
theorem gHi_nonpos_of_vertex_bound (hwf : WellFormed s L U r R) {m : ℝ}
    (hv : ∀ S ⊆ s, 0 < ∑ k ∈ s, vertexMass L U S k →
      Q s (vertexMass L U S) (fun k => R k * vertexMass L U S k) ≤ m) :
    gHi s L U R m ≤ 0 := by
  set S0 := s.filter (fun k => m < R k) with hS0
  have hterm : gHi s L U R m = ∑ k ∈ s, (R k - m) * vertexMass L U S0 k := by
    refine sum_congr rfl fun k hk => ?_
    rw [max_mul_eq_ite (hwf k hk).2.1]
    have hmem : k ∈ S0 ↔ m < R k := by simp [hS0, hk]
    unfold vertexMass
    by_cases hRk : m < R k
    · rw [if_pos (by linarith), if_pos (hmem.mpr hRk)]
    · rcases (not_lt.mp hRk).lt_or_eq with h | h
      · rw [if_neg (by linarith), if_neg (fun h' => hRk (hmem.mp h'))]
      · rw [h, sub_self, zero_mul, zero_mul]
  rw [hterm]
  by_cases hpos : 0 < ∑ k ∈ s, vertexMass L U S0 k
  · have hle := hv S0 (filter_subset _ _) hpos
    rw [Q, div_le_iff₀ hpos] at hle
    rw [sum_sub_mul_eq]
    linarith
  · rw [sum_sub_mul_eq_zero_of_sum_eq_zero (vertexMass_nonneg hwf S0) hpos]

end Vertices

/-- **Theorem 1.3 (existence, lower).** If some region has `U k > 0`, the minimum of `Q` over
the admissible set exists and is attained at a vertex `α k ∈ {L k, U k}` with `β = r α`. -/
theorem exists_isLeast_vertex (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    ∃ α : ι → ℝ, (∀ k ∈ s, α k = L k ∨ α k = U k) ∧
      IsAdm s L U r R α (fun k => r k * α k) ∧
      IsLeast (QSet s L U r R) (Q s α (fun k => r k * α k)) := by
  classical
  set C := s.powerset.filter (fun S => 0 < ∑ k ∈ s, vertexMass L U S k) with hC
  have hsC : s ∈ C := by
    refine mem_filter.mpr ⟨mem_powerset.mpr (subset_refl _), ?_⟩
    obtain ⟨k, hks, hUk⟩ := hne
    refine sum_pos' (vertexMass_nonneg hwf s) ⟨k, hks, ?_⟩
    simp only [vertexMass, if_pos hks]
    exact hUk
  obtain ⟨S, hSC, hmin⟩ := exists_min_image C
    (fun S => Q s (vertexMass L U S) (fun k => r k * vertexMass L U S k)) ⟨s, hsC⟩
  have hSpos : 0 < ∑ k ∈ s, vertexMass L U S k := (mem_filter.mp hSC).2
  refine ⟨vertexMass L U S, fun k _ => vertexMass_mem L U S k,
    vertexMass_isAdm_lo hwf S hSpos, ⟨⟨_, _, vertexMass_isAdm_lo hwf S hSpos, rfl⟩, ?_⟩⟩
  rintro q ⟨α, β, hadm, rfl⟩
  refine le_Q_of_gLo_nonneg (gLo_nonneg_of_vertex_bound hwf fun S' hS' hpos' => ?_) hadm
  exact hmin S' (mem_filter.mpr ⟨mem_powerset.mpr hS', hpos'⟩)

/-- **Theorem 1.3 (existence, upper).** -/
theorem exists_isGreatest_vertex (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    ∃ α : ι → ℝ, (∀ k ∈ s, α k = L k ∨ α k = U k) ∧
      IsAdm s L U r R α (fun k => R k * α k) ∧
      IsGreatest (QSet s L U r R) (Q s α (fun k => R k * α k)) := by
  classical
  set C := s.powerset.filter (fun S => 0 < ∑ k ∈ s, vertexMass L U S k) with hC
  have hsC : s ∈ C := by
    refine mem_filter.mpr ⟨mem_powerset.mpr (subset_refl _), ?_⟩
    obtain ⟨k, hks, hUk⟩ := hne
    refine sum_pos' (vertexMass_nonneg hwf s) ⟨k, hks, ?_⟩
    simp only [vertexMass, if_pos hks]
    exact hUk
  obtain ⟨S, hSC, hmax⟩ := exists_max_image C
    (fun S => Q s (vertexMass L U S) (fun k => R k * vertexMass L U S k)) ⟨s, hsC⟩
  have hSpos : 0 < ∑ k ∈ s, vertexMass L U S k := (mem_filter.mp hSC).2
  refine ⟨vertexMass L U S, fun k _ => vertexMass_mem L U S k,
    vertexMass_isAdm_hi hwf S hSpos, ⟨⟨_, _, vertexMass_isAdm_hi hwf S hSpos, rfl⟩, ?_⟩⟩
  rintro q ⟨α, β, hadm, rfl⟩
  refine Q_le_of_gHi_nonpos (gHi_nonpos_of_vertex_bound hwf fun S' hS' hpos' => ?_) hadm
  exact hmax S' (mem_filter.mpr ⟨mem_powerset.mpr hS', hpos'⟩)

/-- `λ₋ = lamLo` is the minimum of `Q` over the admissible set. -/
theorem isLeast_lamLo (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    IsLeast (QSet s L U r R) (lamLo s L U r R) := by
  obtain ⟨α, -, -, h⟩ := exists_isLeast_vertex hwf hne
  rw [lamLo, h.csInf_eq]
  exact h

/-- `λ₊ = lamHi` is the maximum of `Q` over the admissible set. -/
theorem isGreatest_lamHi (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    IsGreatest (QSet s L U r R) (lamHi s L U r R) := by
  obtain ⟨α, -, -, h⟩ := exists_isGreatest_vertex hwf hne
  rw [lamHi, h.csSup_eq]
  exact h

/-- **Theorem 1.3 (vertex attainment, lower).** -/
theorem lamLo_attained_at_vertex (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    ∃ α : ι → ℝ, (∀ k ∈ s, α k = L k ∨ α k = U k) ∧
      IsAdm s L U r R α (fun k => r k * α k) ∧
      Q s α (fun k => r k * α k) = lamLo s L U r R := by
  obtain ⟨α, hv, hadm, h⟩ := exists_isLeast_vertex hwf hne
  exact ⟨α, hv, hadm, by rw [lamLo, h.csInf_eq]⟩

/-- **Theorem 1.3 (vertex attainment, upper).** -/
theorem lamHi_attained_at_vertex (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    ∃ α : ι → ℝ, (∀ k ∈ s, α k = L k ∨ α k = U k) ∧
      IsAdm s L U r R α (fun k => R k * α k) ∧
      Q s α (fun k => R k * α k) = lamHi s L U r R := by
  obtain ⟨α, hv, hadm, h⟩ := exists_isGreatest_vertex hwf hne
  exact ⟨α, hv, hadm, by rw [lamHi, h.csSup_eq]⟩

/-- Certificate form of Theorem 1.1: `t` is a valid lower endpoint iff `t ≤ λ₋` iff
`0 ≤ g₋(t)`. -/
theorem le_lamLo_iff (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) (t : ℝ) :
    t ≤ lamLo s L U r R ↔ 0 ≤ gLo s L U r t := by
  rw [← forall_le_Q_iff hwf]
  obtain ⟨⟨α₀, β₀, hadm₀, hq₀⟩, hlow⟩ := isLeast_lamLo hwf hne
  exact ⟨fun ht α β hadm => ht.trans (hlow ⟨α, β, hadm, rfl⟩),
    fun h => hq₀ ▸ h α₀ β₀ hadm₀⟩

/-- Certificate form of Theorem 1.2: `λ₊ ≤ t` iff `g₊(t) ≤ 0`. -/
theorem lamHi_le_iff (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) (t : ℝ) :
    lamHi s L U r R ≤ t ↔ gHi s L U R t ≤ 0 := by
  rw [← forall_Q_le_iff hwf]
  obtain ⟨⟨α₀, β₀, hadm₀, hq₀⟩, hup⟩ := isGreatest_lamHi hwf hne
  exact ⟨fun ht α β hadm => (hup ⟨α, β, hadm, rfl⟩).trans ht,
    fun h => hq₀ ▸ h α₀ β₀ hadm₀⟩

/-- **Theorem 1.3: `g₋(λ₋) = 0` at the optimum.** -/
theorem gLo_lamLo_eq_zero (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    gLo s L U r (lamLo s L U r R) = 0 := by
  refine le_antisymm ?_ ((le_lamLo_iff hwf hne _).mp le_rfl)
  obtain ⟨⟨α, β, ⟨hk, hpos⟩, hq⟩, _⟩ := isLeast_lamLo hwf hne
  have key : gLo s L U r (lamLo s L U r R) ≤ ∑ k ∈ s, (r k - lamLo s L U r R) * α k :=
    sum_le_sum fun k hks => min_mul_le_mul (hk k hks).2.2.1 (hk k hks).2.2.2.1
  have hrb : ∑ k ∈ s, r k * α k ≤ ∑ k ∈ s, β k := sum_le_sum fun k hks => (hk k hks).2.2.2.2.1
  have hβ : ∑ k ∈ s, β k = lamLo s L U r R * ∑ k ∈ s, α k := by
    rw [← hq, Q, div_mul_cancel₀ _ hpos.ne']
  rw [sum_sub_mul_eq] at key
  linarith

/-- **Theorem 1.3: `g₊(λ₊) = 0` at the optimum.** -/
theorem gHi_lamHi_eq_zero (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    gHi s L U R (lamHi s L U r R) = 0 := by
  refine le_antisymm ((lamHi_le_iff hwf hne _).mp le_rfl) ?_
  obtain ⟨⟨α, β, ⟨hk, hpos⟩, hq⟩, _⟩ := isGreatest_lamHi hwf hne
  have key : ∑ k ∈ s, (R k - lamHi s L U r R) * α k ≤ gHi s L U R (lamHi s L U r R) :=
    sum_le_sum fun k hks => mul_le_max_mul (hk k hks).2.2.1 (hk k hks).2.2.2.1
  have hRb : ∑ k ∈ s, β k ≤ ∑ k ∈ s, R k * α k := sum_le_sum fun k hks => (hk k hks).2.2.2.2.2
  have hβ : ∑ k ∈ s, β k = lamHi s L U r R * ∑ k ∈ s, α k := by
    rw [← hq, Q, div_mul_cancel₀ _ hpos.ne']
  rw [sum_sub_mul_eq] at key
  linarith

/-! ### Theorem 1.4: dominance over separate enclosures -/

/-- The separate lower endpoint `S₋ = Σ r L / Σ U` is a valid lower bound. -/
theorem sepLo_le_Q (hwf : WellFormed s L U r R) {α β : ι → ℝ} (h : IsAdm s L U r R α β) :
    sepLo s L U r ≤ Q s α β := by
  obtain ⟨hk, hpos⟩ := h
  have hrL : 0 ≤ ∑ k ∈ s, r k * L k :=
    sum_nonneg fun k hk' => mul_nonneg (hwf k hk').2.2.1 (hwf k hk').1
  have hαU : ∑ k ∈ s, α k ≤ ∑ k ∈ s, U k := sum_le_sum fun k hk' => (hk k hk').2.2.2.1
  have hLb : ∑ k ∈ s, r k * L k ≤ ∑ k ∈ s, β k := sum_le_sum fun k hk' =>
    (mul_le_mul_of_nonneg_left (hk k hk').2.2.1 (hwf k hk').2.2.1).trans (hk k hk').2.2.2.2.1
  calc sepLo s L U r = (∑ k ∈ s, r k * L k) / ∑ k ∈ s, U k := rfl
    _ ≤ (∑ k ∈ s, r k * L k) / ∑ k ∈ s, α k := div_le_div_of_nonneg_left hrL hpos hαU
    _ ≤ Q s α β := div_le_div_of_nonneg_right hLb hpos.le

/-- The separate upper endpoint `S₊ = Σ R U / Σ L` is a valid upper bound when `Σ L > 0`.
(When `Σ L = 0` the separate endpoint is `+∞`; Lean's `x / 0 = 0` convention makes the
hypothesis necessary.) -/
theorem Q_le_sepHi (hwf : WellFormed s L U r R) (hL : 0 < ∑ k ∈ s, L k) {α β : ι → ℝ}
    (h : IsAdm s L U r R α β) : Q s α β ≤ sepHi s L U R := by
  obtain ⟨hk, hpos⟩ := h
  have hRU : 0 ≤ ∑ k ∈ s, R k * U k := sum_nonneg fun k hk' =>
    mul_nonneg ((hwf k hk').2.2.1.trans (hwf k hk').2.2.2) ((hwf k hk').1.trans (hwf k hk').2.1)
  have hLα : ∑ k ∈ s, L k ≤ ∑ k ∈ s, α k := sum_le_sum fun k hk' => (hk k hk').2.2.1
  have hbU : ∑ k ∈ s, β k ≤ ∑ k ∈ s, R k * U k := sum_le_sum fun k hk' =>
    (hk k hk').2.2.2.2.2.trans (mul_le_mul_of_nonneg_left (hk k hk').2.2.2.1
      ((hwf k hk').2.2.1.trans (hwf k hk').2.2.2))
  calc Q s α β ≤ (∑ k ∈ s, R k * U k) / ∑ k ∈ s, α k := div_le_div_of_nonneg_right hbU hpos.le
    _ ≤ sepHi s L U R := div_le_div_of_nonneg_left hRU hL hLα

/-- Any `m ≤ r k` on the regions with `U k > 0` is a valid lower bound. -/
theorem le_Q_of_le_ratios {m : ℝ} (hm : ∀ k ∈ s, 0 < U k → m ≤ r k) {α β : ι → ℝ}
    (h : IsAdm s L U r R α β) : m ≤ Q s α β := by
  obtain ⟨hk, hpos⟩ := h
  rw [Q, le_div_iff₀ hpos, mul_sum]
  refine sum_le_sum fun k hks => ?_
  obtain ⟨hα0, hβ0, _, hαU, hrβ, _⟩ := hk k hks
  by_cases hU : 0 < U k
  · exact (mul_le_mul_of_nonneg_right (hm k hks hU) hα0).trans hrβ
  · have hα : α k = 0 := le_antisymm (hαU.trans (not_lt.mp hU)) hα0
    rw [hα, mul_zero]
    exact hβ0

/-- Any `M ≥ R k` on the regions with `U k > 0` is a valid upper bound. -/
theorem Q_le_of_ratios_le {M : ℝ} (hM : ∀ k ∈ s, 0 < U k → R k ≤ M) {α β : ι → ℝ}
    (h : IsAdm s L U r R α β) : Q s α β ≤ M := by
  obtain ⟨hk, hpos⟩ := h
  rw [Q, div_le_iff₀ hpos, mul_sum]
  refine sum_le_sum fun k hks => ?_
  obtain ⟨hα0, _, _, hαU, _, hβR⟩ := hk k hks
  by_cases hU : 0 < U k
  · exact hβR.trans (mul_le_mul_of_nonneg_right (hM k hks hU) hα0)
  · have hα : α k = 0 := le_antisymm (hαU.trans (not_lt.mp hU)) hα0
    rw [hα, mul_zero] at hβR ⊢
    exact hβR

/-- **Theorem 1.4 (lower).** `λ₋ ≥ S₋`. -/
theorem sepLo_le_lamLo (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) :
    sepLo s L U r ≤ lamLo s L U r R := by
  obtain ⟨⟨α, β, hadm, hq⟩, _⟩ := isLeast_lamLo hwf hne
  rw [← hq]
  exact sepLo_le_Q hwf hadm

/-- **Theorem 1.4 (upper).** `λ₊ ≤ S₊` (when `Σ L > 0`). -/
theorem lamHi_le_sepHi (hwf : WellFormed s L U r R) (hL : 0 < ∑ k ∈ s, L k) :
    lamHi s L U r R ≤ sepHi s L U R := by
  obtain ⟨k, hks, hLk⟩ := exists_pos_of_sum_pos hL
  obtain ⟨⟨α, β, hadm, hq⟩, _⟩ :=
    isGreatest_lamHi hwf ⟨k, hks, hLk.trans_le (hwf k hks).2.1⟩
  rw [← hq]
  exact Q_le_sepHi hwf hL hadm

/-- **Theorem 1.4 (lower).** `λ₋ ≥ min { r k : U k > 0 }`. -/
theorem minRatio_le_lamLo (hwf : WellFormed s L U r R)
    (hne : (s.filter (fun k => 0 < U k)).Nonempty) :
    (s.filter (fun k => 0 < U k)).inf' hne r ≤ lamLo s L U r R := by
  obtain ⟨k, hk⟩ := hne
  obtain ⟨⟨α, β, hadm, hq⟩, _⟩ :=
    isLeast_lamLo hwf ⟨k, (mem_filter.mp hk).1, (mem_filter.mp hk).2⟩
  rw [← hq]
  refine le_Q_of_le_ratios ?_ hadm
  intro j hj hU
  exact inf'_le _ (mem_filter.mpr ⟨hj, hU⟩)

/-- **Theorem 1.4 (upper).** `λ₊ ≤ max { R k : U k > 0 }`. -/
theorem lamHi_le_maxRatio (hwf : WellFormed s L U r R)
    (hne : (s.filter (fun k => 0 < U k)).Nonempty) :
    lamHi s L U r R ≤ (s.filter (fun k => 0 < U k)).sup' hne R := by
  obtain ⟨k, hk⟩ := hne
  obtain ⟨⟨α, β, hadm, hq⟩, _⟩ :=
    isGreatest_lamHi hwf ⟨k, (mem_filter.mp hk).1, (mem_filter.mp hk).2⟩
  rw [← hq]
  refine Q_le_of_ratios_le ?_ hadm
  intro j hj hU
  exact le_sup' R (mem_filter.mpr ⟨hj, hU⟩)

/-- **Theorem 1.4 (combined lower).** `λ₋ ≥ max(S₋, min{r_k : U_k > 0})`. -/
theorem coupled_dominates_lower (hwf : WellFormed s L U r R)
    (hne : (s.filter (fun k => 0 < U k)).Nonempty) :
    max (sepLo s L U r) ((s.filter (fun k => 0 < U k)).inf' hne r) ≤ lamLo s L U r R := by
  obtain ⟨k, hk⟩ := hne
  exact max_le (sepLo_le_lamLo hwf ⟨k, (mem_filter.mp hk).1, (mem_filter.mp hk).2⟩)
    (minRatio_le_lamLo hwf _)

/-- **Theorem 1.4 (combined upper).** `λ₊ ≤ min(S₊, max{R_k : U_k > 0})` when `Σ L > 0`. -/
theorem coupled_dominates_upper (hwf : WellFormed s L U r R) (hL : 0 < ∑ k ∈ s, L k)
    (hne : (s.filter (fun k => 0 < U k)).Nonempty) :
    lamHi s L U r R ≤ min (sepHi s L U R) ((s.filter (fun k => 0 < U k)).sup' hne R) :=
  le_min (lamHi_le_sepHi hwf hL) (lamHi_le_maxRatio hwf hne)

/-! ### Theorem 1.5: exactly when coupling helps -/

/-- **Theorem 1.5 (lower).** If `S₋ > 0` then `λ₋ = S₋ ↔ ∀ k, r_k (U_k - L_k) = 0`. -/
theorem lamLo_eq_sepLo_iff (hwf : WellFormed s L U r R) (hS : 0 < sepLo s L U r) :
    lamLo s L U r R = sepLo s L U r ↔ ∀ k ∈ s, r k * (U k - L k) = 0 := by
  have hU0 : 0 ≤ ∑ k ∈ s, U k := sum_nonneg fun k hk => (hwf k hk).1.trans (hwf k hk).2.1
  have hU : 0 < ∑ k ∈ s, U k := by
    rcases hU0.lt_or_eq with h | h
    · exact h
    · rw [sepLo, ← h, div_zero] at hS
      exact absurd hS (lt_irrefl 0)
  have hP : 0 < ∑ k ∈ s, r k * L k := by
    have := mul_pos hS hU
    rwa [sepLo, div_mul_cancel₀ _ hU.ne'] at this
  have hne : ∃ k ∈ s, 0 < U k := exists_pos_of_sum_pos hU
  constructor
  · intro heq
    obtain ⟨⟨α, β, ⟨hk, hpos⟩, hq⟩, _⟩ := isLeast_lamLo hwf hne
    rw [heq, Q, sepLo, div_eq_div_iff hpos.ne' hU.ne'] at hq
    have hX : ∑ k ∈ s, r k * L k ≤ ∑ k ∈ s, r k * α k := sum_le_sum fun k hks =>
      mul_le_mul_of_nonneg_left (hk k hks).2.2.1 (hwf k hks).2.2.1
    have hB : ∑ k ∈ s, r k * α k ≤ ∑ k ∈ s, β k := sum_le_sum fun k hks => (hk k hks).2.2.2.2.1
    have hA : ∑ k ∈ s, α k ≤ ∑ k ∈ s, U k := sum_le_sum fun k hks => (hk k hks).2.2.2.1
    -- `B V = P A` with `P ≤ B`, `A ≤ V` and `P, V > 0` forces `A = V` and `B = P`.
    have h1 : (∑ k ∈ s, r k * L k) * ∑ k ∈ s, U k ≤ (∑ k ∈ s, β k) * ∑ k ∈ s, U k :=
      mul_le_mul_of_nonneg_right (hX.trans hB) hU.le
    have h2 : (∑ k ∈ s, r k * L k) * ∑ k ∈ s, α k ≤ (∑ k ∈ s, r k * L k) * ∑ k ∈ s, U k :=
      mul_le_mul_of_nonneg_left hA hP.le
    have hAV : ∑ k ∈ s, α k = ∑ k ∈ s, U k :=
      le_antisymm hA (le_of_mul_le_mul_left (by linarith) hP)
    have hBP : ∑ k ∈ s, β k = ∑ k ∈ s, r k * L k :=
      le_antisymm (le_of_mul_le_mul_right (by linarith) hU) (hX.trans hB)
    have hαU : ∀ k ∈ s, α k = U k :=
      (sum_eq_sum_iff_of_le fun k hks => (hk k hks).2.2.2.1).mp hAV
    have hrα : ∀ k ∈ s, r k * L k = r k * α k :=
      (sum_eq_sum_iff_of_le fun k hks =>
        mul_le_mul_of_nonneg_left (hk k hks).2.2.1 (hwf k hks).2.2.1).mp (by linarith)
    intro k hks
    rw [mul_sub, ← hαU k hks, ← hrα k hks, sub_self]
  · intro hz
    have hadm : IsAdm s L U r R U (fun k => r k * U k) := by
      refine ⟨fun k hk => ?_, hU⟩
      obtain ⟨h0, h1, h2, h3⟩ := hwf k hk
      have hU0 : 0 ≤ U k := h0.trans h1
      exact ⟨hU0, mul_nonneg h2 hU0, h1, le_rfl, le_rfl, mul_le_mul_of_nonneg_right h3 hU0⟩
    have hQ : Q s U (fun k => r k * U k) = sepLo s L U r := by
      rw [Q, sepLo]
      congr 1
      exact sum_congr rfl fun k hk => by linarith [hz k hk, mul_sub (r k) (U k) (L k)]
    exact le_antisymm (hQ ▸ (isLeast_lamLo hwf hne).2 ⟨U, _, hadm, rfl⟩)
      (sepLo_le_lamLo hwf hne)

/-- **Theorem 1.5 (upper).** If `Σ L > 0` then `λ₊ = S₊ ↔ ∀ k, R_k (U_k - L_k) = 0`. -/
theorem lamHi_eq_sepHi_iff (hwf : WellFormed s L U r R) (hL : 0 < ∑ k ∈ s, L k) :
    lamHi s L U r R = sepHi s L U R ↔ ∀ k ∈ s, R k * (U k - L k) = 0 := by
  have hne : ∃ k ∈ s, 0 < U k := by
    obtain ⟨k, hks, hLk⟩ := exists_pos_of_sum_pos hL
    exact ⟨k, hks, hLk.trans_le (hwf k hks).2.1⟩
  have hR0 : ∀ k ∈ s, 0 ≤ R k := fun k hk => (hwf k hk).2.2.1.trans (hwf k hk).2.2.2
  constructor
  · intro heq
    obtain ⟨⟨α, β, ⟨hk, hpos⟩, hq⟩, _⟩ := isGreatest_lamHi hwf hne
    rw [heq, Q, sepHi, div_eq_div_iff hpos.ne' hL.ne'] at hq
    have hY : ∑ k ∈ s, β k ≤ ∑ k ∈ s, R k * α k := sum_le_sum fun k hks => (hk k hks).2.2.2.2.2
    have hW : ∑ k ∈ s, R k * α k ≤ ∑ k ∈ s, R k * U k := sum_le_sum fun k hks =>
      mul_le_mul_of_nonneg_left (hk k hks).2.2.2.1 (hR0 k hks)
    have hA : ∑ k ∈ s, L k ≤ ∑ k ∈ s, α k := sum_le_sum fun k hks => (hk k hks).2.2.1
    have hW0 : 0 ≤ ∑ k ∈ s, R k * U k := sum_nonneg fun k hks =>
      mul_nonneg (hR0 k hks) ((hwf k hks).1.trans (hwf k hks).2.1)
    -- `B Λ = W A` with `B ≤ W`, `Λ ≤ A`, `Λ > 0`, `W ≥ 0`.
    have h1 : (∑ k ∈ s, β k) * ∑ k ∈ s, L k ≤ (∑ k ∈ s, R k * U k) * ∑ k ∈ s, L k :=
      mul_le_mul_of_nonneg_right (hY.trans hW) hL.le
    have h2 : (∑ k ∈ s, R k * U k) * ∑ k ∈ s, L k ≤ (∑ k ∈ s, R k * U k) * ∑ k ∈ s, α k :=
      mul_le_mul_of_nonneg_left hA hW0
    have hBW : ∑ k ∈ s, β k = ∑ k ∈ s, R k * U k :=
      le_antisymm (hY.trans hW) (le_of_mul_le_mul_right (by linarith) hL)
    have hRα : ∀ k ∈ s, R k * α k = R k * U k :=
      (sum_eq_sum_iff_of_le fun k hks =>
        mul_le_mul_of_nonneg_left (hk k hks).2.2.2.1 (hR0 k hks)).mp (by linarith)
    intro k hks
    rcases hW0.lt_or_eq with hWpos | hWzero
    · have hAL : ∑ k ∈ s, L k = ∑ k ∈ s, α k :=
        le_antisymm hA (le_of_mul_le_mul_left (by linarith) hWpos)
      have hαL : ∀ k ∈ s, L k = α k := (sum_eq_sum_iff_of_le fun k hks => (hk k hks).2.2.1).mp hAL
      rw [mul_sub, ← hRα k hks, hαL k hks, sub_self]
    · have hz := (sum_eq_zero_iff_of_nonneg fun k hks =>
        mul_nonneg (hR0 k hks) ((hwf k hks).1.trans (hwf k hks).2.1)).mp hWzero.symm k hks
      have hRL : R k * L k = 0 := le_antisymm
        (hz ▸ mul_le_mul_of_nonneg_left (hwf k hks).2.1 (hR0 k hks))
        (mul_nonneg (hR0 k hks) (hwf k hks).1)
      rw [mul_sub, hz, hRL, sub_self]
  · intro hz
    have hadm : IsAdm s L U r R L (fun k => R k * L k) := by
      refine ⟨fun k hk => ?_, hL⟩
      obtain ⟨h0, h1, h2, h3⟩ := hwf k hk
      exact ⟨h0, mul_nonneg (h2.trans h3) h0, le_rfl, h1, mul_le_mul_of_nonneg_right h3 h0, le_rfl⟩
    have hQ : Q s L (fun k => R k * L k) = sepHi s L U R := by
      rw [Q, sepHi]
      congr 1
      exact sum_congr rfl fun k hk => by linarith [hz k hk, mul_sub (R k) (U k) (L k)]
    exact le_antisymm (lamHi_le_sepHi hwf hL)
      (hQ ▸ (isGreatest_lamHi hwf hne).2 ⟨L, _, hadm, rfl⟩)

/-! ### Theorem 1.6: uniform-ratio exactness -/

/-- If `r k = R k = ρ` on every region with `U k > 0`, every admissible `Q` equals `ρ`. -/
theorem Q_eq_of_uniform {ρ : ℝ} (hρ : ∀ k ∈ s, 0 < U k → r k = ρ ∧ R k = ρ) {α β : ι → ℝ}
    (h : IsAdm s L U r R α β) : Q s α β = ρ := by
  obtain ⟨hk, hpos⟩ := h
  rw [Q, div_eq_iff hpos.ne', mul_sum]
  refine sum_congr rfl fun k hks => ?_
  obtain ⟨hα0, _, _, hαU, hrβ, hβR⟩ := hk k hks
  by_cases hU : 0 < U k
  · obtain ⟨hr, hR⟩ := hρ k hks hU
    rw [hr] at hrβ
    rw [hR] at hβR
    linarith
  · have hα : α k = 0 := le_antisymm (hαU.trans (not_lt.mp hU)) hα0
    rw [hα, mul_zero] at hrβ hβR ⊢
    linarith

/-- **Theorem 1.6.** Uniform regional ratios make the enclosure exact: `λ₋ = λ₊ = ρ`. -/
theorem lamLo_lamHi_of_uniform (hwf : WellFormed s L U r R) (hne : ∃ k ∈ s, 0 < U k) {ρ : ℝ}
    (hρ : ∀ k ∈ s, 0 < U k → r k = ρ ∧ R k = ρ) :
    lamLo s L U r R = ρ ∧ lamHi s L U r R = ρ := by
  obtain ⟨⟨α, β, hadm, hq⟩, _⟩ := isLeast_lamLo hwf hne
  obtain ⟨⟨α', β', hadm', hq'⟩, _⟩ := isGreatest_lamHi hwf hne
  exact ⟨hq ▸ Q_eq_of_uniform hρ hadm, hq' ▸ Q_eq_of_uniform hρ hadm'⟩

/-! ### Theorem 1.7: monotone refinement -/

section Refinement

variable {L' U' r' R' : ι → ℝ}

/-- Shrinking the intervals shrinks the admissible set. -/
theorem isAdm_mono (hL : ∀ k ∈ s, L k ≤ L' k) (hU : ∀ k ∈ s, U' k ≤ U k)
    (hr : ∀ k ∈ s, r k ≤ r' k) (hR : ∀ k ∈ s, R' k ≤ R k) {α β : ι → ℝ}
    (h : IsAdm s L' U' r' R' α β) : IsAdm s L U r R α β := by
  obtain ⟨hk, hpos⟩ := h
  refine ⟨fun k hks => ?_, hpos⟩
  obtain ⟨h0, h1, h2, h3, h4, h5⟩ := hk k hks
  exact ⟨h0, h1, (hL k hks).trans h2, h3.trans (hU k hks),
    (mul_le_mul_of_nonneg_right (hr k hks) h0).trans h4,
    h5.trans (mul_le_mul_of_nonneg_right (hR k hks) h0)⟩

/-- **Theorem 1.7(a), set form.** -/
theorem QSet_mono (hL : ∀ k ∈ s, L k ≤ L' k) (hU : ∀ k ∈ s, U' k ≤ U k)
    (hr : ∀ k ∈ s, r k ≤ r' k) (hR : ∀ k ∈ s, R' k ≤ R k) :
    QSet s L' U' r' R' ⊆ QSet s L U r R := by
  rintro q ⟨α, β, h, rfl⟩
  exact ⟨α, β, isAdm_mono hL hU hr hR h, rfl⟩

/-- **Theorem 1.7(a), lower endpoint.** Shrinking intervals can only raise `λ₋`.  The refined
data need only have a nonempty admissible set. -/
theorem lamLo_mono (hwf : WellFormed s L U r R) (hne' : (QSet s L' U' r' R').Nonempty)
    (hL : ∀ k ∈ s, L k ≤ L' k) (hU : ∀ k ∈ s, U' k ≤ U k)
    (hr : ∀ k ∈ s, r k ≤ r' k) (hR : ∀ k ∈ s, R' k ≤ R k) :
    lamLo s L U r R ≤ lamLo s L' U' r' R' := by
  have hsub := QSet_mono hL hU hr hR
  obtain ⟨q, hq⟩ := hne'
  obtain ⟨α, β, hadm, -⟩ := hsub hq
  have hne := (adm_nonempty_iff hwf).mp ⟨α, β, hadm⟩
  exact csInf_le_csInf (isLeast_lamLo hwf hne).bddBelow ⟨q, hq⟩ hsub

/-- **Theorem 1.7(a), upper endpoint.** Shrinking intervals can only lower `λ₊`. -/
theorem lamHi_anti (hwf : WellFormed s L U r R) (hne' : (QSet s L' U' r' R').Nonempty)
    (hL : ∀ k ∈ s, L k ≤ L' k) (hU : ∀ k ∈ s, U' k ≤ U k)
    (hr : ∀ k ∈ s, r k ≤ r' k) (hR : ∀ k ∈ s, R' k ≤ R k) :
    lamHi s L' U' r' R' ≤ lamHi s L U r R := by
  have hsub := QSet_mono hL hU hr hR
  obtain ⟨q, hq⟩ := hne'
  obtain ⟨α, β, hadm, -⟩ := hsub hq
  have hne := (adm_nonempty_iff hwf).mp ⟨α, β, hadm⟩
  exact csSup_le_csSup (isGreatest_lamHi hwf hne).bddAbove ⟨q, hq⟩ hsub

/-- **Theorem 1.7(a), certificate form.** A lower certificate for coarse data stays valid for
refined data: `g₋` can only increase when intervals shrink (refined masses nonnegative and
`L' ≤ U'`). -/
theorem gLo_mono (hL'0 : ∀ k ∈ s, 0 ≤ L' k) (hLU' : ∀ k ∈ s, L' k ≤ U' k)
    (hL : ∀ k ∈ s, L k ≤ L' k) (hU : ∀ k ∈ s, U' k ≤ U k) (hr : ∀ k ∈ s, r k ≤ r' k)
    (t : ℝ) : gLo s L U r t ≤ gLo s L' U' r' t := by
  refine sum_le_sum fun k hks => le_min ?_ ?_
  · calc min ((r k - t) * L k) ((r k - t) * U k) ≤ (r k - t) * L' k :=
          min_mul_le_mul (hL k hks) ((hLU' k hks).trans (hU k hks))
      _ ≤ (r' k - t) * L' k :=
          mul_le_mul_of_nonneg_right (by linarith [hr k hks]) (hL'0 k hks)
  · calc min ((r k - t) * L k) ((r k - t) * U k) ≤ (r k - t) * U' k :=
          min_mul_le_mul ((hL k hks).trans (hLU' k hks)) (hU k hks)
      _ ≤ (r' k - t) * U' k :=
          mul_le_mul_of_nonneg_right (by linarith [hr k hks]) ((hL'0 k hks).trans (hLU' k hks))

/-- **Theorem 1.7(a), certificate form (upper).** -/
theorem gHi_anti (hL'0 : ∀ k ∈ s, 0 ≤ L' k) (hLU' : ∀ k ∈ s, L' k ≤ U' k)
    (hL : ∀ k ∈ s, L k ≤ L' k) (hU : ∀ k ∈ s, U' k ≤ U k) (hR : ∀ k ∈ s, R' k ≤ R k)
    (t : ℝ) : gHi s L' U' R' t ≤ gHi s L U R t := by
  refine sum_le_sum fun k hks => max_le ?_ ?_
  · calc (R' k - t) * L' k ≤ (R k - t) * L' k :=
          mul_le_mul_of_nonneg_right (by linarith [hR k hks]) (hL'0 k hks)
      _ ≤ max ((R k - t) * L k) ((R k - t) * U k) :=
          mul_le_max_mul (hL k hks) ((hLU' k hks).trans (hU k hks))
  · calc (R' k - t) * U' k ≤ (R k - t) * U' k :=
          mul_le_mul_of_nonneg_right (by linarith [hR k hks]) ((hL'0 k hks).trans (hLU' k hks))
      _ ≤ max ((R k - t) * L k) ((R k - t) * U k) :=
          mul_le_max_mul ((hL k hks).trans (hLU' k hks)) (hU k hks)

end Refinement

section Splitting

variable [DecidableEq ι] {ι' : Type*} {s' : Finset ι'} {L' U' r' R' : ι' → ℝ}

/-- **Theorem 1.7(b), set form (refinement by splitting).** Fine regions `j ∈ s'` are mapped
to coarse regions `π j ∈ s`.  If the fine mass bounds of each fibre are consistent with the
coarse ones (`L_k ≤ Σ_{π j = k} L'_j`, `Σ_{π j = k} U'_j ≤ U_k`) and every fine ratio
interval lies inside its coarse one, then every refined value of `Q` is a coarse value. -/
theorem QSet_refine (π : ι' → ι) (hπ : ∀ j ∈ s', π j ∈ s)
    (hL : ∀ k ∈ s, L k ≤ ∑ j ∈ s' with π j = k, L' j)
    (hU : ∀ k ∈ s, ∑ j ∈ s' with π j = k, U' j ≤ U k)
    (hr : ∀ j ∈ s', r (π j) ≤ r' j) (hR : ∀ j ∈ s', R' j ≤ R (π j)) :
    QSet s' L' U' r' R' ⊆ QSet s L U r R := by
  rintro q ⟨α', β', ⟨hk, hpos⟩, rfl⟩
  set α : ι → ℝ := fun k => ∑ j ∈ s' with π j = k, α' j with hα
  set β : ι → ℝ := fun k => ∑ j ∈ s' with π j = k, β' j with hβ
  have hsumα : ∑ k ∈ s, α k = ∑ j ∈ s', α' j := sum_fiberwise_of_maps_to hπ α'
  have hsumβ : ∑ k ∈ s, β k = ∑ j ∈ s', β' j := sum_fiberwise_of_maps_to hπ β'
  refine ⟨α, β, ⟨fun k hks => ?_, by rw [hsumα]; exact hpos⟩, ?_⟩
  · have hm : ∀ j ∈ s'.filter (fun j => π j = k), j ∈ s' ∧ π j = k := fun j hj =>
      mem_filter.mp hj
    refine ⟨sum_nonneg fun j hj => (hk j (hm j hj).1).1,
      sum_nonneg fun j hj => (hk j (hm j hj).1).2.1,
      (hL k hks).trans (sum_le_sum fun j hj => (hk j (hm j hj).1).2.2.1),
      (sum_le_sum fun j hj => (hk j (hm j hj).1).2.2.2.1).trans (hU k hks), ?_, ?_⟩
    · rw [hα, hβ, mul_sum]
      refine sum_le_sum fun j hj => ?_
      obtain ⟨hj1, hj2⟩ := hm j hj
      have hrj := hr j hj1
      rw [hj2] at hrj
      exact (mul_le_mul_of_nonneg_right hrj (hk j hj1).1).trans (hk j hj1).2.2.2.2.1
    · rw [hα, hβ, mul_sum]
      refine sum_le_sum fun j hj => ?_
      obtain ⟨hj1, hj2⟩ := hm j hj
      have hRj := hR j hj1
      rw [hj2] at hRj
      exact (hk j hj1).2.2.2.2.2.trans (mul_le_mul_of_nonneg_right hRj (hk j hj1).1)
  · rw [Q, Q, hsumα, hsumβ]

/-- **Theorem 1.7(b), lower endpoint.** Splitting regions can only raise `λ₋`. -/
theorem lamLo_refine (hwf : WellFormed s L U r R) (hne' : (QSet s' L' U' r' R').Nonempty)
    (π : ι' → ι) (hπ : ∀ j ∈ s', π j ∈ s)
    (hL : ∀ k ∈ s, L k ≤ ∑ j ∈ s' with π j = k, L' j)
    (hU : ∀ k ∈ s, ∑ j ∈ s' with π j = k, U' j ≤ U k)
    (hr : ∀ j ∈ s', r (π j) ≤ r' j) (hR : ∀ j ∈ s', R' j ≤ R (π j)) :
    lamLo s L U r R ≤ lamLo s' L' U' r' R' := by
  have hsub := QSet_refine π hπ hL hU hr hR
  obtain ⟨q, hq⟩ := hne'
  obtain ⟨α, β, hadm, -⟩ := hsub hq
  have hne := (adm_nonempty_iff hwf).mp ⟨α, β, hadm⟩
  exact csInf_le_csInf (isLeast_lamLo hwf hne).bddBelow ⟨q, hq⟩ hsub

/-- **Theorem 1.7(b), upper endpoint.** Splitting regions can only lower `λ₊`. -/
theorem lamHi_refine (hwf : WellFormed s L U r R) (hne' : (QSet s' L' U' r' R').Nonempty)
    (π : ι' → ι) (hπ : ∀ j ∈ s', π j ∈ s)
    (hL : ∀ k ∈ s, L k ≤ ∑ j ∈ s' with π j = k, L' j)
    (hU : ∀ k ∈ s, ∑ j ∈ s' with π j = k, U' j ≤ U k)
    (hr : ∀ j ∈ s', r (π j) ≤ r' j) (hR : ∀ j ∈ s', R' j ≤ R (π j)) :
    lamHi s' L' U' r' R' ≤ lamHi s L U r R := by
  have hsub := QSet_refine π hπ hL hU hr hR
  obtain ⟨q, hq⟩ := hne'
  obtain ⟨α, β, hadm, -⟩ := hsub hq
  have hne := (adm_nonempty_iff hwf).mp ⟨α, β, hadm⟩
  exact csSup_le_csSup (isGreatest_lamHi hwf hne).bddAbove ⟨q, hq⟩ hsub

end Splitting

end CoupledRatio
end FoldBounds
