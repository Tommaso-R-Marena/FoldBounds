module

public import FoldBounds.MotionUncertainty

/-!
# FoldBounds.MotionContactChecker

A finite, executable **rational** contact-membership gate. A proposed interval
for squared contact distance determines whether a spring is forced on,
forced off, or optional. A proposed interval for stiffness bounds the active
spring. All untrusted inequalities are checked with exact rational comparisons.

This module proves that accepting actual distances/stiffnesses implies
energy enclosures for every displacement field; no trust is placed in a
floating-point contact selector. It does not verify the upstream conversion
of Cartesian coordinates into squared-distance enclosures.
-/

@[expose] public section

namespace FoldBounds.MotionContactChecker

open FoldBounds.MotionUncertainty

/-- Untrusted finite distance and stiffness intervals. -/
structure RawContacts (m : ℕ) where
  cutoffSq : ℚ
  d2Lower : Fin m → ℚ
  d2Upper : Fin m → ℚ
  kLower : Fin m → ℚ
  kUpper : Fin m → ℚ

/-- Entirely decidable consistency gate over exact rational inputs. -/
def accepts {m : ℕ} (c : RawContacts m) : Bool :=
  decide (0 < c.cutoffSq ∧ ∀ i, 0 ≤ c.d2Lower i ∧
    c.d2Lower i ≤ c.d2Upper i ∧ 0 ≤ c.kLower i ∧
    c.kLower i ≤ c.kUpper i)

/-- Check that the actual contact-network witnesses lie inside the proposed
rational enclosures. No symbolic or real-valued oracle is assumed. -/
def admissible {m : ℕ} (c : RawContacts m)
    (d2 k : Fin m → ℚ) : Prop :=
  ∀ i, c.d2Lower i ≤ d2 i ∧ d2 i ≤ c.d2Upper i ∧
    c.kLower i ≤ k i ∧ k i ≤ c.kUpper i

/-- A guaranteed contact contributes at least its minimum spring stiffness. -/
def lowerWeight {m : ℕ} (c : RawContacts m) (i : Fin m) : ℚ :=
  if c.d2Upper i ≤ c.cutoffSq then c.kLower i else 0

/-- A guaranteed non-contact contributes no stiffness. Every other candidate
is upper-bounded by the maximum active stiffness. -/
def upperWeight {m : ℕ} (c : RawContacts m) (i : Fin m) : ℚ :=
  if c.cutoffSq < c.d2Lower i then 0 else c.kUpper i

/-- Physically selected contact weight for actual squared distance and
stiffness: edges outside the cutoff are absent. -/
def actualWeight {m : ℕ} (c : RawContacts m)
    (d2 k : Fin m → ℚ) (i : Fin m) : ℚ :=
  if d2 i ≤ c.cutoffSq then k i else 0

/-- Execute both the interval-consistency and concrete-witness gates. -/
def check {m : ℕ} (c : RawContacts m)
    (d2 k : Fin m → ℚ) : Bool :=
  accepts c && decide (∀ i, c.d2Lower i ≤ d2 i ∧
    d2 i ≤ c.d2Upper i ∧ c.kLower i ≤ k i ∧ k i ≤ c.kUpper i)

theorem accepts_valid {m : ℕ} (c : RawContacts m)
    (h : accepts c = true) :
    0 < c.cutoffSq ∧ ∀ i, 0 ≤ c.d2Lower i ∧
      c.d2Lower i ≤ c.d2Upper i ∧ 0 ≤ c.kLower i ∧
      c.kLower i ≤ c.kUpper i := by
  exact of_decide_eq_true h

/-- Endpoint bounds for any consistent actual contact realization, including
networks whose contact topology changes at the cutoff. -/
theorem weight_interval_sound {m : ℕ} (c : RawContacts m)
    (d2 k : Fin m → ℚ) (hacc : accepts c = true)
    (hdata : admissible c d2 k) (i : Fin m) :
    lowerWeight c i ≤ actualWeight c d2 k i ∧
      actualWeight c d2 k i ≤ upperWeight c i := by
  have hv := (accepts_valid c hacc).2 i
  have hd := hdata i
  constructor
  · by_cases hAlways : c.d2Upper i ≤ c.cutoffSq
    · have hOn : d2 i ≤ c.cutoffSq := le_trans hd.2.1 hAlways
      simpa [lowerWeight, actualWeight, hAlways, hOn] using hd.2.2.1
    · simp only [lowerWeight, if_neg hAlways]
      by_cases hOn : d2 i ≤ c.cutoffSq
      · simp only [actualWeight, if_pos hOn]
        linarith
      · simp only [actualWeight, if_neg hOn]
        exact le_refl _
  · by_cases hNever : c.cutoffSq < c.d2Lower i
    · have hOff : ¬ d2 i ≤ c.cutoffSq := by linarith [hd.1]
      simp [upperWeight, actualWeight, hNever, hOff]
    · simp only [upperWeight, if_neg hNever]
      by_cases hOn : d2 i ≤ c.cutoffSq
      · simpa [actualWeight, hOn] using hd.2.2.2
      · simp only [actualWeight, if_neg hOn]
        linarith

theorem actual_weight_nonneg {m : ℕ} (c : RawContacts m)
    (d2 k : Fin m → ℚ) (hacc : accepts c = true)
    (hdata : admissible c d2 k) (i : Fin m) :
    0 ≤ actualWeight c d2 k i := by
  have hv := (accepts_valid c hacc).2 i
  have hd := hdata i
  by_cases hOn : d2 i ≤ c.cutoffSq
  · simpa [actualWeight, hOn] using (le_trans hv.2.2.1 hd.2.2.1)
  · simp [actualWeight, hOn]

/-- A valid raw interval automatically produces a coherent endpoint envelope. -/
theorem lower_le_upper {m : ℕ} (c : RawContacts m)
    (hacc : accepts c = true) (i : Fin m) :
    lowerWeight c i ≤ upperWeight c i := by
  have hv := (accepts_valid c hacc).2 i
  by_cases hAlways : c.d2Upper i ≤ c.cutoffSq
  · have hNotOff : ¬ c.cutoffSq < c.d2Lower i := by linarith [hv.2.1]
    simp only [lowerWeight, if_pos hAlways, upperWeight, if_neg hNotOff]
    exact hv.2.2.2
  · simp only [lowerWeight, if_neg hAlways]
    by_cases hOff : c.cutoffSq < c.d2Lower i
    · simp [upperWeight, hOff]
    · simp [upperWeight, hOff]
      linarith

/-- Extract the proved validity of every interval and witness from the
actual executable Boolean gate, without assuming proposer honesty. -/
theorem check_valid {m : ℕ} (c : RawContacts m) (d2 k : Fin m → ℚ)
    (h : check c d2 k = true) :
    accepts c = true ∧ admissible c d2 k := by
  have hc : accepts c = true := by
    cases haccept : accepts c with
    | false => simp [check, haccept] at h
    | true => exact haccept
  have hd : decide (∀ i, c.d2Lower i ≤ d2 i ∧
      d2 i ≤ c.d2Upper i ∧ c.kLower i ≤ k i ∧
      k i ≤ c.kUpper i) = true := by
    simpa [check, hc] using h
  exact ⟨hc, of_decide_eq_true hd⟩

/-- A successful executable check supplies every hypothesis needed for
all possible orientations/displacements of the contact supergraph. -/
theorem check_energy_sound {n m : ℕ}
    (edges : Fin m → Spring n) (c : RawContacts m)
    (d2 k : Fin m → ℚ) (u : Fin n → Fin 3 → ℝ)
    (h : check c d2 k = true) :
    springEnergy edges (fun i => (lowerWeight c i : ℝ)) u ≤
      springEnergy edges (fun i => (actualWeight c d2 k i : ℝ)) u ∧
    springEnergy edges (fun i => (actualWeight c d2 k i : ℝ)) u ≤
      springEnergy edges (fun i => (upperWeight c i : ℝ)) u := by
  have hh := check_valid c d2 k h
  have hbounds := weight_interval_sound c d2 k hh.1 hh.2
  constructor
  · apply springEnergy_mono
    intro i
    exact_mod_cast (hbounds i).1
  · apply springEnergy_mono
    intro i
    exact_mod_cast (hbounds i).2

/-- Universal counterpart: once the untrusted interval certificate passes
the gate, ANY compatible squared distances and spring strengths are enclosed.
The concrete witness does not need to be enumerated by the checker. -/
theorem all_admissible_energy_sound {n m : ℕ}
    (edges : Fin m → Spring n) (c : RawContacts m)
    (u : Fin n → Fin 3 → ℝ) (hacc : accepts c = true) :
    ∀ (d2 k : Fin m → ℚ), admissible c d2 k →
      springEnergy edges (fun i => (lowerWeight c i : ℝ)) u ≤
        springEnergy edges (fun i => (actualWeight c d2 k i : ℝ)) u ∧
      springEnergy edges (fun i => (actualWeight c d2 k i : ℝ)) u ≤
        springEnergy edges (fun i => (upperWeight c i : ℝ)) u := by
  intro d2 k hdata
  have hb := weight_interval_sound c d2 k hacc hdata
  constructor
  · apply springEnergy_mono
    intro i
    exact_mod_cast (hb i).1
  · apply springEnergy_mono
    intro i
    exact_mod_cast (hb i).2

end FoldBounds.MotionContactChecker
