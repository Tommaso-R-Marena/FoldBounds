module

public import FoldBounds.FreeEnergy

/-!
# FoldBounds.Symmetry — counting and symmetry-reduced certificates

If the energy factors through a map `f : C → D` (`E = e ∘ f`), then
`Z = ∑_d m_d · exp (-β e_d)` with multiplicities `m_d = #{c | f c = d}`
(density of states / combinatorial counting).  The prime instance is the orbit map of
an energy-preserving group action, where `D` is the orbit space and `m_d` the orbit size.

Consequence: a certified enclosure of `Z` only needs (i) multiplicities, (ii) exact
energies of the *explored representatives* `R ⊆ D`, and (iii) energy intervals for the
unexplored representatives.  The number of energy evaluations drops from `|C|` to `|D|`.
-/

@[expose] public section

open Finset Real

namespace FoldBounds

variable {C : Type*} [Fintype C] {D : Type*} [Fintype D] [DecidableEq D]

/-- Multiplicity of a reduced state `d`: the number of conformations mapped to `d`. -/
def multiplicity (f : C → D) (d : D) : ℕ := (Finset.univ.filter fun c => f c = d).card

variable (β : ℝ) (E : C → ℝ)

/-- **Counting formula.** If `E = e ∘ f` then `Z = ∑_d m_d · exp (-β e_d)` (any real `β`). -/
theorem Z_eq_sum_multiplicity (f : C → D) (e : D → ℝ) (hE : ∀ c, E c = e (f c)) :
    Z β E = ∑ d, (multiplicity f d : ℝ) * Real.exp (-β * e d) := by
  rw [← ZOn_univ, ← blockBound_eq_ZOn_of_exact E f β e _ (fun c _ => hE c),
    blockBound_eq_sum_card]
  rfl

/-- **Reduced certified enclosure.** Explored reduced states `R ⊆ D` are evaluated exactly;
every unexplored reduced state carries an energy interval `ℓ d ≤ e d ≤ u d`.  Then
`∑_{d∈R} m_d e^{-β e_d} + ∑_{d∉R} m_d e^{-β u_d} ≤ Z ≤ ∑_{d∈R} m_d e^{-β e_d} + ∑_{d∉R} m_d e^{-β ℓ_d}`. -/
theorem Z_mem_Icc_of_reduced (hβ : 0 ≤ β) (f : C → D) (e : D → ℝ)
    (hE : ∀ c, E c = e (f c)) (R : Finset D) (ℓ u : D → ℝ)
    (hℓ : ∀ d ∉ R, ℓ d ≤ e d) (hu : ∀ d ∉ R, e d ≤ u d) :
    Z β E ∈ Set.Icc
      (∑ d ∈ R, (multiplicity f d : ℝ) * Real.exp (-β * e d) +
        ∑ d ∈ Rᶜ, (multiplicity f d : ℝ) * Real.exp (-β * u d))
      (∑ d ∈ R, (multiplicity f d : ℝ) * Real.exp (-β * e d) +
        ∑ d ∈ Rᶜ, (multiplicity f d : ℝ) * Real.exp (-β * ℓ d)) := by
  rw [Z_eq_sum_multiplicity β E f e hE, ← Finset.sum_add_sum_compl R]
  constructor
  · refine add_le_add le_rfl (Finset.sum_le_sum fun d hd => ?_)
    have := hu d (Finset.mem_compl.mp hd)
    exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr (by nlinarith)) (Nat.cast_nonneg _)
  · refine add_le_add le_rfl (Finset.sum_le_sum fun d hd => ?_)
    have := hℓ d (Finset.mem_compl.mp hd)
    exact mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr (by nlinarith)) (Nat.cast_nonneg _)

section GroupAction

open MulAction

variable {G : Type*} [Group G] [MulAction G C]

omit [Fintype C] in
/-- Energy invariance under a group action makes the energy constant on orbits. -/
theorem energy_eq_of_mem_orbit (hinv : ∀ (g : G) c, E (g • c) = E c) {c c' : C}
    (h : c ∈ orbit G c') : E c = E c' := by
  obtain ⟨g, rfl⟩ := h
  exact hinv g c'

omit [Fintype C] in
/-- The energy of a conformation equals the energy of its orbit's chosen representative. -/
theorem energy_eq_energy_out (hinv : ∀ (g : G) c, E (g • c) = E c) (c : C) :
    E c = E (Quotient.out (Quotient.mk (orbitRel G C) c)) := by
  have h : Quotient.out (Quotient.mk (orbitRel G C) c) ∈ orbit G c :=
    (orbitRel_apply).mp (Quotient.mk_out (s := orbitRel G C) c)
  exact (energy_eq_of_mem_orbit E hinv h).symm

variable [Fintype (orbitRel.Quotient G C)] [DecidableEq (orbitRel.Quotient G C)]

/-- **Symmetry reduction.** For an energy-preserving action of `G` on `C`,
`Z = ∑_{orbits ω} |ω| · exp (-β E(rep ω))`, where `|ω|` is the number of conformations in
the orbit and `rep ω` any chosen representative. -/
theorem Z_eq_sum_orbits (hinv : ∀ (g : G) c, E (g • c) = E c) :
    Z β E = ∑ ω : orbitRel.Quotient G C,
      (multiplicity (Quotient.mk (orbitRel G C)) ω : ℝ) * Real.exp (-β * E (Quotient.out ω)) :=
  Z_eq_sum_multiplicity β E _ (fun ω => E (Quotient.out ω)) (energy_eq_energy_out E hinv)

omit [Fintype (orbitRel.Quotient G C)] in
open Classical in
/-- The multiplicity of an orbit in the reduction formula is its size. -/
theorem multiplicity_orbit_eq (ω : orbitRel.Quotient G C) :
    multiplicity (Quotient.mk (orbitRel G C)) ω =
      (Finset.univ.filter fun c => c ∈ orbit G (Quotient.out ω)).card := by
  unfold multiplicity
  congr 1
  refine Finset.filter_congr fun c _ => ?_
  conv_lhs => rw [← Quotient.out_eq ω]
  rw [Quotient.eq]
  exact orbitRel_apply

end GroupAction

end FoldBounds
