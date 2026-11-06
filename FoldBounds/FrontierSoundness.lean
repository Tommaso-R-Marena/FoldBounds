module

public import FoldBounds.FrontierCheck

/-!
# End-to-end soundness of the frontier checker (flagship theorem)

`check_sound`: if `FrontierCheck.check` accepts an untrusted certificate
for the H/P sequence `seq`, rational `q` and threshold `thr`, then the
returned intervals enclose

* the exact partition function `Z = ∑_c q^{k(c)}`,
* the exact event mass `Z_A` of `A = {c | thr ≤ k(c)}`,
* the exact complement mass `Z_{¬A}`, and
* the exact equilibrium probability `P(A) = Z_A / Z`,

all stated with the *pre-existing* definitions `Executable.massQ`,
`LatticeModel.weightQ` and `LatticeModel.contactEvent` of the Lean HP
model over the full conformation type `LatticeModel.Conformation steps`.
No hypothesis about the certificate is assumed; `q ≥ 1` and the sequence
length are checked by `check` itself, and nonemptiness of the state space
is proved, not assumed.
-/
@[expose] public section

namespace FoldBounds
namespace FrontierCheck

open LatticeModel Frontier Finset

theorem massQ_univ_eq (steps : ℕ) (q : ℚ) (seq : Fin (steps + 2) → Bool) :
    Executable.massQ (weightQ steps q seq) univ = massZ steps (List.ofFn seq) q [] := by
  unfold Executable.massQ weightQ massZ
  rw [ext_nil]
  exact sum_conformation_eq_sum_validWords steps (fun l => q ^ wordContacts (List.ofFn seq) l) _

theorem massQ_event_eq (steps : ℕ) (q : ℚ) (seq : Fin (steps + 2) → Bool) (thr : ℕ) :
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) =
      massA steps (List.ofFn seq) q thr [] := by
  unfold Executable.massQ weightQ massA contactEvent
  rw [ext_nil, Finset.sum_filter]
  exact sum_conformation_eq_sum_validWords steps (fun l =>
    if thr ≤ wordContacts (List.ofFn seq) l then q ^ wordContacts (List.ofFn seq) l else 0) _

theorem massQ_compl_eq (steps : ℕ) (q : ℚ) (seq : Fin (steps + 2) → Bool) (thr : ℕ) :
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ =
      massN steps (List.ofFn seq) q thr [] := by
  unfold Executable.massQ weightQ massN contactEvent
  rw [ext_nil, Finset.compl_filter, Finset.sum_filter]
  have := sum_conformation_eq_sum_validWords steps (fun l =>
    if wordContacts (List.ofFn seq) l < thr then q ^ wordContacts (List.ofFn seq) l else 0)
    Finset.univ
  simp only [not_le]
  exact this

theorem massZ_eq_add (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) :
    massZ steps S q p = massA steps S q thr p + massN steps S q thr p := by
  unfold massZ massA massN
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun l _ => ?_)
  split_ifs with h1 h2 <;> first | omega | simp

theorem massZ_root_pos (steps : ℕ) (S : List Bool) {q : ℚ} (hq : 0 < q) :
    0 < massZ steps S q [] := by
  unfold massZ
  rw [ext_nil]
  have hmem : List.replicate steps (0 : Fin 4) ∈ validWords steps :=
    mem_validWords.2 ⟨by simp, walkFrom_straight_nodup steps⟩
  exact Finset.sum_pos' (fun l _ => (pow_pos hq _).le) ⟨_, hmem, pow_pos hq _⟩

/-- Rational form of Aristotle's coupled probability enclosure. -/
theorem ratio_mem_coupled_rat {a b aL aU bL bU : ℚ} (haL : 0 ≤ aL) (hbL : 0 ≤ bL)
    (ha : a ∈ Set.Icc aL aU) (hb : b ∈ Set.Icc bL bU) (hab : 0 < a + b) :
    a / (a + b) ∈ Set.Icc (aL / (aL + bU)) (aU / (aU + bL)) := by
  have h := FoldBounds.ratio_mem_coupled (a := (a : ℝ)) (b := (b : ℝ)) (aL := (aL : ℝ))
    (aU := (aU : ℝ)) (bL := (bL : ℝ)) (bU := (bU : ℝ)) (by exact_mod_cast haL)
    (by exact_mod_cast hbL) ⟨by exact_mod_cast ha.1, by exact_mod_cast ha.2⟩
    ⟨by exact_mod_cast hb.1, by exact_mod_cast hb.2⟩ (by exact_mod_cast hab)
  exact ⟨by exact_mod_cast h.1, by exact_mod_cast h.2⟩

/-- Acceptance exposes the internal root enclosure. -/
theorem check_eq_some {S : List Bool} {q : ℚ} {thr : ℕ} {t : FTree} {r : FrontierResult}
    (h : check S q thr t = some r) :
    2 ≤ S.length ∧ 1 ≤ q ∧ ∃ e, checkNode (S.length - 2) S q thr [] t = some e ∧
      r = finalize e := by
  unfold check at h
  split_ifs at h with hc
  obtain ⟨e, he, rfl⟩ := Option.map_eq_some_iff.1 h
  exact ⟨hc.1, hc.2, e, he, rfl⟩

/-- **Flagship theorem (end-to-end soundness of the compact frontier
checker).** For every H/P sequence of length `steps + 2`, every rational
`q`, every contact threshold `thr` and every untrusted certificate `t`: if
`check` accepts and returns `r`, then `r` encloses the exact partition
function, event mass, complement mass and equilibrium probability
`P(k ≥ thr)` of the Lean HP model. -/
theorem check_sound (steps : ℕ) (seq : Fin (steps + 2) → Bool) (q : ℚ) (thr : ℕ)
    (t : FTree) (r : FrontierResult) (h : check (List.ofFn seq) q thr t = some r) :
    Executable.massQ (weightQ steps q seq) univ ∈ Set.Icc r.zLo r.zHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) ∈
      Set.Icc r.eLo r.eHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ ∈
      Set.Icc r.nLo r.nHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) /
      Executable.massQ (weightQ steps q seq) univ ∈ Set.Icc r.pLo r.pHi := by
  obtain ⟨_, hq, e, he, rfl⟩ := check_eq_some h
  have hlen : (List.ofFn seq).length - 2 = steps := by simp
  rw [hlen] at he
  have hS : (List.ofFn seq).length = steps + 2 := by simp
  obtain ⟨z1, z2, a1, a2, n1, n2, _, a0, n0⟩ :=
    checkNode_sound steps (List.ofFn seq) hS q hq thr t [] e he
  rw [massQ_univ_eq, massQ_event_eq, massQ_compl_eq]
  refine ⟨⟨z1, z2⟩, ⟨a1, a2⟩, ⟨n1, n2⟩, ?_⟩
  rw [massZ_eq_add steps (List.ofFn seq) q thr []]
  have hpos := massZ_root_pos steps (List.ofFn seq) (q := q) (by linarith)
  rw [massZ_eq_add steps (List.ofFn seq) q thr []] at hpos
  exact ratio_mem_coupled_rat a0 n0 ⟨a1, a2⟩ ⟨n1, n2⟩ hpos

/-- The same probability enclosure, stated over the reals exactly in the
form of the existing `checkHP_coupled_sound`. -/
theorem check_coupled_probability_sound (steps : ℕ) (seq : Fin (steps + 2) → Bool)
    (q : ℚ) (thr : ℕ) (t : FTree) (r : FrontierResult)
    (h : check (List.ofFn seq) q thr t = some r) :
    ((Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) : ℝ) /
      ((Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) : ℝ) +
       (Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ : ℝ))) ∈
      Set.Icc (r.pLo : ℝ) (r.pHi : ℝ) := by
  have hsound := check_sound steps seq q thr t r h
  have hsplit : Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) +
      Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ =
      Executable.massQ (weightQ steps q seq) univ :=
    Finset.sum_add_sum_compl _ _
  have hp := hsound.2.2.2
  rw [← hsplit] at hp
  exact ⟨by exact_mod_cast hp.1, by exact_mod_cast hp.2⟩

/-- **Compactness witness.** A one-node certificate (a single `bound` leaf
at the root) is accepted for *every* admissible sequence and `q ≥ 1`:
acceptance never requires enumerating the conformation space. -/
theorem check_root_bound_accepted (S : List Bool) (q : ℚ) (thr : ℕ) (tr : Bool)
    (hS : 2 ≤ S.length) (hq : 1 ≤ q) :
    (check S q thr (.bound tr [])).isSome := by
  have hn : (walkFrom ([] : List (Fin 4))).Nodup := by decide
  have hw : witnessesOK (S.length - 2) [] [] := ⟨List.nodup_nil, by simp⟩
  simp [check, hS, hq, checkNode, hn, hw]

theorem knownContacts_nil (S : List Bool) : knownContacts S [] = 0 := by
  unfold knownContacts parityContactCount
  refine Finset.sum_eq_zero (fun i _ => Finset.sum_eq_zero (fun j _ => ?_))
  have : (walkFrom ([] : List (Fin 4))).length = 2 := rfl
  rw [this]
  by_cases hc : candidatePair S i j = true
  · have := candidatePair_lt hc
    have hj : ¬ j < 2 := by omega
    simp [hj]
  · simp [hc]

/-- The basic root-bound enclosure: `Z ≤ 3^(n-2) q^{K_max}` without any
enumeration (`K_max` = number of parity-compatible H–H candidate pairs). -/
theorem Z_le_basic_root (steps : ℕ) (seq : Fin (steps + 2) → Bool) {q : ℚ} (hq : 1 ≤ q) :
    Executable.massQ (weightQ steps q seq) univ ≤
      3 ^ steps * q ^ futurePairs (List.ofFn seq) 2 := by
  have hS : (List.ofFn seq).length = steps + 2 := by simp
  have hn : (walkFrom ([] : List (Fin 4))).Nodup := by decide
  have hw : witnessesOK steps [] [] := ⟨List.nodup_nil, by simp⟩
  have h := (boundEnc_sound false [] steps (List.ofFn seq) hS q hq 0 [] (by simp) hn hw).2.1
  rw [massQ_univ_eq]
  refine h.trans (le_of_eq ?_)
  simp only [boundEnc, envelopeUppers, basicZ, kcap, kmax, List.map_nil, List.sum_nil, sub_zero,
    Bool.false_eq_true, if_false, List.length_nil, Nat.sub_zero, zero_add]
  have : knownContacts (List.ofFn seq) [] = 0 := knownContacts_nil _
  rw [this, zero_add]

end FrontierCheck
end FoldBounds
