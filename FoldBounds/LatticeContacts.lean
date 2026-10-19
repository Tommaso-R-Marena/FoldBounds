module

public import FoldBounds.LatticeEquivalence

/-!
# Contact-scoring semantics: checkerboard parity and the Python scorer

The Lean reference scorer `LatticeModel.contactCount` sums over all index
pairs `i + 1 < j` of hydrophobic residues whose sites are grid adjacent.
The Python reference (`_possible_hh_pairs` / `_contacts` in
`src/foldbounds/lattice.py`) instead only inspects *candidate* pairs with
`j ≥ i + 3` and odd gap `j - i`, filtered by `j < len(path)`, and tests
adjacency through the Manhattan distance `|dx| + |dy| = 1`.

This file proves

* the checkerboard parity theorem: site `i` of every decoded walk has
  colour `i mod 2`, so lattice-adjacent sites have chain indices of
  opposite parity (`index_parity_of_adjacent`), and an even index gap
  never produces a contact (`not_adjacent_of_even_gap`);
* exact equality of the two scorers on every decoded walk whose length
  matches the sequence (`contactCount_eq_parityContactCount`);
* the rational weight facts `q ^ k` (equality and monotonicity for `q ≥ 1`).
-/
@[expose] public section

namespace FoldBounds
namespace LatticeModel

open Finset

/-- Python's adjacency test: Manhattan distance exactly one. -/
def manhattan1 (p r : Site) : Bool :=
  decide ((p.1 - r.1).natAbs + (p.2 - r.2).natAbs = 1)

theorem adjacent_eq_manhattan1 (p r : Site) : adjacent p r = manhattan1 p r := by
  obtain ⟨a, b⟩ := p
  obtain ⟨c, d⟩ := r
  simp only [adjacent, manhattan1]
  apply Bool.decide_congr
  omega

/-- Checkerboard colour of a site. -/
def colour (p : Site) : ℤ := (p.1 + p.2) % 2

theorem colour_advance (p : Site) (d : Fin 4) :
    colour (advance p d) = (colour p + 1) % 2 := by
  obtain ⟨a, b⟩ := p
  fin_cases d <;> simp [colour, advance] <;> omega

theorem colour_ne_of_adjacent {p r : Site} (h : adjacent p r = true) :
    colour p ≠ colour r := by
  obtain ⟨a, b⟩ := p
  obtain ⟨c, d⟩ := r
  simp only [adjacent, decide_eq_true_eq] at h
  simp only [colour]
  omega

theorem colour_trace (s : Site) (ds : List (Fin 4)) (k : ℕ)
    (hk : k < (trace s ds).length) :
    colour ((trace s ds)[k]) = (colour s + k + 1) % 2 := by
  induction ds generalizing s k with
  | nil => simp [trace] at hk
  | cons d ds ih =>
      cases k with
      | zero => simp [trace, colour_advance]
      | succ k =>
          simp only [trace, List.getElem_cons_succ]
          rw [ih (advance s d) k (by simpa [trace] using hk), colour_advance]
          push_cast
          omega

/-- **Checkerboard parity theorem.** Site `i` of a decoded anchored walk
has colour `i mod 2`. -/
theorem colour_walkFrom (ds : List (Fin 4)) (i : ℕ) (hi : i < (walkFrom ds).length) :
    colour ((walkFrom ds)[i]) = (i : ℤ) % 2 := by
  rcases i with _ | _ | k
  · simp [walkFrom, colour, origin]
  · simp [walkFrom, colour, first]
  · have hk : k < (trace first ds).length := by
      simp [walkFrom_length, trace_length] at hi ⊢; omega
    simp only [walkFrom, List.cons_append, List.nil_append, List.getElem_cons_succ]
    rw [colour_trace first ds k hk]
    simp [colour, first]
    omega

/-- Lattice-adjacent sites of a decoded walk have chain indices of
opposite parity. -/
theorem index_parity_of_adjacent (ds : List (Fin 4)) {i j : ℕ}
    (hi : i < (walkFrom ds).length) (hj : j < (walkFrom ds).length)
    (h : adjacent ((walkFrom ds)[i]) ((walkFrom ds)[j]) = true) :
    i % 2 ≠ j % 2 := by
  have h1 := colour_walkFrom ds i hi
  have h2 := colour_walkFrom ds j hj
  have := colour_ne_of_adjacent h
  omega

/-- Residues with an even index separation are never grid adjacent; this
justifies the parity filter of Python's `_possible_hh_pairs`. -/
theorem not_adjacent_of_even_gap (ds : List (Fin 4)) {i j : ℕ}
    (hi : i < (walkFrom ds).length) (hj : j < (walkFrom ds).length)
    (heven : (j - i) % 2 = 0) (hij : i ≤ j) :
    adjacent ((walkFrom ds)[i]) ((walkFrom ds)[j]) = false := by
  by_contra hc
  have := index_parity_of_adjacent ds hi hj (by simpa using hc)
  omega

/-- A two-step separation can never form a contact. -/
theorem not_adjacent_of_gap_two (ds : List (Fin 4)) {i : ℕ}
    (hi : i + 2 < (walkFrom ds).length) :
    adjacent ((walkFrom ds)[i]) ((walkFrom ds)[i + 2]) = false :=
  not_adjacent_of_even_gap ds (by omega) hi (by omega) (by omega)

/-- Python's candidate contact pairs (`_possible_hh_pairs`): both residues
hydrophobic, index gap at least three, and odd. -/
def candidatePair (S : List Bool) (i j : ℕ) : Bool :=
  decide (i + 3 ≤ j) && decide ((j - i) % 2 = 1) && isHydrophobic S i && isHydrophobic S j

/-- Python's scorer `_contacts(path, _possible_hh_pairs(seq))`: candidate
pairs whose larger index is already placed (`j < len(path)`) and whose
sites are at Manhattan distance one.  It is also meaningful for partial
paths (prefixes). -/
def parityContactCount (S : List Bool) (W : List Site) : ℕ :=
  ∑ i ∈ range S.length, ∑ j ∈ range S.length,
    if candidatePair S i j && decide (j < W.length) &&
       manhattan1 ((W[i]?).getD origin) ((W[j]?).getD origin)
    then 1 else 0

/-- **Contact-count equivalence.** On every decoded walk whose length
matches the sequence, the direct Lean scorer and the parity-filtered
Python-style scorer agree exactly. -/
theorem contactCount_eq_parityContactCount (S : List Bool) (ds : List (Fin 4))
    (hlen : S.length = ds.length + 2) :
    contactCount S (walkFrom ds) = parityContactCount S (walkFrom ds) := by
  unfold contactCount parityContactCount
  have hW : (walkFrom ds).length = S.length := by rw [walkFrom_length]; omega
  refine Finset.sum_congr rfl (fun i hi => Finset.sum_congr rfl (fun j hj => ?_))
  simp only [Finset.mem_range] at hi hj
  have hiW : i < (walkFrom ds).length := by omega
  have hjW : j < (walkFrom ds).length := by omega
  rw [List.getElem?_eq_getElem hiW, List.getElem?_eq_getElem hjW]
  simp only [Option.getD_some, ← adjacent_eq_manhattan1]
  by_cases ha : adjacent ((walkFrom ds)[i]) ((walkFrom ds)[j]) = true
  · have hp := index_parity_of_adjacent ds hiW hjW ha
    simp only [ha, candidatePair, hjW, decide_true, Bool.and_true]
    congr 1
    cases isHydrophobic S i <;> cases isHydrophobic S j <;> simp
    omega
  · simp only [Bool.not_eq_true] at ha
    simp [ha]

/-- Contact-count equivalence for canonical geometric paths. -/
theorem saw_contactCount_eq_parity {steps : ℕ} (S : List Bool)
    (hlen : S.length = steps + 2) (P : CanonicalSAW steps) :
    contactCount S P.1 = parityContactCount S P.1 := by
  have h := walkFrom_pathDirections P.2.2.1 P.2.2.2.1
  rw [← h]
  apply contactCount_eq_parityContactCount
  rw [hlen, pathDirections_length]

/-! ### Rational contact weights -/

/-- Equal contact counts give equal exact rational weights. -/
theorem weightQ_eq_of_contacts_eq {steps : ℕ} (q : ℚ) (seq : Fin (steps + 2) → Bool)
    (c c' : Conformation steps)
    (h : contactCount (List.ofFn seq) (walkFrom (List.ofFn c.1)) =
      contactCount (List.ofFn seq) (walkFrom (List.ofFn c'.1))) :
    weightQ steps q seq c = weightQ steps q seq c' := by
  unfold weightQ
  rw [h]

/-- The weight computed with the Python-style parity scorer coincides with
the Lean weight `q ^ contactCount`. -/
theorem weightQ_eq_parity {steps : ℕ} (q : ℚ) (seq : Fin (steps + 2) → Bool)
    (c : Conformation steps) :
    weightQ steps q seq c =
      q ^ parityContactCount (List.ofFn seq) (walkFrom (List.ofFn c.1)) := by
  rw [weightQ, contactCount_eq_parityContactCount]
  simp

/-- Monotonicity of `q ^ k` in the contact count, for `q ≥ 1`. -/
theorem qpow_mono {q : ℚ} (hq : 1 ≤ q) {k₁ k₂ : ℕ} (h : k₁ ≤ k₂) : q ^ k₁ ≤ q ^ k₂ :=
  pow_le_pow_right₀ hq h

/-- Weights are positive for `q ≥ 1` (indeed for any `q > 0`). -/
theorem qpow_pos {q : ℚ} (hq : 1 ≤ q) (k : ℕ) : 0 < q ^ k :=
  pow_pos (by linarith) k

/-- Weight monotonicity between conformations: more contacts, larger weight. -/
theorem weightQ_mono {steps : ℕ} {q : ℚ} (hq : 1 ≤ q) (seq : Fin (steps + 2) → Bool)
    (c c' : Conformation steps)
    (h : contactCount (List.ofFn seq) (walkFrom (List.ofFn c.1)) ≤
      contactCount (List.ofFn seq) (walkFrom (List.ofFn c'.1))) :
    weightQ steps q seq c ≤ weightQ steps q seq c' :=
  qpow_mono hq h

end LatticeModel
end FoldBounds
