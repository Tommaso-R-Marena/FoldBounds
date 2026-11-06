module

public import FoldBounds.MotionVerifiedInverse

/-!
# FoldBounds.MotionGNMIntegerInverse

An untrusted, scaled integer inverse witness `(N,D)` proves an *exact*
rational inverse of an integer graph Laplacian `A` when all entries of
`A*N=D*I` are accepted by the executable finite-integer checker. The
witness may come from fraction-free elimination, modular arithmetic,
rank-one contact updates or arbitrary untrusted code.

This closes the numerical residual interval for a FIXED specified graph.
It does not establish that graph topology or protein bytes are correct;
those are separate checked input gates. In particular, this file does not
parse the original PDB/NPZ data inside the Lean kernel.
-/

@[expose] public section
namespace FoldBounds.MotionGNMIntegerInverse

open FoldBounds.MotionHessian
open FoldBounds.MotionVerifiedInverse

/-- Finite executable integer witness: positive scale, EVERY entry of
    A*N = D*I. Proposed N is entirely untrusted. -/
def checkScaledInverse {n : ℕ}
    (A N : Matrix (Fin n) (Fin n) ℤ) (D : ℤ) : Bool :=
  decide (0 < D ∧ ∀ i j : Fin n,
    (∑ k : Fin n, A i k * N k j) =
      if i = j then D else 0)

/-- Every accepted integer witness has its claimed matrix identity. -/
theorem checkScaledInverse_valid {n : ℕ}
    (A N : Matrix (Fin n) (Fin n) ℤ) (D : ℤ)
    (h : checkScaledInverse A N D = true) :
    0 < D ∧ ∀ i j : Fin n,
      (∑ k : Fin n, A i k * N k j) =
        if i = j then D else 0 := by
  exact of_decide_eq_true h

/-- Cast an integer network matrix into rationals. -/
def castQ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℤ) :
    Matrix (Fin n) (Fin n) ℚ :=
  fun i j => (A i j : ℚ)

/-- The exact rational inverse defined by a scaled integer witness. -/
def inverseQ {n : ℕ} (N : Matrix (Fin n) (Fin n) ℤ)
    (D : ℤ) : Matrix (Fin n) (Fin n) ℚ :=
  fun i j => (N i j : ℚ) / (D : ℚ)

/-- No numerical residual: an integer-checked scaled inverse gives all
    right-inverse identities over rationals. -/
theorem checkScaledInverse_to_rational {n : ℕ}
    (A N : Matrix (Fin n) (Fin n) ℤ) (D : ℤ)
    (h : checkScaledInverse A N D = true) :
    ∀ i j : Fin n,
      (∑ k : Fin n, castQ A i k * inverseQ N D k j) =
        if i = j then (1 : ℚ) else 0 := by
  obtain ⟨hD, hentries⟩ := checkScaledInverse_valid A N D h
  have hDq : (D : ℚ) ≠ 0 := by
    exact_mod_cast (ne_of_gt hD)
  intro i j
  have hsum :
      (∑ k : Fin n, (A i k : ℚ) * (N k j : ℚ)) =
        if i = j then (D : ℚ) else 0 := by
    exact_mod_cast (hentries i j)
  change (∑ k : Fin n, (A i k : ℚ) * ((N k j : ℚ) / (D : ℚ))) = _
  calc
    (∑ k : Fin n, (A i k : ℚ) * ((N k j : ℚ) / (D : ℚ))) =
        (∑ k : Fin n, (A i k : ℚ) * (N k j : ℚ)) / (D : ℚ) := by
          rw [Finset.sum_div]
          apply Finset.sum_congr rfl
          intro k _
          ring
    _ = (if i = j then (D : ℚ) else 0) / (D : ℚ) := by
          rw [hsum]
    _ = if i = j then (1 : ℚ) else 0 := by
          split_ifs <;> simp [hDq]

/-- Kernel-checked rational matrix gate follows from the integer gate. -/
theorem checkScaledInverse_rightInverse {n : ℕ}
    (A N : Matrix (Fin n) (Fin n) ℤ) (D : ℤ)
    (h : checkScaledInverse A N D = true) :
    checkRightInverse (castQ A) (inverseQ N D) = true := by
  simpa only [checkRightInverse, decide_eq_true_eq] using
    (checkScaledInverse_to_rational A N D h)

/-- A successful finite integer check constructs the exact response for
    *every* rational force, with no approximation or existence premise. -/
theorem checkScaledInverse_solves {n : ℕ}
    (A N : Matrix (Fin n) (Fin n) ℤ) (D : ℤ)
    (h : checkScaledInverse A N D = true)
    (z : Fin n → ℚ) :
    matVecQ (castQ A) (matVecQ (inverseQ N D) z) = z := by
  exact checked_inverse_solves (castQ A) (inverseQ N D)
    (checkScaledInverse_rightInverse A N D h) z

end FoldBounds.MotionGNMIntegerInverse