module

public import FoldBounds.Certificate

/-!
# FoldBounds.Executable — an executable rational certificate-checking kernel

This checks finite-model witness membership and individual upper-weight
obligations, then returns exact rational event/complement enclosures.
The proof is about the actual decidable checker, not an assumed oracle.
It still enumerates the explicit finite universe to validate caps and
does NOT prove correctness of the separate Python lattice implementation.
-/

@[expose] public section

open Finset
namespace FoldBounds
namespace Executable

variable {C : Type*} [DecidableEq C]

/-- Exact rational mass for an explicitly finite collection. -/
def massQ (w : C → ℚ) (T : Finset C) : ℚ := ∑ c ∈ T, w c

def lowerQ (w : C → ℚ) (S : Finset C) : ℚ := massQ w S

def upperQ (w cap : C → ℚ) (S T : Finset C) : ℚ :=
  massQ w S + massQ cap (T \ S)

/-- Sound subset-sum lower bounds and cap-based upper bounds. -/
theorem massQ_mem_Icc (w cap : C → ℚ) {S T : Finset C} (hST : S ⊆ T)
    (hpos : ∀ c ∈ T, 0 ≤ w c)
    (hcap : ∀ c ∈ T \ S, w c ≤ cap c) :
    massQ w T ∈ Set.Icc (lowerQ w S) (upperQ w cap S T) := by
  have hlo : (∑ c ∈ S, w c) ≤ ∑ c ∈ T, w c :=
    Finset.sum_le_sum_of_subset_of_nonneg hST (fun c hc _ => hpos c hc)
  have htail : (∑ c ∈ T \ S, w c) ≤ ∑ c ∈ T \ S, cap c :=
    Finset.sum_le_sum (fun c hc => hcap c hc)
  have hsplit : (∑ c ∈ T \ S, w c) + ∑ c ∈ S, w c = ∑ c ∈ T, w c :=
    Finset.sum_sdiff hST
  constructor
  · simpa only [lowerQ, massQ] using hlo
  · change (∑ c ∈ T, w c) ≤
      (∑ c ∈ S, w c) + ∑ c ∈ T \ S, cap c
    linarith

/-- All certificate fields are untrusted inputs from a proposer. -/
structure RawCertificate (C : Type*) where
  eventWitness : Finset C
  otherWitness : Finset C
  cap : C → ℚ

variable [Fintype C]

/-- Executable rejection gate: subsets, positivity, and caps are decided
against the explicit universe, not trusted on the prover's authority. -/
def accepts (w : C → ℚ) (A : Finset C) (cert : RawCertificate C) : Bool :=
  decide (cert.eventWitness ⊆ A ∧
    cert.otherWitness ⊆ Aᶜ ∧
    (∀ c ∈ A, 0 ≤ w c) ∧
    (∀ c ∈ Aᶜ, 0 ≤ w c) ∧
    (∀ c ∈ A \ cert.eventWitness, w c ≤ cert.cap c) ∧
    (∀ c ∈ Aᶜ \ cert.otherWitness, w c ≤ cert.cap c) ∧
    (∀ c : C, 0 < w c))

structure RationalBounds where
  eventLower : ℚ
  eventUpper : ℚ
  otherLower : ℚ
  otherUpper : ℚ
  deriving DecidableEq, Repr

def endpoints (w : C → ℚ) (A : Finset C) (cert : RawCertificate C) : RationalBounds :=
  ⟨lowerQ w cert.eventWitness,
   upperQ w cert.cap cert.eventWitness A,
   lowerQ w cert.otherWitness,
   upperQ w cert.cap cert.otherWitness Aᶜ⟩

/-- Fully executable rational finite-model certificate checker. -/
def check (w : C → ℚ) (A : Finset C) (cert : RawCertificate C) :
    Option RationalBounds :=
  if accepts w A cert then some (endpoints w A cert) else none

/-- Every accepted certificate soundly encloses both disjoint masses. -/
theorem accepts_sound (w : C → ℚ) (A : Finset C) (cert : RawCertificate C)
    (h : accepts w A cert = true) :
    massQ w A ∈ Set.Icc (endpoints w A cert).eventLower (endpoints w A cert).eventUpper ∧
    massQ w Aᶜ ∈ Set.Icc (endpoints w A cert).otherLower (endpoints w A cert).otherUpper := by
  have hc : cert.eventWitness ⊆ A ∧
      cert.otherWitness ⊆ Aᶜ ∧
      (∀ c ∈ A, 0 ≤ w c) ∧
      (∀ c ∈ Aᶜ, 0 ≤ w c) ∧
      (∀ c ∈ A \ cert.eventWitness, w c ≤ cert.cap c) ∧
      (∀ c ∈ Aᶜ \ cert.otherWitness, w c ≤ cert.cap c) ∧
      (∀ c : C, 0 < w c) := of_decide_eq_true h
  rcases hc with ⟨hS, hR, hposA, hposB, hcapA, hcapB, _⟩
  exact ⟨massQ_mem_Icc w cert.cap hS hposA hcapA,
         massQ_mem_Icc w cert.cap hR hposB hcapB⟩

/-- End-to-end theorem about actual Option-returning executable acceptance. -/
theorem check_sound (w : C → ℚ) (A : Finset C) (cert : RawCertificate C)
    (bounds : RationalBounds) (h : check w A cert = some bounds) :
    bounds = endpoints w A cert ∧
    massQ w A ∈ Set.Icc bounds.eventLower bounds.eventUpper ∧
    massQ w Aᶜ ∈ Set.Icc bounds.otherLower bounds.otherUpper := by
  unfold check at h
  by_cases haccept : accepts w A cert = true
  · simp [haccept] at h
    subst bounds
    exact ⟨rfl, (accepts_sound w A cert haccept).1,
      (accepts_sound w A cert haccept).2⟩
  · have hfalse : accepts w A cert = false := by
      cases hval : accepts w A cert <;> simp_all
    simp [hfalse] at h


/-- An accepted rational certificate has positive weight at every modeled state. -/
theorem accepts_weights_pos (w : C → ℚ) (A : Finset C)
    (cert : RawCertificate C) (h : accepts w A cert = true) (c : C) : 0 < w c := by
  have hc : cert.eventWitness ⊆ A ∧
      cert.otherWitness ⊆ Aᶜ ∧
      (∀ c ∈ A, 0 ≤ w c) ∧
      (∀ c ∈ Aᶜ, 0 ≤ w c) ∧
      (∀ c ∈ A \ cert.eventWitness, w c ≤ cert.cap c) ∧
      (∀ c ∈ Aᶜ \ cert.otherWitness, w c ≤ cert.cap c) ∧
      (∀ c : C, 0 < w c) := of_decide_eq_true h
  rcases hc with ⟨_, _, _, _, _, _, hw⟩
  exact hw c

/-- The accepted rational mass certificate also implies an exact-real
coupled probability enclosure. Unlike separately bounding numerator and Z,
this reuses the numerator in its denominator (Aristotle's ratio theorem).
The rational-to-real casts introduce no floating-point computation. -/
theorem accepts_coupled_probability_sound [Nonempty C] (w : C → ℚ)
    (A : Finset C) (cert : RawCertificate C)
    (h : accepts w A cert = true) :
    (massQ w A : ℝ) / ((massQ w A : ℝ) + (massQ w Aᶜ : ℝ)) ∈
      Set.Icc
        ((endpoints w A cert).eventLower /
          ((endpoints w A cert).eventLower + (endpoints w A cert).otherUpper) : ℝ)
        ((endpoints w A cert).eventUpper /
          ((endpoints w A cert).eventUpper + (endpoints w A cert).otherLower) : ℝ) := by
  have masses := accepts_sound w A cert h
  have hA : (massQ w A : ℝ) ∈
      Set.Icc ((endpoints w A cert).eventLower : ℝ)
        ((endpoints w A cert).eventUpper : ℝ) := by
    exact ⟨by exact_mod_cast masses.1.1, by exact_mod_cast masses.1.2⟩
  have hB : (massQ w Aᶜ : ℝ) ∈
      Set.Icc ((endpoints w A cert).otherLower : ℝ)
        ((endpoints w A cert).otherUpper : ℝ) := by
    exact ⟨by exact_mod_cast masses.2.1, by exact_mod_cast masses.2.2⟩
  have hLowerA : (0 : ℚ) ≤ (endpoints w A cert).eventLower := by
    change (0 : ℚ) ≤ ∑ c ∈ cert.eventWitness, w c
    exact Finset.sum_nonneg (fun c _ => (accepts_weights_pos w A cert h c).le)
  have hLowerB : (0 : ℚ) ≤ (endpoints w A cert).otherLower := by
    change (0 : ℚ) ≤ ∑ c ∈ cert.otherWitness, w c
    exact Finset.sum_nonneg (fun c _ => (accepts_weights_pos w A cert h c).le)
  have hZ : (0 : ℚ) < massQ w (Finset.univ : Finset C) := by
    unfold massQ
    exact Finset.sum_pos (fun c _ => accepts_weights_pos w A cert h c)
      Finset.univ_nonempty
  have hSplit : massQ w A + massQ w Aᶜ =
      massQ w (Finset.univ : Finset C) := by
    unfold massQ
    exact Finset.sum_add_sum_compl A w
  have hDenomQ : (0 : ℚ) < massQ w A + massQ w Aᶜ := by
    rw [hSplit]
    exact hZ
  have hDenom : (0 : ℝ) < (massQ w A : ℝ) + (massQ w Aᶜ : ℝ) := by
    exact_mod_cast hDenomQ
  exact FoldBounds.ratio_mem_coupled (by exact_mod_cast hLowerA)
    (by exact_mod_cast hLowerB) hA hB hDenom


/-- The actual executable returned a certificate only if the full rational
decidable gate accepted it. No extra trusted premise about the untrusted input is required. -/
theorem check_implies_accepts (w : C → ℚ) (A : Finset C)
    (cert : RawCertificate C) (bounds : RationalBounds)
    (h : check w A cert = some bounds) :
    accepts w A cert = true := by
  unfold check at h
  cases hval : accepts w A cert with
  | false => simp [hval] at h
  | true => exact hval

/-- Strongest end-to-end finite-model result: an actual successful call to
the runnable checker guarantees the coupled event-probability enclosure.
It does not assume the validity of the untrusted certificate fields. -/
theorem check_coupled_probability_sound [Nonempty C] (w : C → ℚ)
    (A : Finset C) (cert : RawCertificate C) (bounds : RationalBounds)
    (h : check w A cert = some bounds) :
    (massQ w A : ℝ) / ((massQ w A : ℝ) + (massQ w Aᶜ : ℝ)) ∈
      Set.Icc
        ((bounds.eventLower / (bounds.eventLower + bounds.otherUpper)) : ℝ)
        ((bounds.eventUpper / (bounds.eventUpper + bounds.otherLower)) : ℝ) := by
  have haccept := check_implies_accepts w A cert bounds h
  have hEq := (check_sound w A cert bounds h).1
  subst bounds
  exact accepts_coupled_probability_sound w A cert haccept


end Executable
end FoldBounds
