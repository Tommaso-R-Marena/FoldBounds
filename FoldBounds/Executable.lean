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

end Executable
end FoldBounds
