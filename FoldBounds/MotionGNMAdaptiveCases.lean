module

public import FoldBounds.MotionGNMTopologyRobust
public import FoldBounds.MotionRobustConcordance

/-!
# FoldBounds.MotionGNMAdaptiveCases

Two distinct mathematical statements needed by adaptive protein-GNM contact
verification. No dataset or file-parser claims are made by these theorems.

(1) For any Boolean realization of a contact network, fixing an arbitrary
selected subset of contacts yields a lower and an upper contact-weight graph
that sandwich that *same realization*. Rayleigh monotonicity then encloses
the realized scalar graph response. In particular, 2^k cases can cover an
arbitrarily large uncertain-contact graph: all unsplit contacts are universally
quantified, not sampled or presumed statistically independent.

(2) If the endpoint intervals in *every* selected-contact case certify a
lower concordance bound, that bound applies to the signed concordance in
*every* full graph realization. This is a finite-combinatorial theorem with
no assumption that per-pair score signs can vary independently.

The protein-coordinate box classifier, graph witness parser, exact integer
interval arithmetic and 2AAS observation hashes are separately rechecked by
independent Python. The Lean kernel does not decode those bytes here.
-/

@[expose] public section
namespace FoldBounds.MotionGNMAdaptiveCases

open FoldBounds.MotionGNMTopologyRobust
open FoldBounds.MotionRobustConcordance
open FoldBounds.MotionHessian

/-- Convert a contact's Boolean presence into a nonnegative rational spring. -/
def binaryWeight (b : Bool) : ℚ := if b then 1 else 0

/-- Which selected contacts are present in this actual full graph? -/
def selectedState {m : ℕ} (chosen actual : Fin m → Bool) : Fin m → Bool :=
  fun e => chosen e && actual e

/-- Contact network in a particular full structural world. -/
def actualWeight {m : ℕ}
    (certain optional actual : Fin m → Bool) : Fin m → ℚ :=
  fun e => binaryWeight (certain e || (optional e && actual e))

/-- Every unselected optional edge is OFF; selected optional edges take
    their selected-case value, while certain contacts are always ON. -/
def caseLowWeight {m : ℕ}
    (certain optional chosen caseMask : Fin m → Bool) : Fin m → ℚ :=
  fun e => binaryWeight
    (certain e || (optional e && chosen e && caseMask e))

/-- Every unselected optional edge is ON; selected optional edges still
    take their selected-case value, while certain contacts are always ON. -/
def caseHighWeight {m : ℕ}
    (certain optional chosen caseMask : Fin m → Bool) : Fin m → ℚ :=
  fun e => binaryWeight
    (certain e || (optional e && (!chosen e || caseMask e)))

/-- A branch's extremal contact networks include its actual world.
    This holds for ANY choice of split edges and ANY full edge assignment. -/
theorem case_contact_sandwich {m : ℕ}
    (certain optional chosen actual : Fin m → Bool) (e : Fin m) :
    caseLowWeight certain optional chosen (selectedState chosen actual) e ≤
      actualWeight certain optional actual e ∧
    actualWeight certain optional actual e ≤
      caseHighWeight certain optional chosen (selectedState chosen actual) e := by
  cases hc : certain e <;> cases ho : optional e <;>
    cases hs : chosen e <;> cases ha : actual e <;>
    simp [caseLowWeight, caseHighWeight, actualWeight,
      selectedState, binaryWeight, hc, ho, hs, ha]

/-- The Rayleigh lower/upper sandwich holds for arbitrary full contact states
    after the selected edges are fixed. Endpoint optimal-response premises
    are discharged by the separately checked exact graph certificates. -/
theorem case_response_sandwich {m : ℕ} {I : Type*}
    [Fintype I] [DecidableEq I]
    (b : Fin m → I → ℚ)
    (certain optional chosen actual : Fin m → Bool)
    (z : I → ℚ) (qOff qActual qOn : ℚ)
    (hOff : ∀ u, scalarResponseDual b
      (caseLowWeight certain optional chosen (selectedState chosen actual))
      z u ≤ qOff)
    (hActual : ∀ u, scalarResponseDual b
      (actualWeight certain optional actual) z u ≤ qActual)
    (hActualAttains : ∃ u, scalarResponseDual b
      (actualWeight certain optional actual) z u = qActual)
    (hOnAttains : ∃ u, scalarResponseDual b
      (caseHighWeight certain optional chosen (selectedState chosen actual))
      z u = qOn) :
    qOn ≤ qActual ∧ qActual ≤ qOff := by
  exact scalarResponse_extremes b
    (caseLowWeight certain optional chosen (selectedState chosen actual))
    (actualWeight certain optional actual)
    (caseHighWeight certain optional chosen (selectedState chosen actual))
    z qOff qActual qOn
    (fun e => case_contact_sandwich certain optional chosen actual e)
    hOff hActual hActualAttains hOnAttains

/-- Universal signed-concordance certificate for a family of full contact
    networks. For each state, the selected-edge mask determines an interval
    for every oriented score comparison. If EVERY case has a certified
    minimum signed count, then EVERY full network meets that same minimum.

    Counts are unnormalized; dividing by a fixed positive number of target-
    non-tied comparisons preserves the bound. -/
theorem all_contact_worlds_concordance {m k comparisons : ℕ}
    (pick : Fin k → Fin m)
    (scoreA scoreB : (Fin m → Bool) → Fin comparisons → ℚ)
    (aLo aHi bLo bHi : (Fin k → Bool) → Fin comparisons → ℚ)
    (guaranteed : ℚ)
    (hA : ∀ world i,
      aLo (fun j => world (pick j)) i ≤ scoreA world i ∧
      scoreA world i ≤ aHi (fun j => world (pick j)) i)
    (hB : ∀ world i,
      bLo (fun j => world (pick j)) i ≤ scoreB world i ∧
      scoreB world i ≤ bHi (fun j => world (pick j)) i)
    (hAllCases : ∀ mask : Fin k → Bool,
      guaranteed ≤
        ∑ i : Fin comparisons,
          if aHi mask i < bLo mask i then (1 : ℚ) else -1) :
    ∀ world : Fin m → Bool,
      guaranteed ≤ ∑ i : Fin comparisons,
        alignedSign (scoreA world i) (scoreB world i) := by
  intro world
  let mask : Fin k → Bool := fun j => world (pick j)
  have hscore := sum_alignedSign_interval
    (aLo mask) (aHi mask) (bLo mask) (bHi mask)
    (scoreA world) (scoreB world)
    (fun i => hA world i) (fun i => hB world i)
  exact le_trans (hAllCases mask) hscore.1

end FoldBounds.MotionGNMAdaptiveCases