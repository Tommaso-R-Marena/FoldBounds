import Mathlib

/-!
# Exact reporter-code geometry and bipartite optimality inference

These Lean 4 theorems check the reusable algebraic and combinatorial
core of a robust two-fold geometric reporter certificate.

Input PDB identities, finite conformer coordinates, integer distances,
graph vertex sets and actual selected reporter edges are verified by
separate Python code, not imported as Lean proof terms here.
-/

namespace FoldBounds.ReporterCode

/-- If two distance bands have a strict gap and one measurement has
strictly less than half-gap error, the correct low-side classification
cannot be overturned. -/
theorem shorter_band_correct {upper lower observed eps : ℚ}
    (hBand : observed ≤ upper + eps)
    (hNoise : 2 * eps < lower - upper) :
    observed < (upper + lower) / 2 := by
  linarith

/-- The symmetric gap theorem for the longer-distance fold. -/
theorem longer_band_correct {upper lower observed eps : ℚ}
    (hBand : lower - eps ≤ observed)
    (hNoise : 2 * eps < lower - upper) :
    (upper + lower) / 2 < observed := by
  linarith

/-- A pairwise distance error at most twice a site perturbation bound
fits the half-gap condition when 4*r is below the gap. The separately
required Euclidean-distance reverse triangle inequality is classical. -/
theorem site_budget_implies_pair_budget {gap radius : ℚ}
    (h : 4 * radius < gap) :
    2 * radius < gap / 2 := by
  linarith

/-- At most f arbitrary erroneous decisions in a 2f+1-bit code means
the correct votes strictly outnumber the wrong votes. -/
theorem majority_decodes {f correct wrong : ℕ}
    (hTotal : correct + wrong = 2 * f + 1)
    (hWrong : wrong ≤ f) :
    wrong < correct := by
  omega

variable {L R : Type*} [DecidableEq L] [DecidableEq R]

/-- A proposed cover endpoint is chosen deterministically from the left
when possible, or from the right if the left endpoint is not in cover. -/
def pickedCovered (cover : Finset (Sum L R)) (edge : L × R) : Sum L R :=
  if Sum.inl edge.1 ∈ cover then Sum.inl edge.1 else Sum.inr edge.2

/-- Each edge satisfying the cover premise has an endpoint inside cover. -/
theorem picked_in_cover (cover : Finset (Sum L R))
    (edge : L × R)
    (h : Sum.inl edge.1 ∈ cover ∨ Sum.inr edge.2 ∈ cover) :
    pickedCovered cover edge ∈ cover := by
  by_cases hl : Sum.inl edge.1 ∈ cover
  · simp [pickedCovered, hl]
  · have hr : Sum.inr edge.2 ∈ cover := h.resolve_left hl
    simpa [pickedCovered, hl] using hr

/-- A collection of edges with NO repeated left or right endpoints,
each intersecting a cover, has cardinality at most the cover's size.

The proof constructs an explicit injection from distinct matching
edges to cover vertices; it does not assume an optimization solver
was correct. -/
theorem matching_card_le_cover
    (selected : Finset (L × R))
    (cover : Finset (Sum L R))
    (hLeft : ∀ e ∈ selected, ∀ d ∈ selected, e.1 = d.1 → e = d)
    (hRight : ∀ e ∈ selected, ∀ d ∈ selected, e.2 = d.2 → e = d)
    (hCover : ∀ e ∈ selected,
      Sum.inl e.1 ∈ cover ∨ Sum.inr e.2 ∈ cover) :
    selected.card ≤ cover.card := by
  have memPick : ∀ e ∈ selected, pickedCovered cover e ∈ cover := by
    intro e he
    exact picked_in_cover cover e (hCover e he)
  have injPick : ∀ e ∈ selected, ∀ d ∈ selected,
      pickedCovered cover e = pickedCovered cover d → e = d := by
    intro e he d hd h
    by_cases ha : Sum.inl e.1 ∈ cover
    · by_cases hb : Sum.inl d.1 ∈ cover
      · simp [pickedCovered, ha, hb] at h
        exact hLeft e he d hd h
      · simp [pickedCovered, ha, hb] at h
        cases h
    · by_cases hb : Sum.inl d.1 ∈ cover
      · simp [pickedCovered, ha, hb] at h
        cases h
      · simp [pickedCovered, ha, hb] at h
        exact hRight e he d hd h
  let assignment : {e // e ∈ selected} → {v // v ∈ cover} :=
    fun e => ⟨pickedCovered cover e.1, memPick e.1 e.2⟩
  have hInj : Function.Injective assignment := by
    intro e d heq
    apply Subtype.ext
    exact injPick e.1 e.2 d.1 d.2 (congrArg Subtype.val heq)
  have hCard :
      Fintype.card {e // e ∈ selected} ≤
        Fintype.card {v // v ∈ cover} :=
    Fintype.card_le_of_injective assignment hInj
  simpa using hCard

/-- A verified cover with fewer than K endpoints certifies that NO
K-edge residue-disjoint matching exists in the declared stronger
margin graph. This is a mathematical optimality obstruction. -/
theorem no_k_matching_if_small_cover
    (selected : Finset (L × R))
    (cover : Finset (Sum L R)) (K : ℕ)
    (hLeft : ∀ e ∈ selected, ∀ d ∈ selected, e.1 = d.1 → e = d)
    (hRight : ∀ e ∈ selected, ∀ d ∈ selected, e.2 = d.2 → e = d)
    (hCover : ∀ e ∈ selected,
      Sum.inl e.1 ∈ cover ∨ Sum.inr e.2 ∈ cover)
    (hSize : selected.card = K) (hSmall : cover.card < K) :
    False := by
  have hBound := matching_card_le_cover selected cover hLeft hRight hCover
  omega

#print axioms FoldBounds.ReporterCode.shorter_band_correct
#print axioms FoldBounds.ReporterCode.longer_band_correct
#print axioms FoldBounds.ReporterCode.site_budget_implies_pair_budget
#print axioms FoldBounds.ReporterCode.majority_decodes
#print axioms FoldBounds.ReporterCode.matching_card_le_cover
#print axioms FoldBounds.ReporterCode.no_k_matching_if_small_cover

end FoldBounds.ReporterCode
