import Mathlib

/-!
# Local fold-rewiring geometric core: multiset occupancy and injectivity

This module isolates the finite-list obligations behind a checked
eight-site cube rewire. If a proposed replacement is a permutation
of a site's eight visited coordinates, then replacing that sublist
does not introduce a repeated lattice site in an originally self-
avoiding global walk. A fixed pre/suffix splice also preserves
distinguishability of different local variants, so each may count
once in a constructive lower partition mass.

The premise that a particular replacement has the same endpoints
and traverses only adjacent unit cubic edges is NOT derived here.
Those geometric conditions must be checked separately, as in the
Python eight-site motif verifier. Nor does this module establish
independence of two interacting cubic neighborhoods.
-/

namespace FoldBounds.CubicSplice

variable {α : Type*}

/-- Local permutation preserves the entire visited-site multiset
when the pre and suffix remain fixed. -/
theorem splice_perm
    (pre original replacement suffix : List α)
    (h : original.Perm replacement) :
    (pre ++ original ++ suffix).Perm
      (pre ++ replacement ++ suffix) := by
  exact (h.append_left pre).append_right suffix

/-- Formalize the nonintersection argument for one local Hamiltonian
rearrangement: the complete list remains duplicate-free. -/
theorem splice_preserves_nodup
    (pre original replacement suffix : List α)
    (hn : (pre ++ original ++ suffix).Nodup)
    (hp : original.Perm replacement) :
    (pre ++ replacement ++ suffix).Nodup := by
  exact (splice_perm pre original replacement suffix hp).nodup_iff.mp hn

/-- The complete walk visits precisely the same sites before and after
replacing the local fragment, independent of the residue labels. -/
theorem splice_preserves_occupied_sites
    (pre original replacement suffix : List α)
    (hp : original.Perm replacement)
    (site : α) :
    site ∈ pre ++ original ++ suffix ↔
      site ∈ pre ++ replacement ++ suffix := by
  exact (splice_perm pre original replacement suffix hp).mem_iff

/-- Two different candidate local paths cannot become the same
whole conformational coordinate list when inserted between identical
global prefixes and suffixes. -/
theorem splice_injective
    (pre suffix : List α) :
    Function.Injective (fun patch : List α => pre ++ patch ++ suffix) := by
  intro first second h
  have hWithoutSuffix :
      pre ++ first = pre ++ second := List.append_cancel_right h
  exact List.append_cancel_left hWithoutSuffix

/-- Distinct checked local variants give exactly as many distinct
complete rooted conformations: no hidden multiplicity factor. -/
theorem splice_family_card
    [DecidableEq α]
    (pre suffix : List α)
    (variants : Finset (List α)) :
    (variants.image
      (fun patch => pre ++ patch ++ suffix)).card = variants.card := by
  exact Finset.card_image_of_injective variants (splice_injective pre suffix)

/-- A permutation replacement cannot alter overall chain length. -/
theorem splice_preserves_length
    (pre original replacement suffix : List α)
    (hp : original.Perm replacement) :
    (pre ++ original ++ suffix).length =
      (pre ++ replacement ++ suffix).length := by
  simp [hp.length_eq]

#print axioms FoldBounds.CubicSplice.splice_perm
#print axioms FoldBounds.CubicSplice.splice_preserves_nodup
#print axioms FoldBounds.CubicSplice.splice_preserves_occupied_sites
#print axioms FoldBounds.CubicSplice.splice_injective
#print axioms FoldBounds.CubicSplice.splice_family_card
#print axioms FoldBounds.CubicSplice.splice_preserves_length

end FoldBounds.CubicSplice
