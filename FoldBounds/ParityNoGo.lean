import Mathlib

/-!
# Abstract bipartite-graph contact obstruction

A graph edge joins two different colors. If every H residue is
assigned to the same one of the two colors, no H-H edge can exist.

This is a Lean-kernel-checked abstract logical implication. In the
Python checker the graph is the simple cubic lattice and the coloring
is (x+y+z) mod 2, which flips along every unit lattice step.
The geometry, sequence-index-to-lattice-color relation and global
3^(n-2) SAW lower construction are NOT yet formalized in Lean.
-/

namespace FoldBounds.ParityNoGo

theorem same_color_hydrophobes_cannot_contact
    {V : Type*} (isH : V → Prop) (color : V → Bool)
    (adjacent : V → V → Prop)
    (hHydrophobic : ∀ x, isH x → color x = true)
    (hBipartite : ∀ x y, adjacent x y → color x ≠ color y)
    (x y : V) (hx : isH x) (hy : isH y) :
    ¬ adjacent x y := by
  intro edge
  have same : color x = color y :=
    (hHydrophobic x hx).trans (hHydrophobic y hy).symm
  exact (hBipartite x y edge) same

/-- Even when the chosen H-color is not specified in advance, the
absence of same-color contact follows from the bipartite premise. -/
theorem fixed_color_class_has_no_internal_edges
    {V : Type*} (isH : V → Prop) (color : V → Bool)
    (adjacent : V → V → Prop) (chosen : Bool)
    (hHydrophobic : ∀ x, isH x → color x = chosen)
    (hBipartite : ∀ x y, adjacent x y → color x ≠ color y)
    (x y : V) (hx : isH x) (hy : isH y) :
    ¬ adjacent x y := by
  intro edge
  have same : color x = color y :=
    (hHydrophobic x hx).trans (hHydrophobic y hy).symm
  exact (hBipartite x y edge) same

#print axioms FoldBounds.ParityNoGo.same_color_hydrophobes_cannot_contact
#print axioms FoldBounds.ParityNoGo.fixed_color_class_has_no_internal_edges

end FoldBounds.ParityNoGo
