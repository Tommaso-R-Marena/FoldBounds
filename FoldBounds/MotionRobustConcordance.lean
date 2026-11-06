module

public import FoldBounds.MotionRobustOrdering

/-!
# FoldBounds.MotionRobustConcordance

Conditional, kernel-checked certification of a *global directional concordance*
from coordinate-world-uniform pair score intervals.

For each ordered comparison, the first score should be below the second,
as specified by a frozen target-label ordering. The contribution is +1 for
correct order, -1 for reverse order, and 0 for an exact tie.

Disjoint score intervals fix that contribution for every possible world.
Otherwise the checker conservatively assigns the pair [-1,+1].
Summing these bounds certifies total concordance without assuming that all
uncertain score rankings can be realized independently. In particular, the
aggregate interval may be conservative rather than globally sharp.

This theorem does NOT verify Python coordinate-box enclosures, target NMR
ranks, parser/byte provenance, or physical protein dynamics.
-/

@[expose] public section
namespace FoldBounds.MotionRobustConcordance

/-- Target-order-aligned sign: +1 if x<y, -1 if x>y, 0 if tied. -/
def alignedSign (x y : ℚ) : ℚ :=
  if x < y then 1 else if y < x then -1 else 0

theorem alignedSign_bounds (x y : ℚ) :
    -1 ≤ alignedSign x y ∧ alignedSign x y ≤ 1 := by
  unfold alignedSign
  split_ifs <;> norm_num

theorem alignedSign_of_lt {x y : ℚ} (h : x < y) :
    alignedSign x y = 1 := by
  simp [alignedSign, h]

theorem alignedSign_of_gt {x y : ℚ} (h : y < x) :
    alignedSign x y = -1 := by
  have hn : ¬ x < y := not_lt.mpr (le_of_lt h)
  simp [alignedSign, hn, h]

/-- A fail-closed pair enclosure valid even if contact networks switch. -/
theorem alignedSign_interval
    (loA hiA loB hiB x y : ℚ)
    (ha : loA ≤ x ∧ x ≤ hiA)
    (hb : loB ≤ y ∧ y ≤ hiB) :
    (if hiA < loB then (1 : ℚ) else -1) ≤ alignedSign x y ∧
    alignedSign x y ≤ (if hiB < loA then (-1 : ℚ) else 1) := by
  constructor
  · split_ifs with h
    · have hxy : x < y := lt_of_le_of_lt ha.2 (lt_of_lt_of_le h hb.1)
      simpa [alignedSign_of_lt hxy]
    · exact (alignedSign_bounds x y).1
  · split_ifs with h
    · have hyx : y < x := lt_of_le_of_lt hb.2 (lt_of_lt_of_le h ha.1)
      simpa [alignedSign_of_gt hyx]
    · exact (alignedSign_bounds x y).2

/-- Global concordance bounds from any finite family of intervals. Each
index represents one pair of score comparisons already oriented according
to the frozen target's desired order. No independence is assumed. -/
theorem sum_alignedSign_interval {m : ℕ}
    (loA hiA loB hiB x y : Fin m → ℚ)
    (ha : ∀ i, loA i ≤ x i ∧ x i ≤ hiA i)
    (hb : ∀ i, loB i ≤ y i ∧ y i ≤ hiB i) :
    ((∑ i : Fin m, if hiA i < loB i then (1 : ℚ) else -1) ≤
        (∑ i : Fin m, alignedSign (x i) (y i))) ∧
    ((∑ i : Fin m, alignedSign (x i) (y i)) ≤
        (∑ i : Fin m, if hiB i < loA i then (-1 : ℚ) else 1)) := by
  constructor
  · apply Finset.sum_le_sum
    intro i _
    exact (alignedSign_interval (loA i) (hiA i) (loB i) (hiB i)
      (x i) (y i) (ha i) (hb i)).1
  · apply Finset.sum_le_sum
    intro i _
    exact (alignedSign_interval (loA i) (hiA i) (loB i) (hiB i)
      (x i) (y i) (ha i) (hb i)).2

end FoldBounds.MotionRobustConcordance
