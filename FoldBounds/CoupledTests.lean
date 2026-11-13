import FoldBounds.CoupledBinding

/-! Concrete accept/reject tests for the CEC checkers.  Every Boolean or `Option` fact below is
checked by kernel evaluation (`decide +kernel`); the elaborator's default reduction cannot
unfold the well-founded `Nat.gcd` used by `ℚ` normalisation, so the check is delegated to the
kernel, which evaluates it directly.  No compiled code is trusted. -/

open Finset
open FoldBounds.CoupledRatio
open FoldBounds.CoupledSeparator

/-! ### §1: coupled linear-fractional checker -/

/-- Two regions with exact ratios 1 and 2 and masses in `[1, 2]`.
Sharp enclosure `[4/3, 5/3]`; separate enclosure `[Σ r L / Σ U, Σ R U / Σ L] = [3/4, 3]`. -/
private def twoRegions : List (ℚ × ℚ × ℚ × ℚ) := [(1, 2, 1, 1), (1, 2, 2, 2)]

-- the sharp endpoints are accepted
example : checkLower twoRegions (4/3) = true := by decide +kernel
example : checkUpper twoRegions (5/3) = true := by decide +kernel
example : checkEnclosure twoRegions (4/3) (5/3) = some (4/3, 5/3) := by decide +kernel
-- the (looser) separate enclosure is also valid
example : checkEnclosure twoRegions (3/4) 3 = some (3/4, 3) := by decide +kernel
-- adversarial: endpoints tightened past the sharp values are rejected
example : checkLower twoRegions (4/3 + 1/1000) = false := by decide +kernel
example : checkUpper twoRegions (5/3 - 1/1000) = false := by decide +kernel
example : checkEnclosure twoRegions (4/3 + 1/1000) (5/3) = none := by decide +kernel
example : checkEnclosure twoRegions (4/3) (5/3 - 1/1000) = none := by decide +kernel
-- adversarial: an inverted interval is rejected
example : checkEnclosure twoRegions 2 1 = none := by decide +kernel

/-- Theorem 1.6: a uniform ratio makes the enclosure a point, whatever the mass bounds. -/
private def uniformRegions : List (ℚ × ℚ × ℚ × ℚ) := [(1, 3, 1/2, 1/2), (0, 5, 1/2, 1/2)]
example : checkEnclosure uniformRegions (1/2) (1/2) = some (1/2, 1/2) := by decide +kernel
example : checkEnclosure uniformRegions (1/2 + 1/10^6) (1/2) = none := by decide +kernel

-- fail closed on malformed input
example : checkLower [(2, 1, 0, 0)] (-5) = false := by decide +kernel      -- L > U
example : checkLower [(-1, 1, 1, 1)] 0 = false := by decide +kernel        -- L < 0
example : checkUpper [(1, 1, 2, 1)] 5 = false := by decide +kernel         -- r > R
example : checkUpper [(1, 1, -1, 1)] 5 = false := by decide +kernel        -- r < 0
example : checkEnclosure [(1, 2, 1, 1), (2, 1, 1, 1)] 0 10 = none := by decide +kernel
-- vacuous data: empty admissible set (every `U = 0`) or no regions
example : checkLower [(0, 0, 1, 1)] 7 = true := by decide +kernel          -- vacuously sound
example : checkEnclosure [(0, 0, 1, 1)] 7 0 = none := by decide +kernel    -- but no enclosure
example : checkEnclosure [] 0 0 = none := by decide +kernel

/-- End-to-end use of the soundness theorem on real configurations. -/
example (α β : Fin twoRegions.length → ℝ) (h : ListAdm twoRegions α β) :
    (4/3 : ℝ) ≤ Q univ α β ∧ Q univ α β ≤ 5/3 := by
  obtain ⟨hb, -, hall⟩ := checkEnclosure_sound (show checkEnclosure twoRegions (4/3) (5/3) =
    some (4/3, 5/3) by decide +kernel)
  have := hall α β h
  norm_num [Set.mem_Icc] at this
  exact this

/-! ### §2: separator checker -/

-- local ratios 2, 1/2, 1 over three boundary states
example : checkSeparator [(1, 2), (2, 1), (1, 1)] = some (1/2, 2) := by decide +kernel
-- an inert state (Φ_A = Φ_B = 0) is ignored
example : checkSeparator [(0, 0), (2, 1)] = some (1/2, 1/2) := by decide +kernel
-- support violation Φ_A = 0 < Φ_B: the sharp upper bound is +∞, so reject
example : checkSeparator [(0, 1), (1, 1)] = none := by decide +kernel
-- negative weights are rejected
example : checkSeparator [(1, -1), (1, 1)] = none := by decide +kernel
example : checkSeparator [(-1, 1)] = none := by decide +kernel
-- empty support T⁺ or no boundary states: no enclosure
example : checkSeparator [(0, 0)] = none := by decide +kernel
example : checkSeparator [] = none := by decide +kernel
-- adversarial: a claimed tighter interval is not what the checker returns
example : checkSeparator [(1, 2), (2, 1), (1, 1)] ≠ some (1/2, 3/2) := by decide +kernel
