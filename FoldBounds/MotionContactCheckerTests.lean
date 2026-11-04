import FoldBounds.MotionContactChecker

open FoldBounds.MotionContactChecker

private def exampleCert : RawContacts 3 :=
  { cutoffSq := 100
    d2Lower := ![49, 98, 121]
    d2Upper := ![64, 102, 144]
    kLower := ![1, 2, 3]
    kUpper := ![1, 3, 3] }

private def exampleD2 : Fin 3 → ℚ := ![60, 99, 125]
private def exampleK : Fin 3 → ℚ := ![1, 2, 3]

-- Exact rational acceptance and three-way contact classification:
-- 0 always on, 1 uncertain, 2 always off.
example : accepts exampleCert = true := by decide
example : check exampleCert exampleD2 exampleK = true := by decide
example : lowerWeight exampleCert 0 = 1 := by decide
example : upperWeight exampleCert 0 = 1 := by decide
example : lowerWeight exampleCert 1 = 0 := by decide
example : upperWeight exampleCert 1 = 3 := by decide
example : lowerWeight exampleCert 2 = 0 := by decide
example : upperWeight exampleCert 2 = 0 := by decide

-- Fail-closed adversarial invalid inputs.
example : check exampleCert ![60, 110, 125] exampleK = false := by decide
example : check exampleCert exampleD2 ![1, -1, 3] = false := by decide
example : accepts { exampleCert with d2Lower := ![70, 98, 121] } = false := by decide
example : accepts { exampleCert with kLower := ![-1, 2, 3] } = false := by decide
example : accepts { exampleCert with cutoffSq := 0 } = false := by decide

-- Checker's soundness is about all inputs, not only these concrete examples.
example (c : RawContacts 3) (d2 k : Fin 3 → ℚ)
    (h : check c d2 k = true) :
    ∀ i, lowerWeight c i ≤ actualWeight c d2 k i ∧
      actualWeight c d2 k i ≤ upperWeight c i := by
  have hh := check_valid c d2 k h
  exact weight_interval_sound c d2 k hh.1 hh.2
