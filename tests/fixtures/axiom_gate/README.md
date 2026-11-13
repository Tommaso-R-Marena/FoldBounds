# Axiom-gate regression fixtures

`real_exploit_reports.log` contains five lines of **real** `#print axioms` output from Lean 4.28.0
(Mathlib v4.28.0), captured on 2026-11-13 from a throwaway build that was never committed. That
build added three files:

* `FoldBounds/Sub/Hole.lean` held a `sorry` in a subdirectory, which the non-recursive source
  scan does not see.
* `FoldBounds/Sub/Forged.lean` held `axiom «Classical.choice» : False`, a single name component
  containing a dot.
* Two clean-looking theorems in scanned files used them.

The first report wraps across lines, with `sorryAx` on a continuation line (finding F3e). The
second report lists the forged `«Classical.choice»` (finding F3f). The previous
`check_lean_sources.py` accepted both, printing "standard axioms only". The tests in
`tests/test_check_lean_sources.py` require the hardened script to reject them.
