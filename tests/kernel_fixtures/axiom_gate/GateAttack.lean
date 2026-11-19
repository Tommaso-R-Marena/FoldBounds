import FoldBounds.FrontierKernelAudit

/-!
Regression fixture for the strict Lean-side audit command.  This file is
**expected to fail**: it plants the three ways a proof can escape a textual
audit (a proof hole, a compiler-trusted `native_decide`, and a forged
dependency whose name is spelled like a standard one) and asks the audit
command to accept them.  `scripts/check_kernel_gate_attack.py` checks that
each planted theorem is rejected with the expected message and that the
honest control is accepted.  It is not part of the `FoldBounds` library.
-/

open FoldBounds.FrontierKernelAudit

theorem attack_hole : (2 : Nat) + 2 = 5 := by sorry

theorem attack_native : (10 : Nat) ^ 3 = 1000 := by native_decide

axiom «Classical.choice» : False

theorem attack_lookalike : (0 : Nat) = 1 := «Classical.choice».elim

theorem honest_control : (1 : Nat) + 1 = 2 := rfl

#assert_standard_axioms honest_control
#assert_standard_axioms attack_hole
#assert_standard_axioms attack_native
#assert_standard_axioms attack_lookalike
