#!/usr/bin/env python3
"""Check that the strict Lean-side audit rejects planted proof bypasses.

Runs `lake env lean tests/kernel_fixtures/axiom_gate/GateAttack.lean`, which
must FAIL, and verifies that each planted theorem is rejected for the right
reason while the honest control is accepted.  Exit status 0 means the audit
command behaves as required.
"""
from __future__ import annotations

import subprocess
import sys

FIXTURE = "tests/kernel_fixtures/axiom_gate/GateAttack.lean"
EXPECTED = (
    "attack_hole depends on non-standard axioms: [sorryAx]",
    "attack_native depends on non-standard axioms: [Lean.ofReduceBool",
    "attack_lookalike depends on non-standard axioms: [«Classical.choice»]",
)


def main() -> int:
    result = subprocess.run(["lake", "env", "lean", FIXTURE], capture_output=True, text=True,
                            timeout=900)
    output = result.stdout + result.stderr
    problems = []
    if result.returncode == 0:
        problems.append("the attack fixture elaborated without error")
    for needle in EXPECTED:
        if needle not in output:
            problems.append(f"missing expected rejection: {needle}")
    if "honest_control depends on" in output:
        problems.append("the honest control was rejected")
    if "standard axioms only: 1 listed declarations" not in output:
        problems.append("the honest control was not reported as accepted")
    if problems:
        print(output)
        for problem in problems:
            print("GATE ATTACK CHECK FAILED:", problem)
        return 1
    print("Strict audit rejects planted sorry, native_decide and look-alike dependencies; "
          "honest control accepted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
