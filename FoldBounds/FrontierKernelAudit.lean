import Lean
import FoldBounds.FrontierKernelTests

/-!
# Strict audit of the frontier-kernel development

This file fails to elaborate (and therefore fails `lake build` and CI)
unless **every** theorem and definition of the listed modules depends only
on the three standard axioms `propext`, `Classical.choice` and `Quot.sound`.

Differences from the textual `#print axioms` gate in `scripts/`:

* the check runs inside Lean on `Lean.collectAxioms` and compares `Name`s
  structurally, so it cannot be confused by wrapped reports or by look-alike
  names such as a single-component name spelled like `Classical.choice`
  (the gate findings F3e/F3f fixed for the Python script in draft PR #97);
* it enumerates the modules' declarations itself instead of relying on a
  hand-maintained list, so a new theorem cannot escape the audit;
* `sorryAx` and `Lean.ofReduceBool` (the trace of `native_decide`) are
  rejected like any other non-standard dependency.

The upstream lemmas used by the kernel proofs are covered transitively:
`collectAxioms` follows every dependency of the audited declarations.
-/

open Lean Elab Command

namespace FoldBounds.FrontierKernelAudit

/-- The only dependencies accepted by the audit. -/
def allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]

/-- Audit one declaration; returns the dependencies that are not allowed. -/
def offending (name : Name) : CoreM (Array Name) := do
  let deps ← collectAxioms name
  return deps.filter (fun a => !allowed.contains a)

/-- `#assert_standard_axioms d₁ d₂ …` fails unless every listed declaration
depends only on the standard axioms. -/
elab "#assert_standard_axioms " ids:ident+ : command => do
  for id in ids do
    let name ← liftCoreM <| realizeGlobalConstNoOverloadWithInfo id
    let bad ← liftCoreM <| offending name
    unless bad.isEmpty do
      throwError "{name} depends on non-standard axioms: {bad}"
  logInfo m!"standard axioms only: {ids.size} listed declarations"

/-- `#assert_modules_standard_axioms M₁ M₂ …` audits every theorem and
definition declared in the listed modules. -/
elab "#assert_modules_standard_axioms " mods:ident+ : command => do
  let env ← getEnv
  for m in mods do
    let some idx := env.getModuleIdx? m.getId
      | throwError "unknown module {m.getId}"
    let mut checked := 0
    for (name, info) in env.constants.map₁.toList do
      if env.getModuleIdxFor? name == some idx then
        match info with
        | .thmInfo _ | .defnInfo _ | .opaqueInfo _ =>
            let bad ← liftCoreM <| offending name
            unless bad.isEmpty do
              throwError "{name} depends on non-standard axioms: {bad}"
            checked := checked + 1
        | .axiomInfo _ => throwError "{name} is declared as an axiom"
        | _ => pure ()
    if checked == 0 then
      throwError "module {m.getId} contributed no declarations to the audit"
    logInfo m!"standard axioms only: {checked} declarations of {m.getId}"

end FoldBounds.FrontierKernelAudit

open FoldBounds.FrontierKernelAudit

#assert_modules_standard_axioms FoldBounds.FrontierKernel FoldBounds.FrontierKernelSound
  FoldBounds.FrontierKernelTests

#assert_standard_axioms
  FoldBounds.Kernel.verify_sound
  FoldBounds.Kernel.verifyBytes_sound
  FoldBounds.Kernel.verifyBytes_exact
  FoldBounds.Kernel.certifiesExact_sound
  FoldBounds.Kernel.certifiesInterval_sound
  FoldBounds.Kernel.encEntry_le_basic
  FoldBounds.Kernel.runTree_sound
  FoldBounds.Kernel.enumGo_sound
  FoldBounds.Kernel.encEntry_sound
  FoldBounds.Kernel.boundEntry_ok
  FoldBounds.FrontierKernelTests.hhhh_probability_exact
  FoldBounds.FrontierKernelTests.h8_probability_exact
  FoldBounds.FrontierKernelTests.h8_probability_interval
