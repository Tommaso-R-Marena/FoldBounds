import FoldBounds.Executable

open Finset
open FoldBounds.Executable

private def demoWeights : Fin 3 → ℚ := ![1, 2, 4]
private def demoEvent : Finset (Fin 3) := {0, 1}
private def demoCertificate : RawCertificate (Fin 3) :=
  ⟨{0}, {2}, fun _ => 3⟩
private def invalidCertificate : RawCertificate (Fin 3) :=
  ⟨{0}, {2}, fun _ => 1⟩
private def duplicateClassCertificate : RawCertificate (Fin 3) :=
  ⟨{0, 2}, {2}, fun _ => 10⟩
private def impossibleWeightCertificate : RawCertificate (Fin 3) :=
  ⟨{0}, {2}, fun _ => 0⟩

example : accepts demoWeights demoEvent demoCertificate = true := by decide
example : check demoWeights demoEvent demoCertificate =
    some (endpoints demoWeights demoEvent demoCertificate) := by decide
example : accepts demoWeights demoEvent invalidCertificate = false := by decide
example : check demoWeights demoEvent invalidCertificate = none := by decide
example : check demoWeights demoEvent duplicateClassCertificate = none := by decide
example : check demoWeights demoEvent impossibleWeightCertificate = none := by decide

/-- Concrete, kernel-checked rational finite-universe example. -/
example : (massQ demoWeights demoEvent : ℚ) = 3 := by decide
example : (massQ demoWeights demoEventᶜ : ℚ) = 4 := by decide
