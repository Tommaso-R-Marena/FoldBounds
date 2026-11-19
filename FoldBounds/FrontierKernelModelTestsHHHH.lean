import FoldBounds.FrontierKernelSound

/-!
Exact `HHHH` probability bridge: `P(k ≥ 1) = 4/11` at `q = 2`.
-/

open FoldBounds FoldBounds.Kernel FoldBounds.LatticeModel

namespace FoldBounds.FrontierKernelModelTests

def allH4 : Fin 4 → Bool := fun _ => true

theorem hhhh_probability_exact :
    Executable.massQ (weightQ 2 2 allH4) (contactEvent 2 allH4 1) /
      Executable.massQ (weightQ 2 2 allH4) Finset.univ = 4 / 11 :=
  certifiesExact_sound allH4
    ("FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1\ntree X\nend\n").toUTF8
    2 1 (4 / 11) (by decide +kernel)

end FoldBounds.FrontierKernelModelTests
