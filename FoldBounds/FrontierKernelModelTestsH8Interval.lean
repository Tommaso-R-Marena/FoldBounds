import FoldBounds.FrontierKernelSound

/-!
Certified enclosure for `HHHHHHHH` at `q = 2`, `t = 2`, plus a numeric
sanity check that the exact value lies in the enclosure.
-/

open FoldBounds FoldBounds.Kernel FoldBounds.LatticeModel

namespace FoldBounds.FrontierKernelModelTests

def allH8' : Fin 8 → Bool := fun _ => true

theorem h8_probability_interval :
    Executable.massQ (weightQ 6 2 allH8') (contactEvent 6 allH8' 2) /
      Executable.massQ (weightQ 6 2 allH8') Finset.univ ∈
        Set.Icc (196/733 : ℚ) (989/1028) :=
  certifiesInterval_sound allH8'
    ("FOLDBOUNDS-KERNEL 1\nsequence HHHHHHHH\nquery 2 2\nquery 3/2 1\ntree NNNB-BB-NBXB-NBX-B-NNB-BXN-XX-NBBB--NNB-XBN-X-X-NBB-B\nend\n").toUTF8
    2 2 (196/733) (989/1028) (by decide +kernel)

example : (196/733 : ℚ) ≤ 72/139 ∧
    (72/139 : ℚ) ≤ 989/1028 := by norm_num

end FoldBounds.FrontierKernelModelTests
