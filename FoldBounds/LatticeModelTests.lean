import FoldBounds.LatticeModel

open FoldBounds
open FoldBounds.LatticeModel

example : walkFrom [] = [origin, first] := by decide
example : walkFrom [⟨2, by decide⟩, ⟨1, by decide⟩] =
    [(0,0),(1,0),(1,1),(0,1)] := by decide
example : contactCount [true,true,true,true]
    [(0,0),(1,0),(1,1),(0,1)] = 1 := by decide
example : contactCount [true,false,true,true]
    [(0,0),(1,0),(1,1),(0,1)] = 1 := by decide
example : contactCount [false,true,true,true]
    [(0,0),(1,0),(1,1),(0,1)] = 0 := by decide
example : ¬ (walkFrom [⟨1, by decide⟩]).Nodup := by decide

private def allH : Fin 2 → Bool := fun _ => true
private def demoCertificate : Executable.RawCertificate (Conformation 0) :=
  ⟨Finset.univ, ∅, fun _ => 0⟩
#eval checkHP 0 2 allH 0 demoCertificate
example : checkHP 0 (1/2) allH 0 demoCertificate = none := by
  norm_num [checkHP]

private def wrongEvent : Executable.RawCertificate (Conformation 0) :=
  ⟨Finset.univ, ∅, fun _ => 2⟩
example : checkHP 0 2 allH 1 wrongEvent = none := by decide
