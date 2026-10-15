module

public import FoldBounds.Executable

/-!
A concrete executable 2D HP lattice model. Directions after a fixed
first bond form a finite word. Its subtype enforces self-avoidance.
All mass arithmetic is exact rational. Upper-cap checking is still
enumerative and no equivalence to the Python model is yet claimed.
-/
@[expose] public section

namespace FoldBounds
namespace LatticeModel

abbrev Site := ℤ × ℤ

def origin : Site := (0, 0)
def first : Site := (1, 0)

def advance (p : Site) (d : Fin 4) : Site :=
  if d.val = 0 then (p.1 + 1, p.2)
  else if d.val = 1 then (p.1 - 1, p.2)
  else if d.val = 2 then (p.1, p.2 + 1)
  else (p.1, p.2 - 1)

def trace (p : Site) : List (Fin 4) → List Site
  | [] => []
  | d :: ds =>
      let next := advance p d
      next :: trace next ds

def walkFrom (ds : List (Fin 4)) : List Site :=
  [origin, first] ++ trace first ds

theorem trace_length (p : Site) (ds : List (Fin 4)) :
    (trace p ds).length = ds.length := by
  induction ds generalizing p with
  | nil => rfl
  | cons d ds ih => simp [trace, ih]

theorem walkFrom_length (ds : List (Fin 4)) :
    (walkFrom ds).length = ds.length + 2 := by
  simp [walkFrom, trace_length]

abbrev Conformation (steps : ℕ) :=
  { directions : Fin steps → Fin 4 //
    (walkFrom (List.ofFn directions)).Nodup }

theorem conformation_nodup {steps : ℕ} (c : Conformation steps) :
    (walkFrom (List.ofFn c.1)).Nodup := c.2

theorem conformation_length {steps : ℕ} (c : Conformation steps) :
    (walkFrom (List.ofFn c.1)).length = steps + 2 := by
  rw [walkFrom_length]
  simp

def adjacent (p r : Site) : Bool :=
  decide ((p.1 + 1 = r.1 ∧ p.2 = r.2) ∨
          (p.1 - 1 = r.1 ∧ p.2 = r.2) ∨
          (p.1 = r.1 ∧ p.2 + 1 = r.2) ∨
          (p.1 = r.1 ∧ p.2 - 1 = r.2))

def isHydrophobic (seq : List Bool) (i : ℕ) : Bool :=
  (seq[i]?).getD false

def contactCount (seq : List Bool) (walk : List Site) : ℕ :=
  ∑ i ∈ Finset.range seq.length,
    ∑ j ∈ Finset.range seq.length,
      if decide (i + 1 < j) &&
         isHydrophobic seq i &&
         isHydrophobic seq j &&
         adjacent ((walk[i]?).getD origin) ((walk[j]?).getD origin)
      then 1 else 0

def weightQ (steps : ℕ) (q : ℚ) (seq : Fin (steps + 2) → Bool)
    (c : Conformation steps) : ℚ :=
  q ^ contactCount (List.ofFn seq) (walkFrom (List.ofFn c.1))

def contactEvent (steps : ℕ) (seq : Fin (steps + 2) → Bool)
    (threshold : ℕ) : Finset (Conformation steps) :=
  Finset.univ.filter (fun c =>
    threshold ≤ contactCount (List.ofFn seq) (walkFrom (List.ofFn c.1)))

def checkHP (steps : ℕ) (q : ℚ) (seq : Fin (steps + 2) → Bool)
    (threshold : ℕ) (cert : Executable.RawCertificate (Conformation steps)) :
    Option Executable.RationalBounds :=
  if 1 ≤ q then
    Executable.check (weightQ steps q seq) (contactEvent steps seq threshold) cert
  else none

theorem checkHP_coupled_sound (steps : ℕ) [Nonempty (Conformation steps)]
    (q : ℚ) (seq : Fin (steps + 2) → Bool) (threshold : ℕ)
    (cert : Executable.RawCertificate (Conformation steps))
    (bounds : Executable.RationalBounds)
    (h : checkHP steps q seq threshold cert = some bounds) :
    ((Executable.massQ (weightQ steps q seq) (contactEvent steps seq threshold) : ℝ) /
      ((Executable.massQ (weightQ steps q seq) (contactEvent steps seq threshold) : ℝ) +
       (Executable.massQ (weightQ steps q seq)
         (contactEvent steps seq threshold)ᶜ : ℝ))) ∈
      Set.Icc
        ((bounds.eventLower / (bounds.eventLower + bounds.otherUpper)) : ℝ)
        ((bounds.eventUpper / (bounds.eventUpper + bounds.otherLower)) : ℝ) := by
  by_cases hq : (1 : ℚ) ≤ q
  · have checked :
        Executable.check (weightQ steps q seq) (contactEvent steps seq threshold) cert =
          some bounds := by simpa [checkHP, hq] using h
    exact Executable.check_coupled_probability_sound _ _ _ _ checked
  · simp [checkHP, hq] at h

end LatticeModel
end FoldBounds
