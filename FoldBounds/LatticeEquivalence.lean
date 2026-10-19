module

public import FoldBounds.LatticeModel

/-!
# Representation equivalence for anchored square-lattice conformations

`LatticeModel.Conformation steps` stores a conformation as a word of
`steps` absolute directions following the fixed first bond
`(0,0) → (1,0)`.  This file introduces the *geometric* type
`CanonicalSAW steps` of anchored self-avoiding site paths and proves an
explicit, computable bijection between the two representations
(`conformationEquivSAW`).  Consequently partition sums over either
representation coincide for *every* weight function, with the correct
multiplicity (no symmetry quotient is taken on either side).

It also proves constructively that the conformation space is nonempty
(the straight rod), removing the need for an external `[Nonempty]`
assumption downstream.
-/
@[expose] public section

namespace FoldBounds
namespace LatticeModel

/-! ### Endpoint of a decoded walk and appending directions -/

/-- The site reached from `s` after following the directions `ds`. -/
def endFrom (s : Site) : List (Fin 4) → Site
  | [] => s
  | d :: ds => endFrom (advance s d) ds

theorem trace_append (s : Site) (ds es : List (Fin 4)) :
    trace s (ds ++ es) = trace s ds ++ trace (endFrom s ds) es := by
  induction ds generalizing s with
  | nil => rfl
  | cons d ds ih => simp [trace, endFrom, ih]

theorem endFrom_append (s : Site) (ds es : List (Fin 4)) :
    endFrom s (ds ++ es) = endFrom (endFrom s ds) es := by
  induction ds generalizing s with
  | nil => rfl
  | cons d ds ih => simp [endFrom, ih]

theorem endFrom_mem (s : Site) (ds : List (Fin 4)) :
    endFrom s ds ∈ s :: trace s ds := by
  induction ds generalizing s with
  | nil => simp [endFrom]
  | cons d ds ih =>
      have := ih (advance s d)
      simp only [endFrom, trace, List.mem_cons] at this ⊢
      tauto

/-- Appending one direction appends exactly one new site to the walk. -/
theorem walkFrom_append_single (ds : List (Fin 4)) (d : Fin 4) :
    walkFrom (ds ++ [d]) = walkFrom ds ++ [advance (endFrom first ds) d] := by
  simp [walkFrom, trace_append, trace]

theorem walkFrom_append (ds es : List (Fin 4)) :
    walkFrom (ds ++ es) = walkFrom ds ++ trace (endFrom first ds) es := by
  simp [walkFrom, trace_append]

theorem endFrom_mem_walkFrom (ds : List (Fin 4)) :
    endFrom first ds ∈ walkFrom ds := by
  have := endFrom_mem first ds
  simp only [walkFrom, List.cons_append, List.nil_append, List.mem_cons] at this ⊢
  tauto

/-- Decoding is prefix-monotone: a prefix word decodes to a prefix walk. -/
theorem walkFrom_prefix {ds es : List (Fin 4)} (h : ds <+: es) :
    walkFrom ds <+: walkFrom es := by
  obtain ⟨r, rfl⟩ := h
  rw [walkFrom_append]
  exact List.prefix_append _ _

/-! ### Decoding geometric paths back to directions -/

/-- The direction of a single lattice step `a → b` (meaningful when the
two sites are adjacent; it is the unique such direction). -/
def stepDir (a b : Site) : Fin 4 :=
  if b = (a.1 + 1, a.2) then 0
  else if b = (a.1 - 1, a.2) then 1
  else if b = (a.1, a.2 + 1) then 2
  else 3

theorem advance_stepDir {a b : Site} (h : adjacent a b = true) :
    advance a (stepDir a b) = b := by
  obtain ⟨a1, a2⟩ := a
  obtain ⟨b1, b2⟩ := b
  simp only [adjacent, decide_eq_true_eq] at h
  rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ | ⟨h1, h2⟩ <;>
    subst h1 h2
  · simp [stepDir, advance]
  · have : a1 - 1 ≠ a1 + 1 := by omega
    simp [stepDir, advance, this]
  · simp [stepDir, advance]
  · have : a2 - 1 ≠ a2 + 1 := by omega
    simp [stepDir, advance, this]

theorem stepDir_advance (a : Site) (d : Fin 4) : stepDir a (advance a d) = d := by
  apply advance_injective a
  exact advance_stepDir (advance_adjacent a d)

/-- Read off the directions of the successive steps of a site path. -/
def decodeSites : List Site → List (Fin 4)
  | a :: b :: rest => stepDir a b :: decodeSites (b :: rest)
  | _ => []

theorem decodeSites_trace (s : Site) (ds : List (Fin 4)) :
    decodeSites (s :: trace s ds) = ds := by
  induction ds generalizing s with
  | nil => rfl
  | cons d ds ih =>
      simp only [trace, decodeSites, stepDir_advance]
      rw [ih]

theorem trace_decodeSites (s : Site) (rest : List Site)
    (h : adjacentChain (s :: rest) = true) :
    trace s (decodeSites (s :: rest)) = rest := by
  induction rest generalizing s with
  | nil => rfl
  | cons b rest ih =>
      simp only [adjacentChain, Bool.and_eq_true] at h
      simp only [decodeSites, trace, advance_stepDir h.1]
      rw [ih b h.2]

/-- Directions of an anchored path (the first, fixed bond is skipped). -/
def pathDirections (path : List Site) : List (Fin 4) := decodeSites (path.drop 1)

/-! ### Canonical anchored self-avoiding paths -/

/-- Canonical anchored self-avoiding square-lattice site paths with
`steps + 2` sites (`steps` free directions after the fixed first bond):
the path starts `(0,0), (1,0)`, consecutive sites are grid neighbours and
no site repeats.  This is exactly the Python state convention: reflection
is *not* quotiented out. -/
def CanonicalSAW (steps : ℕ) : Type :=
  { path : List Site //
    path.length = steps + 2 ∧ path.take 2 = [origin, first] ∧
    adjacentChain path = true ∧ path.Nodup }

instance (steps : ℕ) : DecidableEq (CanonicalSAW steps) :=
  inferInstanceAs (DecidableEq (Subtype _))

theorem canonical_shape {path : List Site} (h : path.take 2 = [origin, first]) :
    ∃ rest, path = origin :: first :: rest := by
  rcases path with _ | ⟨a, _ | ⟨b, rest⟩⟩
  · simp at h
  · simp at h
  · simp only [List.take_succ_cons, List.take_zero, List.cons.injEq, and_true] at h
    exact ⟨rest, by rw [h.1, h.2]⟩

theorem walkFrom_pathDirections {path : List Site} (h2 : path.take 2 = [origin, first])
    (hadj : adjacentChain path = true) :
    walkFrom (pathDirections path) = path := by
  obtain ⟨rest, rfl⟩ := canonical_shape h2
  simp only [adjacentChain, Bool.and_eq_true] at hadj
  simp [walkFrom, pathDirections, trace_decodeSites first rest hadj.2]

theorem pathDirections_walkFrom (ds : List (Fin 4)) :
    pathDirections (walkFrom ds) = ds := by
  simp [pathDirections, walkFrom, decodeSites_trace]

theorem walkFrom_take_two (ds : List (Fin 4)) :
    (walkFrom ds).take 2 = [origin, first] := by
  simp [walkFrom]

/-- Encoding: decode a direction word into its geometric path. -/
def toSAW {steps : ℕ} (c : Conformation steps) : CanonicalSAW steps :=
  ⟨walkFrom (List.ofFn c.1),
    conformation_length c, walkFrom_take_two _, conformation_adjacentChain c, c.2⟩

theorem pathDirections_length {steps : ℕ} (P : CanonicalSAW steps) :
    (pathDirections P.1).length = steps := by
  have h := congrArg List.length (walkFrom_pathDirections P.2.2.1 P.2.2.2.1)
  rw [walkFrom_length, P.2.1] at h
  omega

/-- Decoding: read the direction word of a geometric path. -/
def ofSAW {steps : ℕ} (P : CanonicalSAW steps) : Conformation steps :=
  ⟨fun i => (pathDirections P.1)[i.1]'(by rw [pathDirections_length]; exact i.2), by
    have h : (List.ofFn fun i : Fin steps =>
        (pathDirections P.1)[i.1]'(by rw [pathDirections_length]; exact i.2)) =
        pathDirections P.1 := by
      apply List.ext_getElem <;> simp [pathDirections_length]
    rw [h, walkFrom_pathDirections P.2.2.1 P.2.2.2.1]
    exact P.2.2.2.2⟩

theorem ofFn_ofSAW {steps : ℕ} (P : CanonicalSAW steps) :
    List.ofFn (ofSAW P).1 = pathDirections P.1 := by
  apply List.ext_getElem <;> simp [ofSAW, pathDirections_length]

/-- **Representation theorem.** Direction-word conformations and canonical
anchored self-avoiding site paths are in explicit computable bijection.
The bijection fixes the anchor `(0,0) → (1,0)`, preserves length, chain
order, adjacency and self-avoidance, and identifies no two distinct states
(multiplicities are preserved exactly). -/
def conformationEquivSAW (steps : ℕ) : Conformation steps ≃ CanonicalSAW steps where
  toFun := toSAW
  invFun := ofSAW
  left_inv c := by
    apply Subtype.ext
    apply List.ofFn_injective
    rw [ofFn_ofSAW]
    exact pathDirections_walkFrom _
  right_inv P := by
    apply Subtype.ext
    change walkFrom (List.ofFn (ofSAW P).1) = P.1
    rw [ofFn_ofSAW]
    exact walkFrom_pathDirections P.2.2.1 P.2.2.2.1

theorem conformationEquivSAW_apply {steps : ℕ} (c : Conformation steps) :
    (conformationEquivSAW steps c).1 = walkFrom (List.ofFn c.1) := rfl

/-- Injectivity of the decoding map `Conformation → CanonicalSAW`. -/
theorem toSAW_injective (steps : ℕ) :
    Function.Injective (toSAW : Conformation steps → CanonicalSAW steps) :=
  (conformationEquivSAW steps).injective

/-- Surjectivity: every canonical anchored self-avoiding path arises
from a (unique) valid direction word. -/
theorem toSAW_surjective (steps : ℕ) :
    Function.Surjective (toSAW : Conformation steps → CanonicalSAW steps) :=
  (conformationEquivSAW steps).surjective

noncomputable instance (steps : ℕ) : Fintype (CanonicalSAW steps) :=
  Fintype.ofEquiv _ (conformationEquivSAW steps)

/-- **Partition-sum preservation** for an arbitrary weight on geometric
paths: summing over direction words (decoded) equals summing over canonical
paths, with identical multiplicities. -/
theorem sum_conformation_eq_sum_saw {M : Type*} [AddCommMonoid M] (steps : ℕ)
    (w : CanonicalSAW steps → M) :
    ∑ c : Conformation steps, w (toSAW c) = ∑ P : CanonicalSAW steps, w P :=
  Fintype.sum_equiv (conformationEquivSAW steps) _ _ (fun _ => rfl)

/-- The number of valid direction words equals the number of canonical paths. -/
theorem card_conformation_eq_card_saw (steps : ℕ) :
    Fintype.card (Conformation steps) = Fintype.card (CanonicalSAW steps) :=
  Fintype.card_congr (conformationEquivSAW steps)

/-- Contact weight of a geometric path (the same `contactCount` scorer). -/
def sawWeightQ (steps : ℕ) (q : ℚ) (seq : Fin (steps + 2) → Bool)
    (P : CanonicalSAW steps) : ℚ :=
  q ^ contactCount (List.ofFn seq) P.1

/-- Specialisation to FoldBounds' exact rational contact weights: the
partition function over direction words equals the one over geometric
canonical paths. -/
theorem massQ_weightQ_eq_sum_saw (steps : ℕ) (q : ℚ) (seq : Fin (steps + 2) → Bool) :
    Executable.massQ (weightQ steps q seq) Finset.univ =
      ∑ P : CanonicalSAW steps, sawWeightQ steps q seq P :=
  sum_conformation_eq_sum_saw steps (sawWeightQ steps q seq)

/-! ### Constructive nonemptiness: the straight rod -/

theorem trace_straight (a b : ℤ) (k : ℕ) :
    trace (a, b) (List.replicate k 0) =
      (List.range k).map (fun i : ℕ => (a + (i : ℤ) + 1, b)) := by
  induction k generalizing a with
  | zero => rfl
  | succ k ih =>
      rw [List.replicate_succ, trace, List.range_succ_eq_map]
      simp only [advance, Fin.isValue, Fin.val_zero, if_true, List.map_cons, List.map_map]
      rw [ih]
      simp only [Nat.cast_zero, add_zero, List.cons.injEq, true_and]
      apply List.map_congr_left
      intro i _
      simp only [Function.comp, Nat.cast_succ]
      ext <;> simp
      ring

theorem walkFrom_straight_nodup (k : ℕ) : (walkFrom (List.replicate k 0)).Nodup := by
  rw [walkFrom, first, trace_straight]
  simp only [origin, List.cons_append, List.nil_append]
  refine List.nodup_cons.2 ⟨?_, List.nodup_cons.2 ⟨?_, ?_⟩⟩
  · simp only [List.mem_cons, List.mem_map, List.mem_range, not_or]
    refine ⟨by simp, ?_⟩
    rintro ⟨i, _, h⟩
    simp only [Prod.mk.injEq] at h
    omega
  · simp only [List.mem_map, List.mem_range, not_exists, not_and]
    intro i _ h
    simp only [Prod.mk.injEq] at h
    omega
  · refine List.Nodup.map ?_ List.nodup_range
    intro i j h
    simp only [Prod.mk.injEq] at h
    omega

/-- The straight-rod conformation (all directions `+x`). -/
def straight (steps : ℕ) : Conformation steps :=
  ⟨fun _ => 0, by
    have : List.ofFn (fun _ : Fin steps => (0 : Fin 4)) = List.replicate steps 0 := by
      simp
    rw [this]
    exact walkFrom_straight_nodup steps⟩

/-- The conformation space is nonempty for every chain length. -/
instance conformation_nonempty (steps : ℕ) : Nonempty (Conformation steps) :=
  ⟨straight steps⟩

/-- `checkHP_coupled_sound` with the `[Nonempty]` premise discharged
from the model itself. -/
theorem checkHP_coupled_sound_unconditional (steps : ℕ)
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
        ((bounds.eventUpper / (bounds.eventUpper + bounds.otherLower)) : ℝ) :=
  checkHP_coupled_sound steps q seq threshold cert bounds h

end LatticeModel
end FoldBounds
