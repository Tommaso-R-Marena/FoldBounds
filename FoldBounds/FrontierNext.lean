module

public import FoldBounds.VerifiedMemoryTransfer
public import FoldBounds.FrontierSoundness
public import FoldBounds.VerifiedFamilies

/-!
# The next-generation frontier checker (format version 2)

`NTree` extends the prefix-trie certificates of `FrontierCheck` with two
new leaf kinds; every old leaf kind is still accepted with exactly its old
meaning, so the version-2 checker is a strict superset of version 1.

* `mem k ws` — a valid prefix bounded by the **verified finite-memory
  transfer relaxation** of memory `k` (`MemTransfer`), whose final-state
  envelope uses
  - the exact cross score,
  - the history-certified unplaced–unplaced contact window `[ffLo, ffHi]`,
  - the degree cap `known + futureCap`,
  - the **geometric perimeter cap** `geomCap` of the prefix, and the
    endpoint- and history-conditioned geometric cap of each final state;
  together with optional exact witnesses (upper-bound tightening) and a
  **disjoint** family + witness lower bound: witnesses that are not members
  of the certified diagonal family are *added* to the family mass, with a
  kernel-checked disjointness proof.
* `enum` — a valid prefix whose subtree is evaluated exactly by local
  depth-first search inside the checker (the certificate stays one line,
  the verifier pays the subtree size).
* `legacy tr ws` — the old `B`/`T` leaf (`FrontierCheck.boundEnc` +
  `famRaise`), unchanged.
* `exact`, `dead`, `node` — unchanged.

The flagship theorem `check_sound` has exactly the shape of
`FrontierCheck.check_sound`.
-/
@[expose] public section

namespace FoldBounds
namespace FrontierNext

open LatticeModel Frontier Transfer Geometry MemTransfer FrontierCheck Families Finset

/-! ### Certificate format -/

/-- Version-2 untrusted prefix-trie certificate. -/
inductive NTree where
  | exact : NTree
  | dead : NTree
  | enum : NTree
  | legacy (transfer : Bool) (witnesses : List (List (Fin 4))) : NTree
  | mem (memory : ℕ) (witnesses : List (List (Fin 4))) : NTree
  | node (t0 t1 t2 t3 : NTree) : NTree
  deriving Repr, DecidableEq

/-! ### Exhaustive local evaluation -/

/-- Exact local evaluation of a subtree by depth-first search (`r` = number
of remaining free directions). -/
def enumEnc (S : List Bool) (q : ℚ) (thr : ℕ) : ℕ → List (Fin 4) → Enclosure
  | 0, p => if (walkFrom p).Nodup then exactEnc S q thr p else 0
  | r + 1, p =>
      if (walkFrom p).Nodup then
        enumEnc S q thr r (p ++ [0]) + enumEnc S q thr r (p ++ [1]) +
          enumEnc S q thr r (p ++ [2]) + enumEnc S q thr r (p ++ [3])
      else 0

/-! ### Memory-transfer leaves -/

/-- Constant contact cap of a prefix: degree cap and prefix geometric cap. -/
def kcN (steps : ℕ) (S : List Bool) (p : List (Fin 4)) : ℕ :=
  min (kcap true S p) (geomCap (steps + 2) (walkFrom p))

/-- Upper contact exponent of a final memory state. -/
def mExp (steps : ℕ) (S : List Bool) (p : List (Fin 4)) (st : MState) : ℕ :=
  min (min (knownContacts S p + st.cross + st.ffHi) (kcN steps S p))
    (geomCap (steps + 2) (knownSites p st))

/-- Certified lower contact exponent of a final memory state. -/
def mLow (S : List Bool) (p : List (Fin 4)) (st : MState) : ℕ :=
  knownContacts S p + st.cross + st.ffLo

def hZ (steps : ℕ) (S : List Bool) (q : ℚ) (p : List (Fin 4)) (st : MState) : ℚ :=
  q ^ mExp steps S p st

def hA (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) (st : MState) : ℚ :=
  if thr ≤ mExp steps S p st then q ^ mExp steps S p st else 0

def hN (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) (st : MState) : ℚ :=
  if mLow S p st < thr then q ^ mExp steps S p st else 0

/-- Counting × constant-cap bounds (`3^(n-m) q^{kcN}`). -/
def capZ (steps : ℕ) (S : List Bool) (q : ℚ) (p : List (Fin 4)) : ℚ :=
  3 ^ (steps - p.length) * q ^ kcN steps S p

def capA (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) : ℚ :=
  if thr ≤ kcN steps S p then capZ steps S q p else 0

def capN (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) : ℚ :=
  if knownContacts S p < thr then capZ steps S q p else 0

/-- Envelope uppers of a memory-`k` leaf. -/
def memUppers (k steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) :
    ℚ × ℚ × ℚ :=
  let tbl := memTable steps S k p
  (min (capZ steps S q p) (wsum (hZ steps S q p) tbl),
   min (capA steps S q thr p) (wsum (hA steps S q thr p) tbl),
   min (capN steps S q thr p) (wsum (hN steps S q thr p) tbl))

/-- First available diagonal family at `p`. -/
def famPick (p : List (Fin 4)) : Option (Bool × Bool) :=
  if diagOK true true p then some (true, true)
  else if diagOK true false p then some (true, false)
  else if diagOK false true p then some (false, true)
  else if diagOK false false p then some (false, false)
  else none

/-- Is the witness suffix `w` a member of the diagonal family `ab`? -/
def inFam (ab : Bool × Bool) (w : List (Fin 4)) : Bool :=
  w.all (fun d => decide (d = dirX ab.1) || decide (d = dirY ab.2))

/-- A certified constructive family available at a prefix. -/
inductive FamChoice where
  | pd (m : Fin 4)
  | diag (a b : Bool)
  deriving Repr, DecidableEq

/-- Choose a family: a partially directed family if the endpoint is a
strict extreme along some axis, otherwise a diagonal family. -/
def famChoose (p : List (Fin 4)) : Option FamChoice :=
  if pdOK 0 p then some (.pd 0) else if pdOK 1 p then some (.pd 1)
  else if pdOK 2 p then some (.pd 2) else if pdOK 3 p then some (.pd 3)
  else (famPick p).map (fun ab => .diag ab.1 ab.2)

/-- Exact family size for `r` remaining directions (computed by recurrence). -/
def famCount : FamChoice → ℕ → ℕ
  | .pd _, r => pdCount r
  | .diag _ _, r => 2 ^ r

/-- Executable family-membership test of a witness suffix. -/
def famTest : FamChoice → List (Fin 4) → Bool
  | .pd m, w => pdTest m 0 w
  | .diag a b, w => inFam (a, b) w

/-- The family as a (semantic) finite set of complete words. -/
def famSet : FamChoice → List (Fin 4) → ℕ → Finset (List (Fin 4))
  | .pd m, p, r => pdFamily m p r
  | .diag a b, p, r => family a b p r

/-- Lower endpoints: certified family mass plus the witnesses *outside*
the family (disjoint), or the plain witness sums when no family applies. -/
def memLowers (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4))
    (ws : List (List (Fin 4))) : ℚ × ℚ × ℚ :=
  match famChoose p with
  | none => ((ws.map (witZ S q p)).sum, (ws.map (witA S q thr p)).sum,
      (ws.map (witN S q thr p)).sum)
  | some ch =>
      let F : ℚ := (famCount ch (steps - p.length) : ℚ) * q ^ knownContacts S p
      let rest := ws.filter (fun w => !famTest ch w)
      (F + (rest.map (witZ S q p)).sum,
       if thr ≤ knownContacts S p then F + (rest.map (witA S q thr p)).sum
         else (ws.map (witA S q thr p)).sum,
       if kcN steps S p < thr then F + (rest.map (witN S q thr p)).sum
         else (ws.map (witN S q thr p)).sum)

/-- Enclosure of a memory-`k` leaf with witnesses `ws`. -/
def memEnc (k : ℕ) (ws : List (List (Fin 4))) (steps : ℕ) (S : List Bool) (q : ℚ)
    (thr : ℕ) (p : List (Fin 4)) : Enclosure :=
  let U := memUppers k steps S q thr p
  let L := memLowers steps S q thr p ws
  let wst := fun w => gstate S k p (p ++ w)
  ⟨L.1,
   U.1 - (ws.map (fun w => hZ steps S q p (wst w) - witZ S q p w)).sum,
   L.2.1,
   U.2.1 - (ws.map (fun w => hA steps S q thr p (wst w) - witA S q thr p w)).sum,
   L.2.2,
   U.2.2 - (ws.map (fun w => hN steps S q thr p (wst w) - witN S q thr p w)).sum⟩

/-! ### The checker -/

/-- Recursive version-2 trie checker. -/
def checkNode (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) :
    List (Fin 4) → NTree → Option Enclosure
  | p, .exact =>
      if p.length = steps ∧ (walkFrom p).Nodup then some (exactEnc S q thr p) else none
  | p, .dead => if deadOK steps p then some 0 else none
  | p, .enum =>
      if p.length ≤ steps then some (enumEnc S q thr (steps - p.length) p) else none
  | p, .legacy tr ws =>
      if p.length ≤ steps ∧ (walkFrom p).Nodup ∧ witnessesOK steps p ws then
        some (famRaise tr steps S q thr p (boundEnc tr ws steps S q thr p))
      else none
  | p, .mem k ws =>
      if p.length ≤ steps ∧ (walkFrom p).Nodup ∧ witnessesOK steps p ws then
        some (memEnc k ws steps S q thr p)
      else none
  | p, .node t0 t1 t2 t3 =>
      if p.length < steps ∧ (walkFrom p).Nodup then
        match checkNode steps S q thr (p ++ [0]) t0, checkNode steps S q thr (p ++ [1]) t1,
            checkNode steps S q thr (p ++ [2]) t2, checkNode steps S q thr (p ++ [3]) t3 with
        | some a, some b, some c, some d => some (a + b + c + d)
        | _, _, _, _ => none
      else none

/-- **The version-2 frontier checker.** -/
def check (S : List Bool) (q : ℚ) (thr : ℕ) (t : NTree) : Option FrontierResult :=
  if 2 ≤ S.length ∧ 1 ≤ q then (checkNode (S.length - 2) S q thr [] t).map finalize else none

/-! ### Soundness: exhaustive leaves -/

theorem enumEnc_sound (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (q : ℚ)
    (hq : 1 ≤ q) (thr : ℕ) :
    ∀ r (p : List (Fin 4)), p.length + r = steps →
      Encloses steps S q thr p (enumEnc S q thr r p) := by
  intro r
  induction r with
  | zero =>
      intro p hp
      simp only [enumEnc]
      split_ifs with hn
      · exact exactEnc_sound steps S hS q hq thr p (by omega) hn
      · exact masses_of_empty (ext_eq_empty_of_not_nodup hn)
  | succ r ih =>
      intro p hp
      simp only [enumEnc]
      split_ifs with hn
      · exact Encloses_node (by omega) (ih _ (by simp; omega)) (ih _ (by simp; omega))
          (ih _ (by simp; omega)) (ih _ (by simp; omega))
      · exact masses_of_empty (ext_eq_empty_of_not_nodup hn)

/-! ### Soundness: memory leaves -/

theorem famPick_spec {p : List (Fin 4)} {ab : Bool × Bool} (h : famPick p = some ab) :
    diagOK ab.1 ab.2 p = true := by
  unfold famPick at h
  split_ifs at h with h1 h2 h3 h4 <;> cases h <;> assumption

theorem inFam_of_mem_family {a b : Bool} {p w : List (Fin 4)} {r : ℕ}
    (h : p ++ w ∈ family a b p r) : inFam (a, b) w = true := by
  obtain ⟨f, _, hf⟩ := Finset.mem_image.1 h
  have hw : famWord a b f = w := List.append_cancel_left hf
  subst hw
  simp only [inFam, List.all_eq_true, famWord, List.mem_ofFn, Bool.or_eq_true,
    decide_eq_true_eq]
  rintro d ⟨i, rfl⟩
  by_cases h : f i <;> simp [h]

/-- Contact facts for genuine completions under a memory leaf. -/
theorem mem_facts (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (k : ℕ)
    {p l : List (Fin 4)} (hl : l ∈ ext steps p) :
    knownContacts S p ≤ mLow S p (gstate S k p l) ∧
    mLow S p (gstate S k p l) ≤ parityContactCount S (walkFrom l) ∧
    parityContactCount S (walkFrom l) ≤ mExp steps S p (gstate S k p l) ∧
    mExp steps S p (gstate S k p l) ≤ kcN steps S p := by
  have hl' := mem_ext.1 hl
  obtain ⟨b1, b2, b3⟩ := memContacts_bounds steps S hS k hl
  have hdeg := (contacts_le_degree_caps S hl'.2.2 (by omega) hl'.2.1).1
  have hgeo := parityContacts_le_geomCap S hS (mem_ext_valid hl) (K := walkFrom p)
    (fun a ha => (walkFrom_prefix hl'.2.2).subset ha)
  unfold mLow mExp kcN kcap
  simp only [if_true]
  refine ⟨by omega, b1, ?_, ?_⟩
  · exact le_min (le_min b2 (le_min hdeg hgeo)) b3
  · exact (min_le_left _ _).trans (min_le_right _ _)

theorem memEnc_upper_sound (k : ℕ) (ws : List (List (Fin 4))) (steps : ℕ) (S : List Bool)
    (hS : S.length = steps + 2) (q : ℚ) (hq : 1 ≤ q) (thr : ℕ) (p : List (Fin 4))
    (hlen : p.length ≤ steps) (hn : (walkFrom p).Nodup) (hw : witnessesOK steps p ws) :
    massZ steps S q p ≤ (memEnc k ws steps S q thr p).zHi ∧
    massA steps S q thr p ≤ (memEnc k ws steps S q thr p).eHi ∧
    massN steps S q thr p ≤ (memEnc k ws steps S q thr p).nHi := by
  have hq0 : (0 : ℚ) ≤ q := by linarith
  have hpos : ∀ k : ℕ, (0 : ℚ) ≤ q ^ k := fun k => pow_nonneg hq0 k
  set T := ext steps p with hT
  set W := (ws.map (p ++ ·)).toFinset with hWdef
  have hWT : W ⊆ T := wit_subset hw
  set K : List (Fin 4) → ℕ := fun l => parityContactCount S (walkFrom l) with hK
  set g : List (Fin 4) → MState := fun l => gstate S k p l with hg
  have hfacts : ∀ l ∈ T, knownContacts S p ≤ mLow S p (g l) ∧ mLow S p (g l) ≤ K l ∧
      K l ≤ mExp steps S p (g l) ∧ mExp steps S p (g l) ≤ kcN steps S p :=
    fun l hl => mem_facts steps S hS k hl
  let fZ : List (Fin 4) → ℚ := fun l => q ^ K l
  let fA : List (Fin 4) → ℚ := fun l => if thr ≤ K l then q ^ K l else 0
  let fN : List (Fin 4) → ℚ := fun l => if K l < thr then q ^ K l else 0
  let GZ : List (Fin 4) → ℚ := fun l => hZ steps S q p (g l)
  let GA : List (Fin 4) → ℚ := fun l => hA steps S q thr p (g l)
  let GN : List (Fin 4) → ℚ := fun l => hN steps S q thr p (g l)
  have hfZ : ∀ l ∈ T, fZ l ≤ GZ l := fun l hl => qpow_mono hq (hfacts l hl).2.2.1
  have hfA : ∀ l ∈ T, fA l ≤ GA l := by
    intro l hl
    simp only [fA, GA, hA]
    split_ifs with h1 h2
    · exact qpow_mono hq (hfacts l hl).2.2.1
    · exact absurd (h1.trans (hfacts l hl).2.2.1) h2
    · exact hpos _
    · exact le_refl _
  have hfN : ∀ l ∈ T, fN l ≤ GN l := by
    intro l hl
    simp only [fN, GN, hN]
    split_ifs with h1 h2
    · exact qpow_mono hq (hfacts l hl).2.2.1
    · exact absurd (lt_of_le_of_lt (hfacts l hl).2.1 h1) h2
    · exact hpos _
    · exact le_refl _
  have hhZ0 : ∀ st, 0 ≤ hZ steps S q p st := fun st => hpos _
  have hhA0 : ∀ st, 0 ≤ hA steps S q thr p st := fun st => by unfold hA; split_ifs <;> simp [hpos]
  have hhN0 : ∀ st, 0 ≤ hN steps S q thr p st := fun st => by unfold hN; split_ifs <;> simp [hpos]
  -- pointwise constant caps
  have hcZ : ∀ l ∈ T, GZ l ≤ q ^ kcN steps S p := fun l hl => qpow_mono hq (hfacts l hl).2.2.2
  have hcA : ∀ l ∈ T, GA l ≤ if thr ≤ kcN steps S p then q ^ kcN steps S p else 0 := by
    intro l hl
    simp only [GA, hA]
    split_ifs with h1 h2
    · exact qpow_mono hq (hfacts l hl).2.2.2
    · exact absurd (h1.trans (hfacts l hl).2.2.2) h2
    · exact hpos _
    · exact le_refl _
  have hcN : ∀ l ∈ T, GN l ≤ if knownContacts S p < thr then q ^ kcN steps S p else 0 := by
    intro l hl
    simp only [GN, hN]
    split_ifs with h1 h2
    · exact qpow_mono hq (hfacts l hl).2.2.2
    · exact absurd (lt_of_le_of_lt (hfacts l hl).1 h1) h2
    · exact hpos _
    · exact le_refl _
  have hbZ : ∑ l ∈ T, GZ l ≤ capZ steps S q p := sum_le_pow_three_mul GZ _ (hpos _) hcZ
  have hbA : ∑ l ∈ T, GA l ≤ capA steps S q thr p := by
    have := sum_le_pow_three_mul (steps := steps) (p := p) GA
      (if thr ≤ kcN steps S p then q ^ kcN steps S p else 0) (by split_ifs <;> simp [hpos]) hcA
    unfold capA capZ; split_ifs at this ⊢ <;> simpa using this
  have hbN : ∑ l ∈ T, GN l ≤ capN steps S q thr p := by
    have := sum_le_pow_three_mul (steps := steps) (p := p) GN
      (if knownContacts S p < thr then q ^ kcN steps S p else 0)
      (by split_ifs <;> simp [hpos]) hcN
    unfold capN capZ; split_ifs at this ⊢ <;> simpa using this
  have hU : ∑ l ∈ T, GZ l ≤ (memUppers k steps S q thr p).1 ∧
      ∑ l ∈ T, GA l ≤ (memUppers k steps S q thr p).2.1 ∧
      ∑ l ∈ T, GN l ≤ (memUppers k steps S q thr p).2.2 :=
    ⟨le_min hbZ (memCover steps S hS k p hn hlen _ hhZ0),
     le_min hbA (memCover steps S hS k p hn hlen _ hhA0),
     le_min hbN (memCover steps S hS k p hn hlen _ hhN0)⟩
  have hnd := hw.1
  have dZ : (ws.map (fun w => hZ steps S q p (gstate S k p (p ++ w)) - witZ S q p w)).sum =
      ∑ l ∈ W, (GZ l - fZ l) := wit_sum p ws hnd (fun l => GZ l - fZ l)
  have dA : (ws.map (fun w => hA steps S q thr p (gstate S k p (p ++ w)) -
      witA S q thr p w)).sum = ∑ l ∈ W, (GA l - fA l) := wit_sum p ws hnd (fun l => GA l - fA l)
  have dN : (ws.map (fun w => hN steps S q thr p (gstate S k p (p ++ w)) -
      witN S q thr p w)).sum = ∑ l ∈ W, (GN l - fN l) := wit_sum p ws hnd (fun l => GN l - fN l)
  have mZ : massZ steps S q p = ∑ l ∈ T, fZ l := massZ_eq steps S hS q p
  have mA : massA steps S q thr p = ∑ l ∈ T, fA l := massA_eq steps S hS q thr p
  have mN : massN steps S q thr p = ∑ l ∈ T, fN l := massN_eq steps S hS q thr p
  simp only [memEnc]
  rw [dZ, dA, dN, mZ, mA, mN]
  exact ⟨envelope_bound hWT fZ GZ hfZ hU.1, envelope_bound hWT fA GA hfA hU.2.1,
    envelope_bound hWT fN GN hfN hU.2.2⟩

theorem famChoose_spec {steps : ℕ} {p : List (Fin 4)} (hn : (walkFrom p).Nodup)
    (hlen : p.length ≤ steps) {ch : FamChoice} (h : famChoose p = some ch) :
    famSet ch p (steps - p.length) ⊆ ext steps p ∧
    (famSet ch p (steps - p.length)).card = famCount ch (steps - p.length) ∧
    ∀ w, p ++ w ∈ famSet ch p (steps - p.length) → famTest ch w = true := by
  have hpd : ∀ m, pdOK m p = true → famSet (.pd m) p (steps - p.length) ⊆ ext steps p ∧
      (famSet (.pd m) p (steps - p.length)).card = famCount (.pd m) (steps - p.length) ∧
      ∀ w, p ++ w ∈ famSet (.pd m) p (steps - p.length) → famTest (.pd m) w = true :=
    fun m hm => ⟨pdFamily_subset_ext hn hm hlen, card_pdFamily m p _,
      fun w hw => mem_pdFamily_test hw⟩
  unfold famChoose at h
  split_ifs at h with h0 h1 h2 h3
  · cases h; exact hpd 0 h0
  · cases h; exact hpd 1 h1
  · cases h; exact hpd 2 h2
  · cases h; exact hpd 3 h3
  · obtain ⟨ab, hab, rfl⟩ := Option.map_eq_some_iff.1 h
    exact ⟨family_subset_ext hn (famPick_spec hab) hlen, card_family _ _ _ _,
      fun w hw => inFam_of_mem_family hw⟩

/-- Disjoint family + witness lower bound for one nonnegative weight. -/
theorem fam_wit_lower {steps : ℕ} {p : List (Fin 4)} (hn : (walkFrom p).Nodup)
    (hlen : p.length ≤ steps) {ws : List (List (Fin 4))} (hw : witnessesOK steps p ws)
    {ch : FamChoice} (hch : famChoose p = some ch)
    (f : List (Fin 4) → ℚ) (hf0 : ∀ l, 0 ≤ f l) (c : ℚ) (hc : ∀ l ∈ ext steps p, c ≤ f l) :
    (famCount ch (steps - p.length) : ℚ) * c +
      ((ws.filter (fun w => !famTest ch w)).map (fun w => f (p ++ w))).sum
      ≤ ∑ l ∈ ext steps p, f l := by
  obtain ⟨hFm, hcard, htest⟩ := famChoose_spec hn hlen hch
  set Fm := famSet ch p (steps - p.length)
  set rest := ws.filter (fun w => !famTest ch w)
  have hrest : witnessesOK steps p rest :=
    ⟨hw.1.filter _, fun w hw' => hw.2 w (List.mem_of_mem_filter hw')⟩
  set Wr := (rest.map (p ++ ·)).toFinset
  have hWr : Wr ⊆ ext steps p := wit_subset hrest
  have hdisj : Disjoint Fm Wr := by
    rw [Finset.disjoint_left]
    intro l hl hl'
    rw [List.mem_toFinset, List.mem_map] at hl'
    obtain ⟨w, hw', rfl⟩ := hl'
    have h1 := htest w hl
    have h2 := (List.mem_filter.1 hw').2
    simp [h1] at h2
  have e1 : (rest.map (fun w => f (p ++ w))).sum = ∑ l ∈ Wr, f l := wit_sum p rest hrest.1 f
  have e2 : (famCount ch (steps - p.length) : ℚ) * c ≤ ∑ l ∈ Fm, f l := by
    calc (famCount ch (steps - p.length) : ℚ) * c = ∑ _l ∈ Fm, c := by
          rw [Finset.sum_const, hcard, nsmul_eq_mul]
      _ ≤ ∑ l ∈ Fm, f l := Finset.sum_le_sum (fun l hl => hc l (hFm hl))
  rw [e1]
  calc (famCount ch (steps - p.length) : ℚ) * c + ∑ l ∈ Wr, f l
        ≤ ∑ l ∈ Fm, f l + ∑ l ∈ Wr, f l := by linarith
    _ = ∑ l ∈ Fm ∪ Wr, f l := (Finset.sum_union hdisj).symm
    _ ≤ ∑ l ∈ ext steps p, f l :=
        Finset.sum_le_sum_of_subset_of_nonneg (Finset.union_subset hFm hWr) (fun l _ _ => hf0 l)

theorem memEnc_lower_sound (k : ℕ) (ws : List (List (Fin 4))) (steps : ℕ) (S : List Bool)
    (hS : S.length = steps + 2) (q : ℚ) (hq : 1 ≤ q) (thr : ℕ) (p : List (Fin 4))
    (hlen : p.length ≤ steps) (hn : (walkFrom p).Nodup) (hw : witnessesOK steps p ws) :
    (memEnc k ws steps S q thr p).zLo ≤ massZ steps S q p ∧
    (memEnc k ws steps S q thr p).eLo ≤ massA steps S q thr p ∧
    (memEnc k ws steps S q thr p).nLo ≤ massN steps S q thr p ∧
    0 ≤ (memEnc k ws steps S q thr p).zLo ∧ 0 ≤ (memEnc k ws steps S q thr p).eLo ∧
    0 ≤ (memEnc k ws steps S q thr p).nLo := by
  have hq0 : (0 : ℚ) ≤ q := by linarith
  have hpos : ∀ k : ℕ, (0 : ℚ) ≤ q ^ k := fun k => pow_nonneg hq0 k
  set T := ext steps p
  set K : List (Fin 4) → ℕ := fun l => parityContactCount S (walkFrom l) with hK
  let fZ : List (Fin 4) → ℚ := fun l => q ^ K l
  let fA : List (Fin 4) → ℚ := fun l => if thr ≤ K l then q ^ K l else 0
  let fN : List (Fin 4) → ℚ := fun l => if K l < thr then q ^ K l else 0
  have hfZ0 : ∀ l, 0 ≤ fZ l := fun l => hpos _
  have hfA0 : ∀ l, 0 ≤ fA l := fun l => by simp only [fA]; split_ifs <;> simp [hpos]
  have hfN0 : ∀ l, 0 ≤ fN l := fun l => by simp only [fN]; split_ifs <;> simp [hpos]
  have mZ : massZ steps S q p = ∑ l ∈ T, fZ l := massZ_eq steps S hS q p
  have mA : massA steps S q thr p = ∑ l ∈ T, fA l := massA_eq steps S hS q thr p
  have mN : massN steps S q thr p = ∑ l ∈ T, fN l := massN_eq steps S hS q thr p
  have hWT : (ws.map (p ++ ·)).toFinset ⊆ T := wit_subset hw
  have plain : ∀ f : List (Fin 4) → ℚ, (∀ l, 0 ≤ f l) →
      (ws.map (fun w => f (p ++ w))).sum ≤ ∑ l ∈ T, f l ∧
        0 ≤ (ws.map (fun w => f (p ++ w))).sum := by
    intro f hf0
    rw [wit_sum p ws hw.1 f]
    exact ⟨Finset.sum_le_sum_of_subset_of_nonneg hWT (fun l _ _ => hf0 l),
      Finset.sum_nonneg (fun l _ => hf0 l)⟩
  have restNonneg : ∀ (ch : FamChoice) (f : List (Fin 4) → ℚ), (∀ l, 0 ≤ f l) →
      0 ≤ ((ws.filter (fun w => !famTest ch w)).map (fun w => f (p ++ w))).sum :=
    fun ch f hf0 => List.sum_nonneg (by
      intro x hx; obtain ⟨w, _, rfl⟩ := List.mem_map.1 hx; exact hf0 _)
  have hfacts : ∀ l ∈ T, knownContacts S p ≤ K l ∧ K l ≤ kcN steps S p := by
    intro l hl
    obtain ⟨a, b, c, d⟩ := mem_facts steps S hS k (p := p) hl
    exact ⟨a.trans b, c.trans d⟩
  have eZ : (ws.map (witZ S q p)) = ws.map (fun w => fZ (p ++ w)) := rfl
  have eA : (ws.map (witA S q thr p)) = ws.map (fun w => fA (p ++ w)) := rfl
  have eN : (ws.map (witN S q thr p)) = ws.map (fun w => fN (p ++ w)) := rfl
  simp only [memEnc, memLowers]
  rw [mZ, mA, mN]
  cases hpick : famChoose p with
  | none =>
      simp only
      rw [eZ, eA, eN]
      obtain ⟨z1, z0⟩ := plain fZ hfZ0
      obtain ⟨a1, a0⟩ := plain fA hfA0
      obtain ⟨n1, n0⟩ := plain fN hfN0
      exact ⟨z1, a1, n1, z0, a0, n0⟩
  | some ch =>
      have hF0 : (0 : ℚ) ≤ (famCount ch (steps - p.length) : ℚ) * q ^ knownContacts S p := by
        positivity
      simp only
      have eZ' : ((ws.filter (fun w => !famTest ch w)).map (witZ S q p)) =
          (ws.filter (fun w => !famTest ch w)).map (fun w => fZ (p ++ w)) := rfl
      have eA' : ((ws.filter (fun w => !famTest ch w)).map (witA S q thr p)) =
          (ws.filter (fun w => !famTest ch w)).map (fun w => fA (p ++ w)) := rfl
      have eN' : ((ws.filter (fun w => !famTest ch w)).map (witN S q thr p)) =
          (ws.filter (fun w => !famTest ch w)).map (fun w => fN (p ++ w)) := rfl
      rw [eZ', eA', eN']
      simp only [eA, eN]
      refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
      · exact fam_wit_lower hn hlen hw hpick fZ hfZ0 _ (fun l hl => qpow_mono hq (hfacts l hl).1)
      · split_ifs with ht
        · refine fam_wit_lower hn hlen hw hpick fA hfA0 _ (fun l hl => ?_)
          simp only [fA]; rw [if_pos (ht.trans (hfacts l hl).1)]
          exact qpow_mono hq (hfacts l hl).1
        · exact (plain fA hfA0).1
      · split_ifs with ht
        · refine fam_wit_lower hn hlen hw hpick fN hfN0 _ (fun l hl => ?_)
          simp only [fN]; rw [if_pos (lt_of_le_of_lt (hfacts l hl).2 ht)]
          exact qpow_mono hq (hfacts l hl).1
        · exact (plain fN hfN0).1
      · have := restNonneg ch fZ hfZ0; linarith
      · split_ifs
        · have := restNonneg ch fA hfA0; linarith
        · exact (plain fA hfA0).2
      · split_ifs
        · have := restNonneg ch fN hfN0; linarith
        · exact (plain fN hfN0).2

theorem memEnc_sound (k : ℕ) (ws : List (List (Fin 4))) (steps : ℕ) (S : List Bool)
    (hS : S.length = steps + 2) (q : ℚ) (hq : 1 ≤ q) (thr : ℕ) (p : List (Fin 4))
    (hlen : p.length ≤ steps) (hn : (walkFrom p).Nodup) (hw : witnessesOK steps p ws) :
    Encloses steps S q thr p (memEnc k ws steps S q thr p) := by
  obtain ⟨u1, u2, u3⟩ := memEnc_upper_sound k ws steps S hS q hq thr p hlen hn hw
  obtain ⟨l1, l2, l3, l4, l5, l6⟩ := memEnc_lower_sound k ws steps S hS q hq thr p hlen hn hw
  exact ⟨l1, u1, l2, u2, l3, u3, l4, l5, l6⟩

/-! ### Soundness of the version-2 checker -/

theorem checkNode_sound (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2)
    (q : ℚ) (hq : 1 ≤ q) (thr : ℕ) :
    ∀ (t : NTree) (p : List (Fin 4)) (e : Enclosure),
      checkNode steps S q thr p t = some e → Encloses steps S q thr p e := by
  intro t
  induction t with
  | exact =>
      intro p e h
      simp only [checkNode] at h
      split_ifs at h with hc
      cases h
      exact exactEnc_sound steps S hS q hq thr p hc.1 hc.2
  | dead =>
      intro p e h
      simp only [checkNode] at h
      split_ifs at h with hc
      cases h
      exact masses_of_empty (deadOK_ext_empty hc)
  | enum =>
      intro p e h
      simp only [checkNode] at h
      split_ifs at h with hc
      cases h
      exact enumEnc_sound steps S hS q hq thr _ p (by omega)
  | legacy tr ws =>
      intro p e h
      simp only [checkNode] at h
      split_ifs at h with hc
      cases h
      exact famRaise_sound tr steps S hS q hq thr p hc.1 hc.2.1 _
        (boundEnc_sound tr ws steps S hS q hq thr p hc.1 hc.2.1 hc.2.2)
  | mem k ws =>
      intro p e h
      simp only [checkNode] at h
      split_ifs at h with hc
      cases h
      exact memEnc_sound k ws steps S hS q hq thr p hc.1 hc.2.1 hc.2.2
  | node t0 t1 t2 t3 ih0 ih1 ih2 ih3 =>
      intro p e h
      simp only [checkNode] at h
      split_ifs at h with hc
      cases h0 : checkNode steps S q thr (p ++ [0]) t0 <;>
      cases h1 : checkNode steps S q thr (p ++ [1]) t1 <;>
      cases h2 : checkNode steps S q thr (p ++ [2]) t2 <;>
      cases h3 : checkNode steps S q thr (p ++ [3]) t3 <;>
      simp only [h0, h1, h2, h3, reduceCtorEq] at h
      cases h
      exact Encloses_node hc.1 (ih0 _ _ h0) (ih1 _ _ h1) (ih2 _ _ h2) (ih3 _ _ h3)

theorem check_eq_some {S : List Bool} {q : ℚ} {thr : ℕ} {t : NTree} {r : FrontierResult}
    (h : check S q thr t = some r) :
    2 ≤ S.length ∧ 1 ≤ q ∧ ∃ e, checkNode (S.length - 2) S q thr [] t = some e ∧
      r = finalize e := by
  unfold check at h
  split_ifs at h with hc
  obtain ⟨e, he, rfl⟩ := Option.map_eq_some_iff.1 h
  exact ⟨hc.1, hc.2, e, he, rfl⟩

/-- **Flagship theorem (version 2).** If the version-2 checker accepts an
untrusted certificate `t` and returns `r`, then `r` encloses the exact
partition function, event mass, complement mass and equilibrium
probability `P(k ≥ thr)` of the Lean HP model, for every H/P sequence,
every rational `q` and every threshold. -/
theorem check_sound (steps : ℕ) (seq : Fin (steps + 2) → Bool) (q : ℚ) (thr : ℕ)
    (t : NTree) (r : FrontierResult) (h : check (List.ofFn seq) q thr t = some r) :
    Executable.massQ (weightQ steps q seq) univ ∈ Set.Icc r.zLo r.zHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) ∈
      Set.Icc r.eLo r.eHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ ∈
      Set.Icc r.nLo r.nHi ∧
    Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) /
      Executable.massQ (weightQ steps q seq) univ ∈ Set.Icc r.pLo r.pHi := by
  obtain ⟨_, hq, e, he, rfl⟩ := check_eq_some h
  have hlen : (List.ofFn seq).length - 2 = steps := by simp
  rw [hlen] at he
  have hS : (List.ofFn seq).length = steps + 2 := by simp
  obtain ⟨z1, z2, a1, a2, n1, n2, _, a0, n0⟩ :=
    checkNode_sound steps (List.ofFn seq) hS q hq thr t [] e he
  rw [massQ_univ_eq, massQ_event_eq, massQ_compl_eq]
  refine ⟨⟨z1, z2⟩, ⟨a1, a2⟩, ⟨n1, n2⟩, ?_⟩
  rw [massZ_eq_add steps (List.ofFn seq) q thr []]
  have hpos := massZ_root_pos steps (List.ofFn seq) (q := q) (by linarith)
  rw [massZ_eq_add steps (List.ofFn seq) q thr []] at hpos
  exact ratio_mem_coupled_rat a0 n0 ⟨a1, a2⟩ ⟨n1, n2⟩ hpos

/-- The probability enclosure over the reals, in the form of
`LatticeModel.checkHP_coupled_sound`. -/
theorem check_coupled_probability_sound (steps : ℕ) (seq : Fin (steps + 2) → Bool)
    (q : ℚ) (thr : ℕ) (t : NTree) (r : FrontierResult)
    (h : check (List.ofFn seq) q thr t = some r) :
    ((Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) : ℝ) /
      ((Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) : ℝ) +
       (Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ : ℝ))) ∈
      Set.Icc (r.pLo : ℝ) (r.pHi : ℝ) := by
  have hsound := check_sound steps seq q thr t r h
  have hsplit : Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr) +
      Executable.massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ =
      Executable.massQ (weightQ steps q seq) univ :=
    Finset.sum_add_sum_compl _ _
  have hp := hsound.2.2.2
  rw [← hsplit] at hp
  exact ⟨by exact_mod_cast hp.1, by exact_mod_cast hp.2⟩

end FrontierNext
end FoldBounds
