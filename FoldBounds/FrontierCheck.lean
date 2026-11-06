module

public import FoldBounds.Transfer
public import FoldBounds.FrontierDegree
public import FoldBounds.FrontierFamily

/-!
# A compact, non-enumerative frontier certificate checker

An untrusted producer submits a *prefix trie* (`FTree`) whose leaves are

* `exact` — a complete valid conformation, scored exactly;
* `bound tr ws` — a valid prefix standing for its whole (possibly huge)
  completion set, bounded by a counting × contact envelope (`tr = false`)
  or additionally by the verified transfer relaxation (`tr = true`), and
  optionally carrying explicit witness completions `ws` (suffix words)
  that supply exact lower mass and tighten the upper bound;
* `dead` — a prefix proved locally to have no valid completion.

Internal `node t0 t1 t2 t3` entries branch on the four directions.  The
trie shape makes duplicate frontier entries, overlapping ancestor /
descendant entries, hidden descendants and ambiguous encodings
unrepresentable; everything else (geometry, self-avoidance, terminality,
deadness, witness validity and distinctness) is re-checked by `checkNode`.

`checkNode` visits only the submitted trie (plus a bounded local amount
of work per node); it never enumerates the conformation space.

The main theorem `check_sound` states that every accepted certificate
yields rational enclosures of the exact partition function, event mass,
complement mass and equilibrium probability of the Lean HP model.
-/
@[expose] public section

namespace FoldBounds
namespace FrontierCheck

open LatticeModel Frontier Transfer Finset

/-! ### Certificate format -/

/-- An untrusted prefix-trie frontier certificate. -/
inductive FTree where
  | exact : FTree
  | bound (transfer : Bool) (witnesses : List (List (Fin 4))) : FTree
  | dead : FTree
  | node (t0 t1 t2 t3 : FTree) : FTree
  deriving Repr, DecidableEq

/-- Six rational enclosure endpoints for one subtree. -/
structure Enclosure where
  zLo : ℚ
  zHi : ℚ
  eLo : ℚ
  eHi : ℚ
  nLo : ℚ
  nHi : ℚ
  deriving Repr, DecidableEq

instance : Add Enclosure :=
  ⟨fun a b => ⟨a.zLo + b.zLo, a.zHi + b.zHi, a.eLo + b.eLo, a.eHi + b.eHi,
    a.nLo + b.nLo, a.nHi + b.nHi⟩⟩

instance : Zero Enclosure := ⟨⟨0, 0, 0, 0, 0, 0⟩⟩

/-! ### Exact leaves -/

/-- Exact contribution of a complete conformation. -/
def exactEnc (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) : Enclosure :=
  let k := parityContactCount S (walkFrom p)
  let w := q ^ k
  if thr ≤ k then ⟨w, w, w, w, 0, 0⟩ else ⟨w, w, 0, 0, w, w⟩

/-! ### Bounded leaves -/

/-- Basic contact cap `K_max(p)`. -/
def kmax (S : List Bool) (p : List (Fin 4)) : ℕ :=
  knownContacts S p + futurePairs S (p.length + 2)

/-- Contact cap used by a bound leaf: Python's `K_max` in basic mode, the
degree-capped `known + futureCap` in transfer mode. -/
def kcap (tr : Bool) (S : List Bool) (p : List (Fin 4)) : ℕ :=
  if tr then knownContacts S p + futureCap S (p.length + 2) else kmax S p

/-- Upper contact exponent as a function of the realised cross score. -/
def envExp (tr : Bool) (S : List Bool) (p : List (Fin 4)) (k : ℕ) : ℕ :=
  if tr then min (knownContacts S p + k + ffCap S (p.length + 2)) (kcap tr S p) else kcap tr S p

/-- Lower contact exponent as a function of the realised cross score. -/
def lowExp (tr : Bool) (S : List Bool) (p : List (Fin 4)) (k : ℕ) : ℕ :=
  if tr then knownContacts S p + k else knownContacts S p

/-- Envelope for the total weight. -/
def gZ (tr : Bool) (S : List Bool) (q : ℚ) (p : List (Fin 4)) (k : ℕ) : ℚ :=
  q ^ envExp tr S p k

/-- Envelope for the event weight. -/
def gA (tr : Bool) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) (k : ℕ) : ℚ :=
  if thr ≤ envExp tr S p k then q ^ envExp tr S p k else 0

/-- Envelope for the complement weight. -/
def gN (tr : Bool) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) (k : ℕ) : ℚ :=
  if lowExp tr S p k < thr then q ^ envExp tr S p k else 0

/-- Basic counting × contact upper bound `3^(n-m) q^{K_max}`. -/
def basicZ (tr : Bool) (steps : ℕ) (S : List Bool) (q : ℚ) (p : List (Fin 4)) : ℚ :=
  3 ^ (steps - p.length) * q ^ kcap tr S p

def basicA (tr : Bool) (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) : ℚ :=
  if thr ≤ kcap tr S p then basicZ tr steps S q p else 0

def basicN (tr : Bool) (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) : ℚ :=
  if knownContacts S p < thr then basicZ tr steps S q p else 0

/-- Upper bounds on the three envelope sums (total, event, complement). -/
def envelopeUppers (tr : Bool) (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ)
    (p : List (Fin 4)) : ℚ × ℚ × ℚ :=
  if tr then
    let tbl := relaxedTable steps S p
    (min (basicZ tr steps S q p) (wsum (fun st => gZ tr S q p st.2.2.1) tbl),
     min (basicA tr steps S q thr p) (wsum (fun st => gA tr S q thr p st.2.2.1) tbl),
     min (basicN tr steps S q thr p) (wsum (fun st => gN tr S q thr p st.2.2.1) tbl))
  else (basicZ tr steps S q p, basicA tr steps S q thr p, basicN tr steps S q thr p)

/-- Exact contact count of a witness completion `p ++ w`. -/
def witK (S : List Bool) (p w : List (Fin 4)) : ℕ := parityContactCount S (walkFrom (p ++ w))

/-- Cross score of a witness completion. -/
def witC (S : List Bool) (p w : List (Fin 4)) : ℕ :=
  crossContacts S (p.length + 2) (walkFrom (p ++ w))

def witZ (S : List Bool) (q : ℚ) (p w : List (Fin 4)) : ℚ := q ^ witK S p w

def witA (S : List Bool) (q : ℚ) (thr : ℕ) (p w : List (Fin 4)) : ℚ :=
  if thr ≤ witK S p w then q ^ witK S p w else 0

def witN (S : List Bool) (q : ℚ) (thr : ℕ) (p w : List (Fin 4)) : ℚ :=
  if witK S p w < thr then q ^ witK S p w else 0

/-- Enclosure of a bounded subtree, with optional verified witnesses. -/
def boundEnc (tr : Bool) (ws : List (List (Fin 4))) (steps : ℕ) (S : List Bool) (q : ℚ)
    (thr : ℕ) (p : List (Fin 4)) : Enclosure :=
  let U := envelopeUppers tr steps S q thr p
  ⟨(ws.map (witZ S q p)).sum,
   U.1 - (ws.map (fun w => gZ tr S q p (witC S p w) - witZ S q p w)).sum,
   (ws.map (witA S q thr p)).sum,
   U.2.1 - (ws.map (fun w => gA tr S q thr p (witC S p w) - witA S q thr p w)).sum,
   (ws.map (witN S q thr p)).sum,
   U.2.2 - (ws.map (fun w => gN tr S q thr p (witC S p w) - witN S q thr p w)).sum⟩

/-- Certified-family lower bound `2^r q^{known}` (transfer mode only; it is
`0` unless a diagonal family is locally available at `p`). -/
def famLow (tr : Bool) (steps : ℕ) (S : List Bool) (q : ℚ) (p : List (Fin 4)) : ℚ :=
  if tr && famOK p then 2 ^ (steps - p.length) * q ^ knownContacts S p else 0

/-- Raise the lower endpoints of a bound-leaf enclosure by the family lower
bound.  Each new endpoint is the maximum of two independently sound lower
bounds (witnesses and family may overlap, so they are never added). -/
def famRaise (tr : Bool) (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4))
    (e : Enclosure) : Enclosure :=
  let L := famLow tr steps S q p
  ⟨max e.zLo L, e.zHi,
   if thr ≤ knownContacts S p then max e.eLo L else e.eLo, e.eHi,
   if kcap tr S p < thr then max e.nLo L else e.nLo, e.nHi⟩

/-- Witness validity: distinct suffixes, each giving a complete valid
conformation extending `p`. -/
def witnessesOK (steps : ℕ) (p : List (Fin 4)) (ws : List (List (Fin 4))) : Prop :=
  ws.Nodup ∧ ∀ w ∈ ws, (p ++ w).length = steps ∧ (walkFrom (p ++ w)).Nodup

instance (steps : ℕ) (p : List (Fin 4)) (ws : List (List (Fin 4))) :
    Decidable (witnessesOK steps p ws) := by
  unfold witnessesOK; infer_instance

/-! ### Dead leaves -/

/-- A locally checkable proof of emptiness: the prefix itself collides, or
it is nonterminal and all four one-step extensions collide. -/
def deadOK (steps : ℕ) (p : List (Fin 4)) : Prop :=
  ¬ (walkFrom p).Nodup ∨ (p.length < steps ∧ ∀ d : Fin 4, ¬ (walkFrom (p ++ [d])).Nodup)

instance (steps : ℕ) (p : List (Fin 4)) : Decidable (deadOK steps p) := by
  unfold deadOK; infer_instance

/-! ### The checker -/

/-- Recursive trie checker. It visits each submitted trie node once. -/
def checkNode (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) :
    List (Fin 4) → FTree → Option Enclosure
  | p, .exact =>
      if p.length = steps ∧ (walkFrom p).Nodup then some (exactEnc S q thr p) else none
  | p, .dead => if deadOK steps p then some 0 else none
  | p, .bound tr ws =>
      if p.length ≤ steps ∧ (walkFrom p).Nodup ∧ witnessesOK steps p ws then
        some (famRaise tr steps S q thr p (boundEnc tr ws steps S q thr p))
      else none
  | p, .node t0 t1 t2 t3 =>
      if p.length < steps ∧ (walkFrom p).Nodup then
        match checkNode steps S q thr (p ++ [0]) t0, checkNode steps S q thr (p ++ [1]) t1,
            checkNode steps S q thr (p ++ [2]) t2, checkNode steps S q thr (p ++ [3]) t3 with
        | some a, some b, some c, some d => some (a + b + c + d)
        | _, _, _, _ => none
      else none

/-- Final rational result returned by the checker. -/
structure FrontierResult where
  zLo : ℚ
  zHi : ℚ
  eLo : ℚ
  eHi : ℚ
  nLo : ℚ
  nHi : ℚ
  pLo : ℚ
  pHi : ℚ
  deriving Repr, DecidableEq

/-- Coupled probability endpoints (Lean convention `x / 0 = 0`). -/
def finalize (e : Enclosure) : FrontierResult :=
  ⟨e.zLo, e.zHi, e.eLo, e.eHi, e.nLo, e.nHi,
   e.eLo / (e.eLo + e.nHi), e.eHi / (e.eHi + e.nLo)⟩

/-- **The frontier checker.** `S` is the H/P sequence (`true` = H), `q ≥ 1`
the exact rational Boltzmann factor, `thr` the event threshold
(event: at least `thr` H–H contacts), `t` the untrusted certificate. -/
def check (S : List Bool) (q : ℚ) (thr : ℕ) (t : FTree) : Option FrontierResult :=
  if 2 ≤ S.length ∧ 1 ≤ q then (checkNode (S.length - 2) S q thr [] t).map finalize else none

/-! ### Semantics -/

/-- Exact subtree partition function. -/
def massZ (steps : ℕ) (S : List Bool) (q : ℚ) (p : List (Fin 4)) : ℚ :=
  ∑ l ∈ ext steps p, q ^ wordContacts S l

/-- Exact subtree event mass (`≥ thr` contacts). -/
def massA (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) : ℚ :=
  ∑ l ∈ ext steps p, if thr ≤ wordContacts S l then q ^ wordContacts S l else 0

/-- Exact subtree complement mass (`< thr` contacts). -/
def massN (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4)) : ℚ :=
  ∑ l ∈ ext steps p, if wordContacts S l < thr then q ^ wordContacts S l else 0

/-- Soundness predicate for one subtree enclosure. -/
def Encloses (steps : ℕ) (S : List Bool) (q : ℚ) (thr : ℕ) (p : List (Fin 4))
    (e : Enclosure) : Prop :=
  e.zLo ≤ massZ steps S q p ∧ massZ steps S q p ≤ e.zHi ∧
  e.eLo ≤ massA steps S q thr p ∧ massA steps S q thr p ≤ e.eHi ∧
  e.nLo ≤ massN steps S q thr p ∧ massN steps S q thr p ≤ e.nHi ∧
  0 ≤ e.zLo ∧ 0 ≤ e.eLo ∧ 0 ≤ e.nLo

/-! ### Soundness proofs -/

theorem mem_ext_valid {steps : ℕ} {p l : List (Fin 4)} (hl : l ∈ ext steps p) :
    l ∈ validWords steps := by
  rw [mem_ext] at hl
  exact mem_validWords.2 ⟨hl.1, hl.2.1⟩

/-- Partition function written with the Python-style scorer. -/
theorem massZ_eq (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (q : ℚ)
    (p : List (Fin 4)) :
    massZ steps S q p = ∑ l ∈ ext steps p, q ^ parityContactCount S (walkFrom l) :=
  Finset.sum_congr rfl (fun l hl => by rw [wordContacts_eq_parity S hS (mem_ext_valid hl)])

theorem massA_eq (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (q : ℚ) (thr : ℕ)
    (p : List (Fin 4)) :
    massA steps S q thr p = ∑ l ∈ ext steps p,
      if thr ≤ parityContactCount S (walkFrom l) then q ^ parityContactCount S (walkFrom l)
      else 0 :=
  Finset.sum_congr rfl (fun l hl => by rw [wordContacts_eq_parity S hS (mem_ext_valid hl)])

theorem massN_eq (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (q : ℚ) (thr : ℕ)
    (p : List (Fin 4)) :
    massN steps S q thr p = ∑ l ∈ ext steps p,
      if parityContactCount S (walkFrom l) < thr then q ^ parityContactCount S (walkFrom l)
      else 0 :=
  Finset.sum_congr rfl (fun l hl => by rw [wordContacts_eq_parity S hS (mem_ext_valid hl)])

/-- **Envelope bound with exact witnesses.** If `f ≤ G` on `T`, the envelope
sum is at most `U`, and `W ⊆ T`, then the true sum is at most
`U - Σ_W (G - f)`. -/
theorem envelope_bound {ι : Type*} {T W : Finset ι} (hW : W ⊆ T) (f G : ι → ℚ)
    (hfG : ∀ x ∈ T, f x ≤ G x) {U : ℚ} (hU : ∑ x ∈ T, G x ≤ U) :
    ∑ x ∈ T, f x ≤ U - ∑ x ∈ W, (G x - f x) := by
  have h1 : ∑ x ∈ W, (G x - f x) ≤ ∑ x ∈ T, (G x - f x) :=
    Finset.sum_le_sum_of_subset_of_nonneg hW (fun x hx _ => by linarith [hfG x hx])
  have h2 : ∑ x ∈ T, (G x - f x) = ∑ x ∈ T, G x - ∑ x ∈ T, f x := Finset.sum_sub_distrib _ _
  linarith

/-- Witness lists sum like the corresponding finite sets of completions. -/
theorem wit_sum (p : List (Fin 4)) (ws : List (List (Fin 4))) (hnd : ws.Nodup)
    (H : List (Fin 4) → ℚ) :
    (ws.map (fun w => H (p ++ w))).sum = ∑ l ∈ (ws.map (p ++ ·)).toFinset, H l := by
  rw [List.sum_toFinset H (hnd.map (List.append_right_injective p)), List.map_map]
  rfl

theorem wit_subset {steps : ℕ} {p : List (Fin 4)} {ws : List (List (Fin 4))}
    (hw : witnessesOK steps p ws) : (ws.map (p ++ ·)).toFinset ⊆ ext steps p := by
  intro l hl
  rw [List.mem_toFinset, List.mem_map] at hl
  obtain ⟨w, hw', rfl⟩ := hl
  rw [mem_ext]
  exact ⟨(hw.2 w hw').1, (hw.2 w hw').2, List.prefix_append p w⟩

theorem envExp_le_kcap (tr : Bool) (S : List Bool) (p : List (Fin 4)) (k : ℕ) :
    envExp tr S p k ≤ kcap tr S p := by
  unfold envExp; split_ifs <;> omega

/-- Contact facts for genuine completions: `lowExp ≤ k ≤ envExp`. -/
theorem exp_facts (tr : Bool) (S : List Bool) {p l : List (Fin 4)} (h : p <+: l)
    (hlen : S.length = l.length + 2) (hn : (walkFrom l).Nodup) :
    lowExp tr S p (crossContacts S (p.length + 2) (walkFrom l)) ≤
      parityContactCount S (walkFrom l) ∧
    parityContactCount S (walkFrom l) ≤
      envExp tr S p (crossContacts S (p.length + 2) (walkFrom l)) := by
  have h1 := contacts_le_cross_bound S h
  have h2 := contacts_le_kmax S h
  have h3 := knownContacts_le S h
  have h4 := contacts_le_degree_caps S h hlen hn
  unfold lowExp envExp kcap kmax
  cases tr <;> simp only [Bool.false_eq_true, if_false, if_true] <;> omega

theorem sum_le_pow_three_mul {steps : ℕ} {p : List (Fin 4)} (G : List (Fin 4) → ℚ) (C : ℚ)
    (hC : 0 ≤ C) (hG : ∀ l ∈ ext steps p, G l ≤ C) :
    ∑ l ∈ ext steps p, G l ≤ 3 ^ (steps - p.length) * C := by
  calc ∑ l ∈ ext steps p, G l ≤ (ext steps p).card • C := Finset.sum_le_card_nsmul _ _ _ hG
    _ = ((ext steps p).card : ℚ) * C := by rw [nsmul_eq_mul]
    _ ≤ 3 ^ (steps - p.length) * C := by
        apply mul_le_mul_of_nonneg_right _ hC
        exact_mod_cast card_ext_le_pow_three steps p

theorem boundEnc_sound (tr : Bool) (ws : List (List (Fin 4))) (steps : ℕ) (S : List Bool)
    (hS : S.length = steps + 2) (q : ℚ) (hq : 1 ≤ q) (thr : ℕ) (p : List (Fin 4))
    (hlen : p.length ≤ steps) (hn : (walkFrom p).Nodup) (hw : witnessesOK steps p ws) :
    Encloses steps S q thr p (boundEnc tr ws steps S q thr p) := by
  have hq0 : (0 : ℚ) ≤ q := by linarith
  set T := ext steps p with hT
  set W := (ws.map (p ++ ·)).toFinset with hWdef
  have hWT : W ⊆ T := wit_subset hw
  set c : List (Fin 4) → ℕ := fun l => crossContacts S (p.length + 2) (walkFrom l) with hc
  set K : List (Fin 4) → ℕ := fun l => parityContactCount S (walkFrom l) with hK
  have hfacts : ∀ l ∈ T, lowExp tr S p (c l) ≤ K l ∧ K l ≤ envExp tr S p (c l) := by
    intro l hl
    have hl' := mem_ext.1 hl
    exact exp_facts tr S hl'.2.2 (by omega) hl'.2.1
  -- the three true weights and envelopes
  let fZ : List (Fin 4) → ℚ := fun l => q ^ K l
  let fA : List (Fin 4) → ℚ := fun l => if thr ≤ K l then q ^ K l else 0
  let fN : List (Fin 4) → ℚ := fun l => if K l < thr then q ^ K l else 0
  let GZ : List (Fin 4) → ℚ := fun l => gZ tr S q p (c l)
  let GA : List (Fin 4) → ℚ := fun l => gA tr S q thr p (c l)
  let GN : List (Fin 4) → ℚ := fun l => gN tr S q thr p (c l)
  have hpos : ∀ k : ℕ, (0 : ℚ) ≤ q ^ k := fun k => pow_nonneg hq0 k
  have hfZ : ∀ l ∈ T, fZ l ≤ GZ l := fun l hl => qpow_mono hq (hfacts l hl).2
  have hfA : ∀ l ∈ T, fA l ≤ GA l := by
    intro l hl
    simp only [fA, GA, gA]
    split_ifs with h1 h2
    · exact qpow_mono hq (hfacts l hl).2
    · exact absurd (h1.trans (hfacts l hl).2) h2
    · exact hpos _
    · exact le_refl _
  have hfN : ∀ l ∈ T, fN l ≤ GN l := by
    intro l hl
    simp only [fN, GN, gN]
    split_ifs with h1 h2
    · exact qpow_mono hq (hfacts l hl).2
    · exact absurd (lt_of_le_of_lt (hfacts l hl).1 h1) h2
    · exact hpos _
    · exact le_refl _
  have hfZ0 : ∀ l, 0 ≤ fZ l := fun l => hpos _
  have hfA0 : ∀ l, 0 ≤ fA l := fun l => by simp only [fA]; split_ifs <;> simp [hpos]
  have hfN0 : ∀ l, 0 ≤ fN l := fun l => by simp only [fN]; split_ifs <;> simp [hpos]
  have hgZ0 : ∀ k, 0 ≤ gZ tr S q p k := fun k => hpos _
  have hgA0 : ∀ k, 0 ≤ gA tr S q thr p k := fun k => by unfold gA; split_ifs <;> simp [hpos]
  have hgN0 : ∀ k, 0 ≤ gN tr S q thr p k := fun k => by unfold gN; split_ifs <;> simp [hpos]
  -- pointwise basic caps
  have hcapZ : ∀ k, gZ tr S q p k ≤ q ^ kcap tr S p := fun k => qpow_mono hq (envExp_le_kcap _ _ _ _)
  have hcapA : ∀ k, gA tr S q thr p k ≤ if thr ≤ kcap tr S p then q ^ kcap tr S p else 0 := by
    intro k; unfold gA
    split_ifs with h1 h2
    · exact qpow_mono hq (envExp_le_kcap _ _ _ _)
    · exact absurd (h1.trans (envExp_le_kcap _ _ _ _)) h2
    · exact hpos _
    · exact le_refl _
  have hcapN : ∀ k, gN tr S q thr p k ≤
      if knownContacts S p < thr then q ^ kcap tr S p else 0 := by
    intro k; unfold gN
    split_ifs with h1 h2
    · exact qpow_mono hq (envExp_le_kcap _ _ _ _)
    · exfalso; apply h2
      have : knownContacts S p ≤ lowExp tr S p k := by unfold lowExp; split_ifs <;> omega
      omega
    · exact hpos _
    · exact le_refl _
  -- envelope sums are bounded by the computed uppers
  have hbZ : ∑ l ∈ T, GZ l ≤ basicZ tr steps S q p :=
    sum_le_pow_three_mul GZ _ (hpos _) (fun l _ => hcapZ _)
  have hbA : ∑ l ∈ T, GA l ≤ basicA tr steps S q thr p := by
    have := sum_le_pow_three_mul (steps := steps) (p := p) GA
      (if thr ≤ kcap tr S p then q ^ kcap tr S p else 0) (by split_ifs <;> simp [hpos])
      (fun l _ => hcapA _)
    unfold basicA basicZ; split_ifs at this ⊢ <;> simpa using this
  have hbN : ∑ l ∈ T, GN l ≤ basicN tr steps S q thr p := by
    have := sum_le_pow_three_mul (steps := steps) (p := p) GN
      (if knownContacts S p < thr then q ^ kcap tr S p else 0) (by split_ifs <;> simp [hpos])
      (fun l _ => hcapN _)
    unfold basicN basicZ; split_ifs at this ⊢ <;> simpa using this
  have hU : ∑ l ∈ T, GZ l ≤ (envelopeUppers tr steps S q thr p).1 ∧
      ∑ l ∈ T, GA l ≤ (envelopeUppers tr steps S q thr p).2.1 ∧
      ∑ l ∈ T, GN l ≤ (envelopeUppers tr steps S q thr p).2.2 := by
    unfold envelopeUppers
    cases tr
    · exact ⟨hbZ, hbA, hbN⟩
    · simp only [if_true]
      refine ⟨le_min hbZ ?_, le_min hbA ?_, le_min hbN ?_⟩
      · exact sum_ext_le_dp steps S hS p hn hlen _ hgZ0
      · exact sum_ext_le_dp steps S hS p hn hlen _ hgA0
      · exact sum_ext_le_dp steps S hS p hn hlen _ hgN0
  -- witness sums
  have hnd := hw.1
  have eZ : (ws.map (witZ S q p)).sum = ∑ l ∈ W, fZ l := wit_sum p ws hnd fZ
  have eA : (ws.map (witA S q thr p)).sum = ∑ l ∈ W, fA l := wit_sum p ws hnd fA
  have eN : (ws.map (witN S q thr p)).sum = ∑ l ∈ W, fN l := wit_sum p ws hnd fN
  have dZ : (ws.map (fun w => gZ tr S q p (witC S p w) - witZ S q p w)).sum =
      ∑ l ∈ W, (GZ l - fZ l) := wit_sum p ws hnd (fun l => GZ l - fZ l)
  have dA : (ws.map (fun w => gA tr S q thr p (witC S p w) - witA S q thr p w)).sum =
      ∑ l ∈ W, (GA l - fA l) := wit_sum p ws hnd (fun l => GA l - fA l)
  have dN : (ws.map (fun w => gN tr S q thr p (witC S p w) - witN S q thr p w)).sum =
      ∑ l ∈ W, (GN l - fN l) := wit_sum p ws hnd (fun l => GN l - fN l)
  have mZ : massZ steps S q p = ∑ l ∈ T, fZ l := massZ_eq steps S hS q p
  have mA : massA steps S q thr p = ∑ l ∈ T, fA l := massA_eq steps S hS q thr p
  have mN : massN steps S q thr p = ∑ l ∈ T, fN l := massN_eq steps S hS q thr p
  simp only [Encloses, boundEnc]
  rw [eZ, eA, eN, dZ, dA, dN, mZ, mA, mN]
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact Finset.sum_le_sum_of_subset_of_nonneg hWT (fun l _ _ => hfZ0 l)
  · exact envelope_bound hWT fZ GZ hfZ hU.1
  · exact Finset.sum_le_sum_of_subset_of_nonneg hWT (fun l _ _ => hfA0 l)
  · exact envelope_bound hWT fA GA hfA hU.2.1
  · exact Finset.sum_le_sum_of_subset_of_nonneg hWT (fun l _ _ => hfN0 l)
  · exact envelope_bound hWT fN GN hfN hU.2.2
  · exact Finset.sum_nonneg (fun l _ => hfZ0 l)
  · exact Finset.sum_nonneg (fun l _ => hfA0 l)
  · exact Finset.sum_nonneg (fun l _ => hfN0 l)

/-- **Soundness of the family lower bound** for every completion set. -/
theorem famLow_le (tr : Bool) (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (q : ℚ)
    (hq : 1 ≤ q) (thr : ℕ) (p : List (Fin 4)) (hlen : p.length ≤ steps)
    (hn : (walkFrom p).Nodup) :
    famLow tr steps S q p ≤ massZ steps S q p ∧
    (thr ≤ knownContacts S p → famLow tr steps S q p ≤ massA steps S q thr p) ∧
    (kcap tr S p < thr → famLow tr steps S q p ≤ massN steps S q thr p) := by
  have hq0 : (0 : ℚ) ≤ q := by linarith
  unfold famLow
  by_cases hf : (tr && famOK p) = true
  · rw [if_pos hf]
    have hfam : famOK p = true := by simp_all
    have facts : ∀ l ∈ ext steps p, knownContacts S p ≤ parityContactCount S (walkFrom l) ∧
        parityContactCount S (walkFrom l) ≤ kcap tr S p := by
      intro l hl
      have hl' := mem_ext.1 hl
      have hx := exp_facts tr S hl'.2.2 (by omega) hl'.2.1
      refine ⟨le_trans ?_ hx.1, hx.2.trans (envExp_le_kcap _ _ _ _)⟩
      unfold lowExp; split_ifs <;> omega
    refine ⟨?_, ?_, ?_⟩
    · rw [massZ_eq steps S hS]
      exact family_lower hn hlen hfam _ (fun l => pow_nonneg hq0 _) _
        (fun l hl => qpow_mono hq (facts l hl).1)
    · intro ht
      rw [massA_eq steps S hS]
      refine family_lower hn hlen hfam _ (fun l => by split_ifs <;> positivity) _
        (fun l hl => ?_)
      rw [if_pos (ht.trans (facts l hl).1)]
      exact qpow_mono hq (facts l hl).1
    · intro ht
      rw [massN_eq steps S hS]
      refine family_lower hn hlen hfam _ (fun l => by split_ifs <;> positivity) _
        (fun l hl => ?_)
      rw [if_pos (lt_of_le_of_lt (facts l hl).2 ht)]
      exact qpow_mono hq (facts l hl).1
  · rw [if_neg hf]
    refine ⟨?_, fun _ => ?_, fun _ => ?_⟩
    · exact Finset.sum_nonneg (fun l _ => pow_nonneg hq0 _)
    · exact Finset.sum_nonneg (fun l _ => by split_ifs <;> positivity)
    · exact Finset.sum_nonneg (fun l _ => by split_ifs <;> positivity)

theorem famRaise_sound (tr : Bool) (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2)
    (q : ℚ) (hq : 1 ≤ q) (thr : ℕ) (p : List (Fin 4)) (hlen : p.length ≤ steps)
    (hn : (walkFrom p).Nodup) (e : Enclosure) (he : Encloses steps S q thr p e) :
    Encloses steps S q thr p (famRaise tr steps S q thr p e) := by
  obtain ⟨z1, z2, a1, a2, n1, n2, z0, a0, n0⟩ := he
  obtain ⟨fz, fa, fn⟩ := famLow_le tr steps S hS q hq thr p hlen hn
  unfold famRaise
  refine ⟨max_le z1 fz, z2, ?_, a2, ?_, n2, le_max_of_le_left z0, ?_, ?_⟩
  · dsimp only; split_ifs with h
    · exact max_le a1 (fa h)
    · exact a1
  · dsimp only; split_ifs with h
    · exact max_le n1 (fn h)
    · exact n1
  · dsimp only; split_ifs
    · exact le_max_of_le_left a0
    · exact a0
  · dsimp only; split_ifs
    · exact le_max_of_le_left n0
    · exact n0

theorem exactEnc_sound (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2) (q : ℚ)
    (hq : 1 ≤ q) (thr : ℕ) (p : List (Fin 4)) (hlen : p.length = steps)
    (hn : (walkFrom p).Nodup) : Encloses steps S q thr p (exactEnc S q thr p) := by
  have hq0 : (0 : ℚ) ≤ q := by linarith
  have hext := ext_eq_singleton hlen hn
  rw [Encloses, massZ_eq steps S hS, massA_eq steps S hS, massN_eq steps S hS, hext]
  simp only [Finset.sum_singleton, exactEnc]
  have := pow_nonneg hq0 (parityContactCount S (walkFrom p))
  split_ifs with h1 h2 <;> first | omega | (refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
    simp_all)

theorem masses_of_empty {steps : ℕ} {S : List Bool} {q : ℚ} {thr : ℕ} {p : List (Fin 4)}
    (h : ext steps p = ∅) : Encloses steps S q thr p 0 := by
  simp [Encloses, massZ, massA, massN, h, show (0 : Enclosure) = ⟨0, 0, 0, 0, 0, 0⟩ from rfl]

theorem deadOK_ext_empty {steps : ℕ} {p : List (Fin 4)} (h : deadOK steps p) :
    ext steps p = ∅ := by
  rcases h with h | ⟨hlt, hall⟩
  · exact ext_eq_empty_of_not_nodup h
  · rw [ext_eq_biUnion hlt]
    ext l
    simp [ext_eq_empty_of_not_nodup (hall _)]

theorem Encloses_node {steps : ℕ} {S : List Bool} {q : ℚ} {thr : ℕ} {p : List (Fin 4)}
    (hlt : p.length < steps) {a b c d : Enclosure}
    (ha : Encloses steps S q thr (p ++ [0]) a) (hb : Encloses steps S q thr (p ++ [1]) b)
    (hc : Encloses steps S q thr (p ++ [2]) c) (hd : Encloses steps S q thr (p ++ [3]) d) :
    Encloses steps S q thr p (a + b + c + d) := by
  have sZ : massZ steps S q p = massZ steps S q (p ++ [0]) + massZ steps S q (p ++ [1]) +
      massZ steps S q (p ++ [2]) + massZ steps S q (p ++ [3]) := by
    rw [massZ, sum_ext_split hlt, Fin.sum_univ_four]; rfl
  have sA : massA steps S q thr p = massA steps S q thr (p ++ [0]) +
      massA steps S q thr (p ++ [1]) + massA steps S q thr (p ++ [2]) +
      massA steps S q thr (p ++ [3]) := by
    rw [massA, sum_ext_split hlt, Fin.sum_univ_four]; rfl
  have sN : massN steps S q thr p = massN steps S q thr (p ++ [0]) +
      massN steps S q thr (p ++ [1]) + massN steps S q thr (p ++ [2]) +
      massN steps S q thr (p ++ [3]) := by
    rw [massN, sum_ext_split hlt, Fin.sum_univ_four]; rfl
  simp only [Encloses] at ha hb hc hd ⊢
  rw [sZ, sA, sN]
  change _ ≤ _ ∧ _ ≤ (a.zHi + b.zHi + c.zHi + d.zHi) ∧ _
  simp only [show ∀ x y : Enclosure, x + y = ⟨x.zLo + y.zLo, x.zHi + y.zHi, x.eLo + y.eLo,
    x.eHi + y.eHi, x.nLo + y.nLo, x.nHi + y.nHi⟩ from fun _ _ => rfl]
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> linarith

/-- **Soundness of the recursive trie checker.** Every enclosure returned
for prefix `p` encloses the exact subtree partition function, event mass
and complement mass of `p`. -/
theorem checkNode_sound (steps : ℕ) (S : List Bool) (hS : S.length = steps + 2)
    (q : ℚ) (hq : 1 ≤ q) (thr : ℕ) :
    ∀ (t : FTree) (p : List (Fin 4)) (e : Enclosure),
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
  | bound tr ws =>
      intro p e h
      simp only [checkNode] at h
      split_ifs at h with hc
      cases h
      exact famRaise_sound tr steps S hS q hq thr p hc.1 hc.2.1 _
        (boundEnc_sound tr ws steps S hS q hq thr p hc.1 hc.2.1 hc.2.2)
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

end FrontierCheck
end FoldBounds
