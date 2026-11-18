module

/-!
# FoldBounds.FrontierKernel — a dependency-light executable frontier checker

This module contains **only executable definitions** and imports nothing
beyond the Lean core library (no Mathlib).  It is the code that the
`foldkernel` binary runs.  Its soundness against the Lean HP lattice model
(`FoldBounds.LatticeModel`, `FoldBounds.Frontier`) is proved in
`FoldBounds/FrontierKernelSound.lean`, which imports this module together with
Mathlib and the existing semantic library.

Design (format `FOLDBOUNDS-KERNEL 1`, see
`docs/research/CERTIFIED_FRONTIER_KERNEL_SPEC.md`):

* The certificate is a byte string: the preorder encoding of a prefix tree
  over the bytes `N E X B -`.  It carries no numbers.  The checker walks it
  once, maintaining the walk built so far *incrementally* (occupied sites in
  reverse order, the end site, the number of placed residues and the number
  of H–H contacts among them).
* `E` adds one conformation to an exact contact histogram; `X` asks the
  checker to enumerate a whole subtree into that histogram; `B` records one
  *bound entry* `(a, K, cnt, fam)` recomputed from the walk; `N` splits on the
  next direction and requires a `-` marker exactly on occupied directions.
* Nothing in the tree walk depends on the Boltzmann factor `q` or on the
  event threshold: the accumulated histogram and bound entries are integers,
  and `evalQuery` turns them into exact rational enclosures for any `(q, t)`.
* `verify` is the function from the raw input string to the final report.
-/

@[expose] public section

namespace FoldBounds
namespace Kernel

/-! ### Lattice primitives (definitionally the ones of `LatticeModel`) -/

/-- A square-lattice site. -/
abbrev Site := Int × Int

/-- The anchor site `(0,0)`. -/
def site0 : Site := (0, 0)
/-- The second site `(1,0)` of every anchored walk. -/
def site1 : Site := (1, 0)

/-- One lattice step: `0:+x 1:-x 2:+y 3:-y`. -/
def step (p : Site) (d : Fin 4) : Site :=
  if d.val = 0 then (p.1 + 1, p.2)
  else if d.val = 1 then (p.1 - 1, p.2)
  else if d.val = 2 then (p.1, p.2 + 1)
  else (p.1, p.2 - 1)

/-- Manhattan distance exactly one. -/
def manhattan1 (p r : Site) : Bool :=
  decide ((p.1 - r.1).natAbs + (p.2 - r.2).natAbs = 1)

/-! ### Incremental walk state -/

/-- The walk built so far.  `rev` lists the placed sites, most recent first,
each with its hydrophobicity flag; `known` is the number of H–H contacts
among the placed residues. -/
structure State where
  cur : Site
  rev : List (Site × Bool)
  placed : Nat
  known : Nat

/-- Scan the placed sites (most recent first) for the candidate new site `s`.
`gap` is the chain distance from the new residue to the head entry.  Returns
`none` when `s` is occupied, and otherwise the accumulated number of H–H
contacts the new residue makes. -/
def scan (s : Site) (hNew : Bool) : Nat → List (Site × Bool) → Nat → Option Nat
  | _, [], acc => some acc
  | gap, e :: rest, acc =>
      if e.1 = s then none
      else
        scan s hNew (gap + 1) rest
          (if hNew && e.2 && decide (3 ≤ gap) && decide (gap % 2 = 1) && manhattan1 e.1 s
           then acc + 1 else acc)

/-- Extend the walk by direction `d`; `none` when the target site is occupied. -/
def stepState (seq : Array Bool) (st : State) (d : Fin 4) : Option State :=
  let s := step st.cur d
  let hNew := seq.getD st.placed false
  match scan s hNew 1 st.rev 0 with
  | none => none
  | some c =>
      some { cur := s, rev := (s, hNew) :: st.rev, placed := st.placed + 1, known := st.known + c }

/-- The anchored two-site walk `(0,0), (1,0)`. -/
def rootState (seq : Array Bool) : State :=
  { cur := site1
    rev := [(site1, seq.getD 1 false), (site0, seq.getD 0 false)]
    placed := 2
    known := 0 }

/-! ### Exact contact histogram -/

/-- Add one conformation with `k` contacts; fails closed when `k` is outside
the histogram. -/
def bump (h : Array Nat) (k : Nat) : Option (Array Nat) :=
  if hk : k < h.size then some (h.set k (h[k] + 1)) else none

/-- Continue below direction `d` if the target site is free; an occupied
direction contributes nothing. -/
@[inline] def childWith (seq : Array Bool) (st : State) (d : Fin 4) (ex : Array Nat)
    (k : State → Array Nat → Option (Array Nat)) : Option (Array Nat) :=
  match stepState seq st d with
  | none => some ex
  | some st' => k st' ex

/-- Exhaustive depth-first enumeration of all completions (`r` = number of
residues still to place), accumulating the exact contact histogram. -/
def enumGo (seq : Array Bool) : Nat → State → Array Nat → Option (Array Nat)
  | 0, st, ex => bump ex st.known
  | r + 1, st, ex =>
      match childWith seq st 0 ex (enumGo seq r) with
      | none => none
      | some ex =>
        match childWith seq st 1 ex (enumGo seq r) with
        | none => none
        | some ex =>
          match childWith seq st 2 ex (enumGo seq r) with
          | none => none
          | some ex => childWith seq st 3 ex (enumGo seq r)

/-! ### Bound entries -/

/-- Candidate contact pair: both hydrophobic, index gap at least three and odd. -/
def candPair (seq : Array Bool) (i j : Nat) : Bool :=
  decide (i + 3 ≤ j) && decide ((j - i) % 2 = 1) && seq.getD i false && seq.getD j false

/-- Free contact slots of residue `j` (two, or three for the last residue). -/
def slots (n j : Nat) : Nat := if j + 1 < n then 2 else 3

/-- Number of candidate partners `i` of residue `j`. -/
def candCount (seq : Array Bool) (j : Nat) : Nat :=
  ((List.range seq.size).map (fun i => if candPair seq i j then 1 else 0)).sum

/-- Degree-capped number of contacts whose larger index is at least `m`. -/
def futureCap (seq : Array Bool) (m : Nat) : Nat :=
  ((List.range seq.size).map
    (fun j => if m ≤ j then min (slots seq.size j) (candCount seq j) else 0)).sum

/-- Table of `futureCap seq m` for `m = 0 .. seq.size`. -/
def capTable (seq : Array Bool) : Array Nat :=
  Array.ofFn (n := seq.size + 1) (fun m => futureCap seq m.val)

/-- Least `c + r` over rectangles with `c ≥ w` columns, `r ≥ h` rows and at
least `N` cells. -/
def minPerim (N w h : Nat) : Nat :=
  ((List.range (N + 1)).filter (fun c => w ≤ c)).foldr
    (fun c acc => min acc (c + max h ((N + c - 1) / c))) (2 * N + 2)

def minX : List (Site × Bool) → Int
  | [] => 0
  | [e] => e.1.1
  | e :: rest => min e.1.1 (minX rest)

def maxX : List (Site × Bool) → Int
  | [] => 0
  | [e] => e.1.1
  | e :: rest => max e.1.1 (maxX rest)

def minY : List (Site × Bool) → Int
  | [] => 0
  | [e] => e.1.2
  | e :: rest => min e.1.2 (minY rest)

def maxY : List (Site × Bool) → Int
  | [] => 0
  | [e] => e.1.2
  | e :: rest => max e.1.2 (maxY rest)

/-- Coordinate along axis direction `m`. -/
def oF (m : Fin 4) (s : Site) : Int :=
  if m.val = 0 then s.1 else if m.val = 1 then -s.1 else if m.val = 2 then s.2 else -s.2

/-- Diagonal potential `(±x) + (±y)`. -/
def pot (a b : Bool) (s : Site) : Int :=
  (if a then s.1 else -s.1) + (if b then s.2 else -s.2)

/-- The end site is the strict unique extreme of the placed sites along `m`. -/
def pdTest (m : Fin 4) (cur : Site) (rev : List (Site × Bool)) : Bool :=
  rev.all (fun e => e.1 == cur || decide (oF m e.1 < oF m cur))

/-- The end site maximises the diagonal potential over the placed sites. -/
def diagTest (a b : Bool) (cur : Site) (rev : List (Site × Bool)) : Bool :=
  rev.all (fun e => decide (pot a b e.1 ≤ pot a b cur))

/-- Counts `(free, after-up, after-down)` of partially directed words. -/
def pdTriple : Nat → Nat × Nat × Nat
  | 0 => (1, 1, 1)
  | r + 1 =>
      let t := pdTriple r
      (t.1 + t.2.1 + t.2.2, t.1 + t.2.1, t.1 + t.2.2)

/-- Number of partially directed continuations of length `r`. -/
def pdCount (r : Nat) : Nat := (pdTriple r).1

/-- Certified number of completions available without enumeration. -/
def famCount (st : State) (r : Nat) : Nat :=
  if pdTest 0 st.cur st.rev || pdTest 1 st.cur st.rev ||
      pdTest 2 st.cur st.rev || pdTest 3 st.cur st.rev then pdCount r
  else if diagTest true true st.cur st.rev || diagTest true false st.cur st.rev ||
      diagTest false true st.cur st.rev || diagTest false false st.cur st.rev then 2 ^ r
  else 0

/-- A bound entry: every completion has between `a` and `K` contacts, and the
number of completions is between `fam` and `cnt`. -/
structure Entry where
  a : Nat
  K : Nat
  cnt : Nat
  fam : Nat

/-- Recompute the bound entry of the current walk. -/
def boundEntry (seq : Array Bool) (capTab : Array Nat) (st : State) : Option Entry :=
  match capTab[st.placed]? with
  | none => none
  | some fc =>
      let n := seq.size
      let r := n - st.placed
      let w := (maxX st.rev - minX st.rev + 1).toNat
      let h := (maxY st.rev - minY st.rev + 1).toNat
      some { a := st.known
             K := min (st.known + fc) (n + 1 - minPerim n w h)
             cnt := 3 ^ r
             fam := famCount st r }

/-! ### Certificate interpreter -/

/-- Rejection categories.  Positions are byte offsets in the input. -/
inductive Err where
  | badMagic
  | badSequence
  | shortSequence
  | badQuery (pos : Nat)
  | noQuery
  | qBelowOne
  | badTreeLine (pos : Nat)
  | truncated (pos : Nat)
  | badByte (pos : Nat)
  | exactNotTerminal (pos : Nat)
  | splitAtTerminal (pos : Nat)
  | expectedBlocked (pos : Nat)
  | histogramOverflow (pos : Nat)
  | capTableMiss (pos : Nat)
  | badTrailer (pos : Nat)
  deriving Repr, DecidableEq

/-- Everything accumulated by the tree walk.  Only `ex` and `entries` carry
meaning; the counters are instrumentation. -/
structure Acc where
  ex : Array Nat
  entries : List Entry
  splits : Nat
  boundLeaves : Nat
  enumLeaves : Nat
  exactLeaves : Nat

/-- Handle the node encoded at `pos` (everything except the recursion of a
split, which is supplied as `onSplit`). -/
@[inline] def runNodeBody (seq : Array Bool) (capTab : Array Nat) (cert : ByteArray)
    (r : Nat) (st : State) (pos : Nat) (acc : Acc)
    (onSplit : Acc → Except Err (Acc × Nat)) : Except Err (Acc × Nat) :=
  if h : pos < cert.size then
    let b := cert[pos]
    if b = 69 then
      if r = 0 then
        match acc with
        | ⟨ex, entries, splits, boundLeaves, enumLeaves, exactLeaves⟩ =>
          match bump ex st.known with
          | some ex' => .ok (⟨ex', entries, splits, boundLeaves, enumLeaves, exactLeaves + 1⟩, pos + 1)
          | none => .error (.histogramOverflow pos)
      else .error (.exactNotTerminal pos)
    else if b = 88 then
      match acc with
      | ⟨ex, entries, splits, boundLeaves, enumLeaves, exactLeaves⟩ =>
        match enumGo seq r st ex with
        | some ex' => .ok (⟨ex', entries, splits, boundLeaves, enumLeaves + 1, exactLeaves⟩, pos + 1)
        | none => .error (.histogramOverflow pos)
    else if b = 66 then
      match boundEntry seq capTab st with
      | some e =>
        match acc with
        | ⟨ex, entries, splits, boundLeaves, enumLeaves, exactLeaves⟩ =>
          .ok (⟨ex, e :: entries, splits, boundLeaves + 1, enumLeaves, exactLeaves⟩, pos + 1)
      | none => .error (.capTableMiss pos)
    else if b = 78 then onSplit acc
    else .error (.badByte pos)
  else .error (.truncated pos)

/-- Child `d` of a split: an occupied direction must be marked `-`; a free
direction must carry a node, checked by the continuation `k`. -/
@[inline] def runChild (seq : Array Bool) (cert : ByteArray) (st : State) (d : Fin 4)
    (pos : Nat) (acc : Acc)
    (k : State → Nat → Acc → Except Err (Acc × Nat)) : Except Err (Acc × Nat) :=
  match stepState seq st d with
  | none =>
      if h : pos < cert.size then
        if cert[pos] = 45 then .ok (acc, pos + 1) else .error (.expectedBlocked pos)
      else .error (.truncated pos)
  | some st' => k st' pos acc

/-- Walk the certificate tree encoded from byte `pos` on; `r` is the number
of residues still to place.  Returns the updated accumulator and the
position after the subtree. -/
def runTree (seq : Array Bool) (capTab : Array Nat) (cert : ByteArray) :
    Nat → State → Nat → Acc → Except Err (Acc × Nat)
  | 0, st, pos, acc =>
      runNodeBody seq capTab cert 0 st pos acc (fun _ => .error (.splitAtTerminal pos))
  | r + 1, st, pos, acc =>
      runNodeBody seq capTab cert (r + 1) st pos acc (fun acc =>
        match runChild seq cert st 0 (pos + 1) { acc with splits := acc.splits + 1 }
            (runTree seq capTab cert r) with
        | .error e => .error e
        | .ok (acc, pos) =>
          match runChild seq cert st 1 pos acc (runTree seq capTab cert r) with
          | .error e => .error e
          | .ok (acc, pos) =>
            match runChild seq cert st 2 pos acc (runTree seq capTab cert r) with
            | .error e => .error e
            | .ok (acc, pos) => runChild seq cert st 3 pos acc (runTree seq capTab cert r))

/-- Number of histogram slots for a chain of `n` residues. -/
def histSize (n : Nat) : Nat := 2 * n + 2

def emptyAcc (n : Nat) : Acc :=
  { ex := Array.replicate (histSize n) 0
    entries := []
    splits := 0
    boundLeaves := 0
    enumLeaves := 0
    exactLeaves := 0 }

/-- Check a whole tree for the sequence `seq`, starting at byte `start`.
Returns the accumulator and the position after the tree. -/
def runRoot (seq : Array Bool) (cert : ByteArray) (start : Nat) : Except Err (Acc × Nat) :=
  if seq.size < 2 then .error .shortSequence
  else runTree seq (capTable seq) cert (seq.size - 2) (rootState seq) start (emptyAcc seq.size)

/-! ### Evaluation for a Boltzmann factor `q` and a threshold `t` -/

/-- Exactly counted event mass `∑_{k ≥ t} ex[k] q^k`. -/
def histE (q : Rat) (t : Nat) (ex : Array Nat) : Rat :=
  ((List.range ex.size).map (fun k => if t ≤ k then (ex.getD k 0 : Rat) * q ^ k else 0)).sum

/-- Exactly counted complement mass `∑_{k < t} ex[k] q^k`. -/
def histN (q : Rat) (t : Nat) (ex : Array Nat) : Rat :=
  ((List.range ex.size).map (fun k => if k < t then (ex.getD k 0 : Rat) * q ^ k else 0)).sum

def entZLo (q : Rat) (e : Entry) : Rat := (e.fam : Rat) * q ^ e.a
def entZHi (q : Rat) (e : Entry) : Rat := (e.cnt : Rat) * q ^ e.K
def entELo (q : Rat) (t : Nat) (e : Entry) : Rat := if t ≤ e.a then (e.fam : Rat) * q ^ e.a else 0
def entEHi (q : Rat) (t : Nat) (e : Entry) : Rat := if t ≤ e.K then (e.cnt : Rat) * q ^ e.K else 0
def entNLo (q : Rat) (t : Nat) (e : Entry) : Rat := if e.K < t then (e.fam : Rat) * q ^ e.a else 0
def entNHi (q : Rat) (t : Nat) (e : Entry) : Rat :=
  if e.a < t then (e.cnt : Rat) * q ^ (min e.K (t - 1)) else 0

/-- The result reported for one query. -/
structure QueryResult where
  q : Rat
  t : Nat
  zLo : Rat
  zHi : Rat
  eLo : Rat
  eHi : Rat
  nLo : Rat
  nHi : Rat
  pLo : Rat
  pHi : Rat

/-- Exact rational enclosures for the pair `(q, t)`. -/
def evalQuery (acc : Acc) (q : Rat) (t : Nat) : QueryResult :=
  let hE := histE q t acc.ex
  let hN := histN q t acc.ex
  let zLo := hE + hN + (acc.entries.map (entZLo q)).sum
  let zHi := hE + hN + (acc.entries.map (entZHi q)).sum
  let eLo := hE + (acc.entries.map (entELo q t)).sum
  let eHi := hE + (acc.entries.map (entEHi q t)).sum
  let nLo := hN + (acc.entries.map (entNLo q t)).sum
  let nHi := hN + (acc.entries.map (entNHi q t)).sum
  { q := q, t := t, zLo := zLo, zHi := zHi, eLo := eLo, eHi := eHi, nLo := nLo, nHi := nHi
    pLo := if eLo = 0 then 0 else eLo / (eLo + nHi)
    pHi := if eHi = 0 then 0 else if nLo = 0 then 1 else eHi / (eHi + nLo) }

/-! ### Strict text front end -/

/-- ASCII bytes of a literal. -/
def lit (s : String) : List UInt8 := s.toList.map (fun c => c.toNat.toUInt8)

/-- Expect the bytes `l` at `pos`; returns the position after them. -/
def expectLit (b : ByteArray) : List UInt8 → Nat → Option Nat
  | [], pos => some pos
  | c :: cs, pos =>
      if h : pos < b.size then
        if b[pos] = c then expectLit b cs (pos + 1) else none
      else none

/-- Read `H`/`P` bytes up to (not including) a newline. -/
def readSeq (b : ByteArray) : Nat → Nat → Array Bool → Option (Array Bool × Nat)
  | 0, _, _ => none
  | fuel + 1, pos, acc =>
      if h : pos < b.size then
        if b[pos] = 72 then readSeq b fuel (pos + 1) (acc.push true)
        else if b[pos] = 80 then readSeq b fuel (pos + 1) (acc.push false)
        else if b[pos] = 10 then some (acc, pos)
        else none
      else none

def isDigit (c : UInt8) : Bool := 48 ≤ c && c ≤ 57

/-- Accumulate decimal digits. -/
def readDigits (b : ByteArray) : Nat → Nat → Nat → Nat × Nat
  | 0, pos, acc => (acc, pos)
  | fuel + 1, pos, acc =>
      if h : pos < b.size then
        if isDigit b[pos] then readDigits b fuel (pos + 1) (10 * acc + (b[pos].toNat - 48))
        else (acc, pos)
      else (acc, pos)

/-- Canonical natural number: `0`, or digits without a leading zero. -/
def readNat (b : ByteArray) (pos : Nat) : Option (Nat × Nat) :=
  if h : pos < b.size then
    if b[pos] = 48 then
      if h2 : pos + 1 < b.size then
        if isDigit b[pos + 1] then none else some (0, pos + 1)
      else some (0, pos + 1)
    else if isDigit b[pos] then some (readDigits b (b.size - pos) pos 0)
    else none
  else none

/-- Canonical nonnegative rational: a canonical natural, or `num/den` in
lowest terms with `den ≥ 2`. -/
def readRat (b : ByteArray) (pos : Nat) : Option (Rat × Nat) :=
  match readNat b pos with
  | none => none
  | some (num, p1) =>
      if h : p1 < b.size then
        if b[p1] = 47 then
          match readNat b (p1 + 1) with
          | none => none
          | some (den, p2) =>
              if 2 ≤ den ∧ Nat.gcd num den = 1 then some ((num : Rat) / (den : Rat), p2)
              else none
        else some ((num : Rat), p1)
      else some ((num : Rat), p1)

/-- Read `query <q> <t>` lines until a line does not start with `query `. -/
def readQueries (b : ByteArray) : Nat → Nat → List (Rat × Nat) →
    Except Err (List (Rat × Nat) × Nat)
  | 0, pos, _ => .error (.badQuery pos)
  | fuel + 1, pos, acc =>
      match expectLit b (lit "query ") pos with
      | none => .ok (acc.reverse, pos)
      | some p1 =>
        match readRat b p1 with
        | none => .error (.badQuery p1)
        | some (q, p2) =>
          match expectLit b [32] p2 with
          | none => .error (.badQuery p2)
          | some p3 =>
            match readNat b p3 with
            | none => .error (.badQuery p3)
            | some (t, p4) =>
              match expectLit b [10] p4 with
              | none => .error (.badQuery p4)
              | some p5 => readQueries b fuel p5 ((q, t) :: acc)

/-- A parsed input: model parameters, queries, and where the tree starts. -/
structure Doc where
  seq : Array Bool
  queries : List (Rat × Nat)
  bytes : ByteArray
  treeStart : Nat

/-- Parse everything up to the start of the tree. -/
def parseHeader (b : ByteArray) : Except Err Doc :=
  match expectLit b (lit "FOLDBOUNDS-KERNEL 1\nsequence ") 0 with
  | none => .error .badMagic
  | some p0 =>
    match readSeq b (b.size + 1) p0 #[] with
    | none => .error .badSequence
    | some (seq, p1) =>
      match readQueries b (b.size + 1) (p1 + 1) [] with
      | .error e => .error e
      | .ok (queries, p2) =>
        if queries.isEmpty then .error .noQuery
        else
          match expectLit b (lit "tree ") p2 with
          | none => .error (.badTreeLine p2)
          | some p3 => .ok { seq := seq, queries := queries, bytes := b, treeStart := p3 }

/-- The verified report: the echoed model, one result per query, and
instrumentation counters. -/
structure Report where
  seq : List Bool
  results : List QueryResult
  exactHist : List Nat
  entries : Nat
  treeBytes : Nat
  splits : Nat
  boundLeaves : Nat
  enumLeaves : Nat
  exactLeaves : Nat

/-- Check a parsed document. -/
def checkDoc (doc : Doc) : Except Err Report :=
  if doc.queries.any (fun qt => decide (qt.1 < 1)) then .error .qBelowOne
  else
    match runRoot doc.seq doc.bytes doc.treeStart with
    | .error e => .error e
    | .ok (acc, stop) =>
      if expectLit doc.bytes (lit "\nend\n") stop = some doc.bytes.size then
        .ok { seq := doc.seq.toList
              results := doc.queries.map (fun qt => evalQuery acc qt.1 qt.2)
              exactHist := acc.ex.toList
              entries := acc.entries.length
              treeBytes := stop - doc.treeStart
              splits := acc.splits
              boundLeaves := acc.boundLeaves
              enumLeaves := acc.enumLeaves
              exactLeaves := acc.exactLeaves }
      else .error (.badTrailer stop)

/-- **The checker**, as a function of the raw input bytes. -/
def verifyBytes (b : ByteArray) : Except Err Report :=
  match parseHeader b with
  | .error e => .error e
  | .ok doc => checkDoc doc

/-- **The checker**, as a function of the raw input string. -/
def verify (input : String) : Except Err Report := verifyBytes input.toUTF8

end Kernel
end FoldBounds
