module

public import FoldBounds.FrontierFlat

/-!
# Text front end for the frontier checker (UNVERIFIED boundary)

This module contains *no theorems*.  It parses a canonical line-based text
certificate into the model parameters `(S, q, thr)` and a flat frontier,
and prints results.  Soundness of every printed interval relative to the
*parsed* parameters is guaranteed by `checkFlat_sound`; that the parsed
parameters are the ones the user intended is part of the trusted
boundary (the parser echoes them back so they can be inspected).

Format (one item per line, tokens separated by single spaces):

    FOLDBOUNDS-FRONTIER 1
    sequence <H|P string, length ≥ 2>
    q <natural> | <natural>/<positive natural>
    threshold <natural>
    <KIND> <WORD> [<WITNESS-WORD> ...]
    ...
    end

`KIND` is `E` (exact), `B` (basic bound), `T` (transfer bound) or `D`
(dead). Words are strings over `0123` (Lean direction encoding
`0:+x 1:-x 2:+y 3:-y`), with `.` for the empty word. Witness words are
only allowed after `B`/`T`. Anything else is rejected.
-/
@[expose] public section

namespace FoldBounds
namespace FrontierCLI

open FrontierCheck LatticeModel

def parseDir (c : Char) : Option (Fin 4) :=
  if c = '0' then some 0 else if c = '1' then some 1
  else if c = '2' then some 2 else if c = '3' then some 3 else none

def parseWord (s : String) : Option (List (Fin 4)) :=
  if s = "." then some [] else
  if s.isEmpty then none else s.toList.mapM parseDir

def parseSeq (s : String) : Option (List Bool) :=
  let r := s.toList.mapM (fun c => if c = 'H' then some true else if c = 'P' then some false
    else none)
  match r with
  | some l => if 2 ≤ l.length then some l else none
  | none => none

def parseNatStrict (s : String) : Option ℕ :=
  if s.isEmpty || !(s.all Char.isDigit) then none else s.toNat?

def parseRat (s : String) : Option ℚ :=
  match s.splitOn "/" with
  | [a] => (parseNatStrict a).map (fun n => (n : ℚ))
  | [a, b] =>
      match parseNatStrict a, parseNatStrict b with
      | some n, some d => if d = 0 then none else some ((n : ℚ) / (d : ℚ))
      | _, _ => none
  | _ => none

def parseNode (line : String) : Option (List (Fin 4) × NodeKind) :=
  match line.splitOn " " with
  | ["E", w] => (parseWord w).map (fun p => (p, .exact))
  | ["D", w] => (parseWord w).map (fun p => (p, .dead))
  | k :: w :: ws =>
      if k = "B" ∨ k = "T" then
        match parseWord w, ws.mapM parseWord with
        | some p, some wit => some (p, .bound (k = "T") wit)
        | _, _ => none
      else none
  | _ => none

structure Parsed where
  seq : List Bool
  q : ℚ
  thr : ℕ
  nodes : FlatCert

def parseCert (txt : String) : Except String Parsed := do
  let lines := (txt.splitOn "\n").filter (fun l => !l.isEmpty)
  match lines with
  | hdr :: sl :: ql :: tl :: rest =>
      if hdr ≠ "FOLDBOUNDS-FRONTIER 1" then throw "bad header"
      let some seq := (match sl.splitOn " " with
        | ["sequence", s] => parseSeq s | _ => none) | throw "bad sequence line"
      let some q := (match ql.splitOn " " with
        | ["q", s] => parseRat s | _ => none) | throw "bad q line"
      let some thr := (match tl.splitOn " " with
        | ["threshold", s] => parseNatStrict s | _ => none) | throw "bad threshold line"
      match rest.reverse with
      | "end" :: body =>
          let some nodes := body.reverse.mapM parseNode | throw "bad node line"
          return ⟨seq, q, thr, nodes⟩
      | _ => throw "missing end line"
  | _ => throw "truncated certificate"

def ratStr (x : ℚ) : String := s!"{x.num}/{x.den}"

def seqStr (S : List Bool) : String := String.ofList (S.map (fun b => if b then 'H' else 'P'))

def resultJson (P : Parsed) (r : Option FrontierResult) : String :=
  let echo := s!"\"sequence\":\"{seqStr P.seq}\",\"q\":\"{ratStr P.q}\",\"threshold\":{P.thr}," ++
    s!"\"nodes\":{P.nodes.length}"
  match r with
  | none => "{" ++ echo ++ ",\"accepted\":false}"
  | some r =>
      "{" ++ echo ++ s!",\"accepted\":true,\"zLo\":\"{ratStr r.zLo}\",\"zHi\":\"{ratStr r.zHi}\"" ++
      s!",\"eLo\":\"{ratStr r.eLo}\",\"eHi\":\"{ratStr r.eHi}\"" ++
      s!",\"nLo\":\"{ratStr r.nLo}\",\"nHi\":\"{ratStr r.nHi}\"" ++
      s!",\"pLo\":\"{ratStr r.pLo}\",\"pHi\":\"{ratStr r.pHi}\"" ++ "}"

/-- All direction words of length `n` (independent Lean enumeration used
only for differential testing). -/
def allWords : ℕ → List (List (Fin 4))
  | 0 => [[]]
  | n + 1 => (allWords n).flatMap (fun w => [0, 1, 2, 3].map (fun d => w ++ [d]))

/-- Valid conformations with their decoded sites and reference contact
counts, for differential testing against Python. -/
def enumerateJson (S : List Bool) : String :=
  let n := S.length - 2
  let ws := (allWords n).filter (fun w => (walkFrom w).Nodup)
  let items := ws.map (fun w =>
    let sites := (walkFrom w).map (fun s => s!"[{s.1},{s.2}]")
    "{\"sites\":[" ++ ",".intercalate sites ++ "],\"contacts\":" ++
      toString (contactCount S (walkFrom w)) ++ "}")
  "[" ++ ",".intercalate items ++ "]"

end FrontierCLI
end FoldBounds
