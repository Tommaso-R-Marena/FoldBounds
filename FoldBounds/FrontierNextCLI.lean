module

public import FoldBounds.FrontierNextFlat
public import FoldBounds.FrontierCLI

/-!
# Text front end for the version-2 checker (UNVERIFIED boundary)

No theorems.  Parses a canonical line-based certificate into `(S, q, thr)`
and a flat version-2 frontier.  Soundness of every printed interval
relative to the *parsed* parameters is `FrontierNext.checkFlat2_sound`;
that the parsed parameters are the intended ones is part of the trusted
boundary (they are echoed back).

    FOLDBOUNDS-FRONTIER 2            (or 1: only E/D/B/T lines allowed)
    sequence <H|P string, length ≥ 2>
    q <natural> | <natural>/<positive natural>
    threshold <natural>
    <KIND> <WORD> [<WITNESS-WORD> ...]
    ...
    end

Kinds: `E` exact, `D` dead, `X` exhaustive local evaluation (version 2),
`B`/`T` legacy basic/transfer bound, `M<k>` memory-`k` transfer bound
with geometric caps (version 2, e.g. `M4`).  Witness words only after
`B`, `T`, `M<k>`.  Anything else is rejected.
-/
@[expose] public section

namespace FoldBounds
namespace FrontierNextCLI

open FrontierCheck FrontierNext FrontierCLI LatticeModel

/-- Parse one node line; `v2` enables the version-2 kinds. -/
def parseNode2 (v2 : Bool) (line : String) : Option (List (Fin 4) × NTree) :=
  match line.splitOn " " with
  | ["E", w] => (parseWord w).map (fun p => (p, .exact))
  | ["D", w] => (parseWord w).map (fun p => (p, .dead))
  | ["X", w] => if v2 then (parseWord w).map (fun p => (p, .enum)) else none
  | k :: w :: ws =>
      match parseWord w, ws.mapM parseWord with
      | some p, some wit =>
          if k = "B" then some (p, .legacy false wit)
          else if k = "T" then some (p, .legacy true wit)
          else if v2 ∧ k.startsWith "M" then
            (parseNatStrict (k.drop 1).toString).map (fun m => (p, .mem m wit))
          else none
      | _, _ => none
  | _ => none

structure Parsed2 where
  version : ℕ
  seq : List Bool
  q : ℚ
  thr : ℕ
  nodes : FlatCert2

def parseCert2 (txt : String) : Except String Parsed2 := do
  let lines := (txt.splitOn "\n").filter (fun l => !l.isEmpty)
  match lines with
  | hdr :: sl :: ql :: tl :: rest =>
      let v ← (if hdr = "FOLDBOUNDS-FRONTIER 2" then pure 2
        else if hdr = "FOLDBOUNDS-FRONTIER 1" then pure 1 else throw "bad header")
      let some seq := (match sl.splitOn " " with
        | ["sequence", s] => parseSeq s | _ => none) | throw "bad sequence line"
      let some q := (match ql.splitOn " " with
        | ["q", s] => parseRat s | _ => none) | throw "bad q line"
      let some thr := (match tl.splitOn " " with
        | ["threshold", s] => parseNatStrict s | _ => none) | throw "bad threshold line"
      match rest.reverse with
      | "end" :: body =>
          let some nodes := body.reverse.mapM (parseNode2 (v = 2)) | throw "bad node line"
          return ⟨v, seq, q, thr, nodes⟩
      | _ => throw "missing end line"
  | _ => throw "truncated certificate"

def resultJson2 (P : Parsed2) (r : Option FrontierResult) : String :=
  let echo := s!"\"version\":{P.version},\"sequence\":\"{seqStr P.seq}\",\"q\":\"{ratStr P.q}\"," ++
    s!"\"threshold\":{P.thr},\"nodes\":{P.nodes.length}"
  match r with
  | none => "{" ++ echo ++ ",\"accepted\":false}"
  | some r =>
      "{" ++ echo ++ s!",\"accepted\":true,\"zLo\":\"{ratStr r.zLo}\",\"zHi\":\"{ratStr r.zHi}\"" ++
      s!",\"eLo\":\"{ratStr r.eLo}\",\"eHi\":\"{ratStr r.eHi}\"" ++
      s!",\"nLo\":\"{ratStr r.nLo}\",\"nHi\":\"{ratStr r.nHi}\"" ++
      s!",\"pLo\":\"{ratStr r.pLo}\",\"pHi\":\"{ratStr r.pHi}\"" ++ "}"

end FrontierNextCLI
end FoldBounds
