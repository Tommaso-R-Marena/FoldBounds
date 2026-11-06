import FoldBounds.FrontierCLI

/-! Command-line runner (unverified I/O shell around the verified checker).

    lake env lean --run scripts/foldcheck.lean check CERT.txt [CERT2.txt ...]
    lake env lean --run scripts/foldcheck.lean enumerate HHPH...

`check` prints one JSON object per certificate; `enumerate` prints all
valid conformations (Lean decoding) with reference contact counts. -/

open FoldBounds FoldBounds.FrontierCheck FoldBounds.FrontierCLI

def main (args : List String) : IO UInt32 := do
  match args with
  | "check" :: files =>
      for f in files do
        let txt ← IO.FS.readFile f
        match parseCert txt with
        | .error e => IO.println ("{\"file\":\"" ++ f ++ "\",\"parse_error\":\"" ++ e ++ "\"}")
        | .ok P =>
            let t0 ← IO.monoMsNow
            let r := checkFlat P.seq P.q P.thr P.nodes
            let s := resultJson P r
            -- print first (forces evaluation), then report the elapsed time
            IO.print ("{\"file\":\"" ++ f ++ "\"," ++ ((s.drop 1).dropEnd 1).toString)
            let t1 ← IO.monoMsNow
            IO.println (",\"verify_ms\":" ++ toString (t1 - t0) ++ "}")
      return 0
  | ["enumerate", seq] =>
      match parseSeq seq with
      | some S => IO.println (enumerateJson S); return 0
      | none => IO.eprintln "bad sequence"; return 1
  | _ => IO.eprintln "usage: check FILES... | enumerate SEQ"; return 1
