import FoldBounds.FrontierNextCLI

/-! Command-line runner for the version-2 checker (unverified I/O shell
around the verified `FrontierNext.checkFlat2`).  Accepts version-1 and
version-2 text certificates.

    foldcheck2 check CERT.txt [CERT2.txt ...]

Prints one JSON object per certificate. -/

open FoldBounds FoldBounds.FrontierNext FoldBounds.FrontierNextCLI

def main (args : List String) : IO UInt32 := do
  match args with
  | "check" :: files =>
      for f in files do
        let txt ← IO.FS.readFile f
        match parseCert2 txt with
        | .error e => IO.println ("{\"file\":\"" ++ f ++ "\",\"parse_error\":\"" ++ e ++ "\"}")
        | .ok P =>
            let t0 ← IO.monoMsNow
            let r := checkFlat2 P.seq P.q P.thr P.nodes
            let s := resultJson2 P r
            IO.print ("{\"file\":\"" ++ f ++ "\"," ++ ((s.drop 1).dropEnd 1).toString)
            let t1 ← IO.monoMsNow
            IO.println (",\"verify_ms\":" ++ toString (t1 - t0) ++ "}")
      return 0
  | _ => IO.eprintln "usage: foldcheck2 check FILES..."; return 1
