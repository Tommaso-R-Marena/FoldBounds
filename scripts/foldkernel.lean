import FoldBounds.FrontierKernel

/-! Command-line runner for the verified frontier kernel
(`FoldBounds.Kernel.verifyBytes`).  This file is an **unverified I/O shell**:
it reads files, calls the verified function on the raw bytes, and prints the
returned report as JSON.  It imports nothing but the Mathlib-free kernel.

    foldkernel check CERT [CERT ...]

One JSON object is printed per certificate.  Exit status: `0` if every
certificate was accepted, `1` if some certificate was rejected, `2` on usage
or I/O errors. -/

open FoldBounds.Kernel

def ratStr (x : Rat) : String := s!"{x.num}/{x.den}"

def seqStr (S : List Bool) : String := String.ofList (S.map (fun b => if b then 'H' else 'P'))

def jsonEscape (s : String) : String :=
  String.join (s.toList.map (fun c =>
    if c = '"' then "\\\"" else if c = '\\' then "\\\\"
    else if c.toNat < 32 then "?" else c.toString))

def errJson : Err → String
  | .badMagic => "\"reason\":\"bad-magic\""
  | .badSequence => "\"reason\":\"bad-sequence\""
  | .shortSequence => "\"reason\":\"short-sequence\""
  | .badQuery pos => s!"\"reason\":\"bad-query\",\"pos\":{pos}"
  | .noQuery => "\"reason\":\"no-query\""
  | .qBelowOne => "\"reason\":\"q-below-one\""
  | .badTreeLine pos => s!"\"reason\":\"bad-tree-line\",\"pos\":{pos}"
  | .truncated pos => s!"\"reason\":\"tree-truncated\",\"pos\":{pos}"
  | .badByte pos => s!"\"reason\":\"tree-bad-byte\",\"pos\":{pos}"
  | .exactNotTerminal pos => s!"\"reason\":\"exact-not-terminal\",\"pos\":{pos}"
  | .splitAtTerminal pos => s!"\"reason\":\"split-at-terminal\",\"pos\":{pos}"
  | .expectedBlocked pos => s!"\"reason\":\"expected-blocked-marker\",\"pos\":{pos}"
  | .histogramOverflow pos => s!"\"reason\":\"histogram-overflow\",\"pos\":{pos}"
  | .capTableMiss pos => s!"\"reason\":\"cap-table-miss\",\"pos\":{pos}"
  | .badTrailer pos => s!"\"reason\":\"bad-trailer\",\"pos\":{pos}"

def resultJson (r : QueryResult) : String :=
  "{" ++ s!"\"q\":\"{ratStr r.q}\",\"threshold\":{r.t}" ++
  s!",\"zLo\":\"{ratStr r.zLo}\",\"zHi\":\"{ratStr r.zHi}\"" ++
  s!",\"eLo\":\"{ratStr r.eLo}\",\"eHi\":\"{ratStr r.eHi}\"" ++
  s!",\"nLo\":\"{ratStr r.nLo}\",\"nHi\":\"{ratStr r.nHi}\"" ++
  s!",\"pLo\":\"{ratStr r.pLo}\",\"pHi\":\"{ratStr r.pHi}\"" ++ "}"

def trimHist (h : List Nat) : List Nat := (h.reverse.dropWhile (· == 0)).reverse

def reportJson (rep : Report) : String :=
  s!"\"accepted\":true,\"format\":\"FOLDBOUNDS-KERNEL 1\"" ++
  s!",\"sequence\":\"{seqStr rep.seq}\",\"n\":{rep.seq.length}" ++
  s!",\"cert_bytes\":{rep.treeBytes},\"split_nodes\":{rep.splits}" ++
  s!",\"bound_leaves\":{rep.boundLeaves},\"enum_leaves\":{rep.enumLeaves}" ++
  s!",\"exact_leaves\":{rep.exactLeaves},\"exact_count\":{rep.exactHist.foldl (· + ·) 0}" ++
  s!",\"exact_histogram\":[{",".intercalate ((trimHist rep.exactHist).map toString)}]" ++
  s!",\"results\":[{",".intercalate (rep.results.map resultJson)}]"

def main (args : List String) : IO UInt32 := do
  match args with
  | "check" :: files =>
      if files.isEmpty then
        IO.eprintln "usage: foldkernel check CERT [CERT ...]"
        return 2
      let mut status : UInt32 := 0
      for f in files do
        match ← (IO.FS.readBinFile f).toBaseIO with
        | .error e =>
            IO.eprintln s!"cannot read {f}: {e}"
            return 2
        | .ok bytes =>
            let t0 ← IO.monoNanosNow
            match verifyBytes bytes with
            | .error e =>
                let t1 ← IO.monoNanosNow
                status := 1
                IO.println ("{" ++
                  s!"\"file\":\"{jsonEscape f}\",\"accepted\":false,{errJson e}" ++
                  s!",\"verify_ms\":{(t1 - t0) / 1000000}" ++ "}")
            | .ok rep =>
                -- Rendering forces every reported number before the clock is read.
                let body := reportJson rep
                let t1 ← IO.monoNanosNow
                IO.println ("{" ++ s!"\"file\":\"{jsonEscape f}\",{body}" ++
                  s!",\"verify_ms\":{(t1 - t0) / 1000000},\"verify_us\":{(t1 - t0) / 1000}" ++
                  "}")
      return status
  | _ =>
      IO.eprintln "usage: foldkernel check CERT [CERT ...]"
      return 2
