#!/usr/bin/env python3
"""Regenerate FoldBounds/FrontierKernelTests.lean.

The generated Lean file contains only statements that the Lean *kernel*
re-checks (`decide +kernel`): this script merely proposes inputs and the
expected verdicts (obtained by running the compiled `foldkernel` binary).
A wrong expectation makes the Lean file fail to compile, so nothing here is
trusted.

Usage (from the repository root, after `lake build foldkernel`):

    python3 scripts/gen_frontier_kernel_tests.py > FoldBounds/FrontierKernelTests.lean
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

BINARY = Path(".lake/build/bin/foldkernel")
DIRS = [(1, 0), (-1, 0), (0, 1), (0, -1)]
HEADER = "FOLDBOUNDS-KERNEL 1\n"


def known(seq: str, sites: list[tuple[int, int]]) -> int:
    total = 0
    for i in range(len(sites)):
        for j in range(i + 3, len(sites)):
            if seq[i] == "H" and seq[j] == "H":
                if abs(sites[i][0] - sites[j][0]) + abs(sites[i][1] - sites[j][1]) == 1:
                    total += 1
    return total


def tree(seq: str, leaf, sites=None) -> str:
    """Preorder encoding; `leaf(sites)` returns 'B', 'X' or None (split)."""
    if sites is None:
        sites = [(0, 0), (1, 0)]
    choice = leaf(sites)
    if len(sites) == len(seq):
        return choice if choice in ("B", "X") else "E"
    if choice:
        return choice
    out = ["N"]
    x, y = sites[-1]
    for dx, dy in DIRS:
        site = (x + dx, y + dy)
        out.append("-" if site in sites else tree(seq, leaf, sites + [site]))
    return "".join(out)


def doc(seq: str, queries: list[tuple[str, int]], encoded: str) -> str:
    lines = "".join(f"query {q} {t}\n" for q, t in queries)
    return f"{HEADER}sequence {seq}\n{lines}tree {encoded}\nend\n"


def run(data: bytes) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".cert", delete=False) as handle:
        handle.write(data)
        name = handle.name
    result = subprocess.run([str(BINARY), "check", name], capture_output=True, text=True)
    Path(name).unlink()
    if result.returncode not in (0, 1):
        raise SystemExit(f"checker failed on generated input: {result.stderr}")
    return json.loads(result.stdout.strip().splitlines()[-1])


def lean_string(text: str) -> str:
    out = []
    for ch in text:
        if ch == "\n":
            out.append("\\n")
        elif ch == "\r":
            out.append("\\r")
        elif ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        else:
            out.append(ch)
    return '"' + "".join(out) + '"'


def lean_rat(text: str) -> str:
    num, den = text.split("/")
    return num if den == "1" else f"{num}/{den}"


REASONS = {
    "bad-magic": ".badMagic",
    "bad-sequence": ".badSequence",
    "short-sequence": ".shortSequence",
    "bad-query": ".badQuery {pos}",
    "no-query": ".noQuery",
    "q-below-one": ".qBelowOne",
    "bad-tree-line": ".badTreeLine {pos}",
    "tree-truncated": ".truncated {pos}",
    "tree-bad-byte": ".badByte {pos}",
    "exact-not-terminal": ".exactNotTerminal {pos}",
    "split-at-terminal": ".splitAtTerminal {pos}",
    "expected-blocked-marker": ".expectedBlocked {pos}",
    "histogram-overflow": ".histogramOverflow {pos}",
    "cap-table-miss": ".capTableMiss {pos}",
    "bad-trailer": ".badTrailer {pos}",
}


def accept_case(name: str, comment: str, text: str) -> str:
    report = run(text.encode())
    if not report["accepted"]:
        raise SystemExit(f"expected acceptance for {name}: {report}")
    rows = []
    for res in report["results"]:
        keys = ("zLo", "zHi", "eLo", "eHi", "nLo", "nHi", "pLo", "pHi")
        rows.append("[" + ", ".join(lean_rat(res[k]) for k in keys) + "]")
    return (
        f"/-- {comment} -/\n"
        f"theorem {name} :\n"
        f"    endpoints {lean_string(text)}\n"
        f"      = some [{', '.join(rows)}] := by decide +kernel\n"
    )


def reject_case(name: str, comment: str, text: str) -> str:
    report = run(text.encode())
    if report["accepted"]:
        raise SystemExit(f"expected rejection for {name}: {report}")
    err = REASONS[report["reason"]].format(pos=report.get("pos", 0))
    return (
        f"/-- {comment} -/\n"
        f"theorem {name} :\n"
        f"    rejection {lean_string(text)}\n"
        f"      = some ({err}) := by decide +kernel\n"
    )


def reject_bytes_case(name: str, comment: str, data: bytes) -> str:
    report = run(data)
    if report["accepted"]:
        raise SystemExit(f"expected rejection for {name}: {report}")
    err = REASONS[report["reason"]].format(pos=report.get("pos", 0))
    listing = ", ".join(str(b) for b in data)
    return (
        f"/-- {comment} -/\n"
        f"theorem {name} :\n"
        f"    rejectionBytes ⟨#[{listing}]⟩\n"
        f"      = some ({err}) := by decide +kernel\n"
    )


def main() -> None:
    if not BINARY.exists():
        raise SystemExit("build the checker first: lake build foldkernel")
    hhhh_full = tree("HHHH", lambda s: None)
    h8 = "HHHHHHHH"
    h8_hybrid = tree(h8, lambda s: ("X" if known(h8, s) >= 1 else "B") if len(s) >= 5 else None)
    dead_path = [(0, 0), (1, 0), (1, 1), (1, 2), (0, 2), (-1, 2), (-1, 1), (0, 1)]
    h9 = "HHHHHHHHH"
    h9_dead = tree(h9, lambda s: None if s == dead_path[: len(s)] else "B")
    mixed = "HPHHPHPH"
    mixed_tree = tree(mixed, lambda s: "B" if len(s) >= 6 else None)

    out: list[str] = []
    out.append(
        '''import FoldBounds.FrontierKernelSound

/-!
# Kernel-checked tests of the executable frontier kernel

GENERATED by `scripts/gen_frontier_kernel_tests.py`; do not edit by hand.

Every statement below is decided by the Lean **kernel** (`decide +kernel`):
the kernel itself evaluates the string-level checker `Kernel.verify` (header
parser, tree interpreter, enumeration, bound entries and rational
evaluation) on the literal input and compares the result with the stated
value.  No `native_decide`, no compiler trust, no extra axioms.  The
expected values were proposed by running the compiled binary, so these
tests also check that the compiled code and the kernel agree.

The last section turns Boolean checks into theorems about the Lean HP model
through `certifiesExact_sound` and `certifiesInterval_sound`.
-/

open FoldBounds FoldBounds.Kernel FoldBounds.LatticeModel

namespace FoldBounds.FrontierKernelTests

/-- The eight endpoints `zLo zHi eLo eHi nLo nHi pLo pHi` of each reported
result, or `none` when the input is rejected. -/
def endpoints (s : String) : Option (List (List Rat)) :=
  match verify s with
  | .ok rep =>
      some (rep.results.map
        (fun r => [r.zLo, r.zHi, r.eLo, r.eHi, r.nLo, r.nHi, r.pLo, r.pHi]))
  | .error _ => none

/-- The rejection category, or `none` when the input is accepted. -/
def rejection (s : String) : Option Err :=
  match verify s with
  | .ok _ => none
  | .error e => some e

/-- As `rejection`, on raw bytes (inputs that are not valid UTF-8). -/
def rejectionBytes (b : ByteArray) : Option Err :=
  match verifyBytes b with
  | .ok _ => none
  | .error e => some e

/-! ### Accepted certificates -/
'''
    )
    out.append(accept_case(
        "accept_full_exact",
        "Complete explicit tree for `HHHH`: 9 conformations, 2 with one contact; "
        "`q = 2`, `t = 1`: `Z = 11`, `Z_A = 4`, `P = 4/11`.",
        doc("HHHH", [("2", 1)], hhhh_full)))
    out.append(accept_case(
        "accept_root_enum_three_queries",
        "One `X` leaf at the root answers three `(q, t)` queries exactly.",
        doc("HHHH", [("2", 1), ("3/2", 1), ("1", 0)], "X")))
    out.append(accept_case(
        "accept_root_bound",
        "One `B` leaf at the root: `Z ≤ 3^2 · 2^1 = 18`, family lower bound `7`, and the "
        "complement cap uses the exponent `min(K, t - 1) = 0`, giving `9` (the version-1 "
        "basic leaf gives `18`).",
        doc("HHHH", [("2", 1)], "B")))
    out.append(accept_case(
        "accept_split_with_bounds",
        "A split with a blocked marker and three bound leaves.",
        doc("HHHH", [("2", 1)], "NB-BB")))
    out.append(accept_case(
        "accept_threshold_above_range",
        "A threshold above every contact count: event mass and probability are exactly 0.",
        doc("HHHH", [("2", 5)], "X")))
    out.append(accept_case(
        "accept_two_residues",
        "The shortest chain: the root is already complete, `E` is accepted.",
        doc("HP", [("3", 0), ("3", 1)], "E")))
    out.append(accept_case(
        "accept_hybrid_h8",
        "`HHHHHHHH`, depth-5 frontier: prefixes with a known contact are enumerated (`X`), "
        "the others are bounded (`B`); two queries on one tree.",
        doc(h8, [("2", 2), ("3/2", 1)], h8_hybrid)))
    out.append(accept_case(
        "accept_dead_end",
        "A certificate descending to the shortest dead end (8 placed sites): the dead "
        "prefix is encoded as `N----`.",
        doc(h9, [("2", 2)], h9_dead)))
    out.append(accept_case(
        "accept_mixed_sequence",
        "A mixed H/P sequence with a depth-6 bound frontier.",
        doc(mixed, [("2", 1), ("7/3", 2)], mixed_tree)))

    out.append("/-! ### Rejected certificates and malformed inputs -/\n")
    good = doc("HHHH", [("2", 1)], hhhh_full)
    cases = [
        ("reject_omitted_child", "An omitted fourth child of a split.",
         doc("HHHH", [("2", 1)], "NB-B")),
        ("reject_blocked_marker_on_free_direction",
         "A `-` marker on a free direction (an omitted legal branch).",
         doc("HHHH", [("2", 1)], "N--BB")),
        ("reject_node_on_occupied_direction",
         "A node supplied for an occupied direction (an extraneous branch).",
         doc("HHHH", [("2", 1)], "NBBBB")),
        ("reject_exact_at_nonterminal", "`E` at a prefix with residues left.",
         doc("HHHH", [("2", 1)], "E")),
        ("reject_exact_child_at_nonterminal", "`E` one level too early inside a split.",
         doc("HHHH", [("2", 1)], "NE-EE")),
        ("reject_split_at_terminal", "`N` at a complete walk.",
         doc("HP", [("2", 0)], "N----")),
        ("reject_truncated_tree", "The input ends inside the tree.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1\ntree NB-"),
        ("reject_trailing_tree_bytes", "Extra bytes after a complete tree.",
         doc("HHHH", [("2", 1)], "XB")),
        ("reject_unknown_tree_byte", "A byte outside `N E X B -`.",
         doc("HHHH", [("2", 1)], "Q")),
        ("reject_lowercase_tree_byte", "Tree bytes are case sensitive.",
         doc("HHHH", [("2", 1)], "x")),
        ("reject_bad_magic", "A different format version.",
         good.replace("KERNEL 1", "KERNEL 2")),
        ("reject_missing_magic", "No magic line.",
         good.replace(HEADER, "")),
        ("reject_bad_sequence_symbol", "A symbol other than `H`/`P` in the sequence.",
         doc("HHAH", [("2", 1)], "X")),
        ("reject_short_sequence", "A one-residue sequence.",
         doc("H", [("2", 1)], "X")),
        ("reject_sequence_substitution",
         "A valid tree for `HHHH` presented with a longer sequence.",
         doc("HHHHH", [("2", 1)], hhhh_full)),
        ("reject_q_below_one", "A Boltzmann factor below one.",
         doc("HHHH", [("1/2", 1)], "X")),
        ("reject_q_zero", "A zero Boltzmann factor.",
         doc("HHHH", [("0", 1)], "X")),
        ("reject_q_not_lowest_terms", "A rational not in lowest terms.",
         doc("HHHH", [("4/2", 1)], "X")),
        ("reject_q_denominator_one", "An integer written as a fraction.",
         doc("HHHH", [("3/1", 1)], "X")),
        ("reject_q_zero_denominator", "A zero denominator.",
         doc("HHHH", [("1/0", 1)], "X")),
        ("reject_q_leading_zero", "A leading zero.",
         doc("HHHH", [("02", 1)], "X")),
        ("reject_q_sign", "An explicit sign.",
         doc("HHHH", [("+2", 1)], "X")),
        ("reject_q_decimal", "A decimal fraction.",
         doc("HHHH", [("2.0", 1)], "X")),
        ("reject_q_negative", "A negative rational.",
         doc("HHHH", [("-2", 1)], "X")),
        ("reject_threshold_leading_zero", "A threshold with a leading zero.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 01\ntree X\nend\n"),
        ("reject_threshold_missing", "A query without a threshold.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2\ntree X\nend\n"),
        ("reject_threshold_negative", "A negative threshold.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 -1\ntree X\nend\n"),
        ("reject_query_trailing_space", "A trailing space on a query line.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1 \ntree X\nend\n"),
        ("reject_no_query", "No query line.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\ntree X\nend\n"),
        ("reject_duplicate_sequence_line", "Two sequence lines.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nsequence HHHH\nquery 2 1\ntree X\nend\n"),
        ("reject_query_after_tree", "A query line after the tree.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1\ntree X\nquery 2 1\nend\n"),
        ("reject_missing_tree_line", "No tree line.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1\nend\n"),
        ("reject_missing_end", "No `end` line.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1\ntree X\n"),
        ("reject_missing_final_newline", "No newline after `end`.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1\ntree X\nend"),
        ("reject_trailing_blank_line", "A blank line after `end`.",
         "FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1\ntree X\nend\n\n"),
        ("reject_blank_line_in_header", "A blank line inside the header.",
         "FOLDBOUNDS-KERNEL 1\n\nsequence HHHH\nquery 2 1\ntree X\nend\n"),
        ("reject_crlf", "Windows line endings.",
         good.replace("\n", "\r\n")),
        ("reject_tree_separator", "A space inside the tree.",
         doc("HHHH", [("2", 1)], "NB -BB")),
        ("reject_empty_input", "The empty input.", ""),
    ]
    for name, comment, text in cases:
        out.append(reject_case(name, comment, text))
    out.append(reject_bytes_case(
        "reject_invalid_utf8",
        "A byte that is not valid UTF-8 inside the sequence line.",
        b"FOLDBOUNDS-KERNEL 1\nsequence HH\xffH\nquery 2 1\ntree X\nend\n"))

    h8_report = run(doc(h8, [("2", 2), ("3/2", 1)], h8_hybrid).encode())
    first = h8_report["results"][0]
    exact8 = run(doc(h8, [("2", 2)], "X").encode())["results"][0]
    out.append(
        f'''/-! ### From one Boolean check to a theorem about the Lean HP model -/

/-- All-hydrophobic sequence of four residues. -/
def allH4 : Fin 4 → Bool := fun _ => true

/-- All-hydrophobic sequence of eight residues. -/
def allH8 : Fin 8 → Bool := fun _ => true

/-- **Exact probability from the string-level checker.** The kernel
evaluates `verify` on a one-leaf certificate; `certifiesExact_sound` turns
the accepted report into the exact equilibrium probability `P(k ≥ 1) = 4/11`
of the Lean HP model for `HHHH` at `q = 2`. -/
theorem hhhh_probability_exact :
    Executable.massQ (weightQ 2 2 allH4) (contactEvent 2 allH4 1) /
      Executable.massQ (weightQ 2 2 allH4) Finset.univ = 4 / 11 :=
  certifiesExact_sound allH4
    ({lean_string(doc("HHHH", [("2", 1)], "X"))}).toUTF8
    2 1 (4 / 11) (by decide +kernel)

/-- **Exact probability for eight residues** (543 conformations enumerated by
the checker inside the kernel). -/
theorem h8_probability_exact :
    Executable.massQ (weightQ 6 2 allH8) (contactEvent 6 allH8 2) /
      Executable.massQ (weightQ 6 2 allH8) Finset.univ = {lean_rat(exact8["pLo"])} :=
  certifiesExact_sound allH8
    ({lean_string(doc(h8, [("2", 2)], "X"))}).toUTF8
    2 2 ({lean_rat(exact8["pLo"])}) (by decide +kernel)

/-- **A certified enclosure from a partially bounded certificate.** The
tree below enumerates 152 of the 543 conformations and bounds 17 subtrees
without enumerating them; `certifiesInterval_sound` gives a kernel-checked
enclosure of the same probability. -/
theorem h8_probability_interval :
    Executable.massQ (weightQ 6 2 allH8) (contactEvent 6 allH8 2) /
      Executable.massQ (weightQ 6 2 allH8) Finset.univ ∈
        Set.Icc ({lean_rat(first["pLo"])} : ℚ) ({lean_rat(first["pHi"])}) :=
  certifiesInterval_sound allH8
    ({lean_string(doc(h8, [("2", 2), ("3/2", 1)], h8_hybrid))}).toUTF8
    2 2 ({lean_rat(first["pLo"])}) ({lean_rat(first["pHi"])}) (by decide +kernel)

/-- The exact value indeed lies in the certified enclosure. -/
example : ({lean_rat(first["pLo"])} : ℚ) ≤ {lean_rat(exact8["pLo"])} ∧
    ({lean_rat(exact8["pLo"])} : ℚ) ≤ {lean_rat(first["pHi"])} := by norm_num

end FoldBounds.FrontierKernelTests
'''
    )
    sys.stdout.write("\n".join(out))


if __name__ == "__main__":
    main()
