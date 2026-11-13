#!/usr/bin/env python3
"""Fail-closed Lean source hygiene and declaration-driven axiom-output audit.

Usage:
    python3 scripts/check_lean_sources.py            # source audit only
    python3 scripts/check_lean_sources.py AXIOM_LOG  # source audit + axiom-output audit

AXIOM_LOG is the output of `lake env lean FoldBounds/AxiomAudit.lean`.

Hardening (findings F3e and F3f of the 2026-11 assurance audit):

* F3e. Lean wraps long `#print axioms` reports over several lines. The previous parser read
  only the first line of each report, so an axiom printed on a continuation line was never
  examined. Reports are now parsed across lines up to the closing bracket. The following all
  FAIL: a report whose bracket never closes (truncated output), trailing text after the
  bracket, an unconsumed report-like line, and a report in an unrecognised format.
* F3f. A user axiom declared as `axiom «Classical.choice» : False` is a single name component
  that contains a dot, and Lean prints it as `«Classical.choice»`. The previous parser
  extracted identifier-like substrings and read that as the trusted `Classical.choice`. The
  previous source regex `\\baxiom\\s+\\w+` also did not match `axiom «`. Now every axiom
  entry, after splitting on commas, must EXACTLY equal one of the authorised names. Anything
  else FAILS: guillemets, private names such as `_private.M.0.propext`, primed or
  macro-scoped variants, other namespaces. The source scan also matches `axiom` followed by
  `«`.

Known coverage limits. These are documented and deliberately NOT changed here; each needs its
own separately reviewed strengthening:
  1. The source scan covers only `FoldBounds/*.lean` (non-recursive), not subdirectories or
     `scripts/*.lean`.
  2. `native_decide` is not forbidden in source. It is caught only when a declaration that
     uses it is listed in AxiomAudit.lean, via `Lean.ofReduceBool`.
  3. Only declarations listed in `FoldBounds/AxiomAudit.lean` are audited.
  4. The exact-name check relies on Lean's name printer rendering non-standard names
     distinctly. This was verified empirically for guillemet and private names on Lean 4.28.0.
     A Lean-side `collectAxioms` check with `Name` equality would remove this dependency.
"""
from pathlib import Path
import re
import sys

AUTHORISED_AXIOMS = frozenset({"propext", "Classical.choice", "Quot.sound"})

# `axiom` followed, after optional whitespace, by an identifier start or a guillemet. The
# lookahead keeps the plural "axioms" in prose from matching. Every string matched by the
# previous pattern `\baxiom\s+\w+` is still matched.
FORBIDDEN_SOURCE = re.compile(
    r"\b(?:sorry|admit)\b|\baxiom(?=[\s«])\s*[\w«]|@\[\s*implemented_by\b"
)

REPORT_HEADER = re.compile(
    r"^'(?P<name>FoldBounds\.[A-Za-z0-9_.]+)' "
    r"(?P<kind>does not depend on any axioms|depends on axioms:)(?P<rest>.*)$"
)


def fail(message: str) -> None:
    raise SystemExit(message)


def parse_axiom_reports(log: str) -> dict[str, frozenset[str]]:
    """Parse every `#print axioms` report in LOG, failing closed on anything unexpected."""
    lines = log.splitlines()
    reports: dict[str, frozenset[str]] = {}
    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\r")
        match = REPORT_HEADER.match(line)
        if match is None:
            # A report-like line that the strict grammar does not accept must not be skipped.
            if "depends on axioms" in line or "does not depend on any axioms" in line:
                fail(f"Unparseable axiom report (line {i + 1}): {line!r}")
            if line.lstrip().startswith("'FoldBounds."):
                fail(f"Malformed axiom report header (line {i + 1}): {line!r}")
            i += 1
            continue
        name, kind, rest = match.group("name"), match.group("kind"), match.group("rest")
        if name in reports:
            fail(f"Duplicate axiom output: {name}")
        if kind == "does not depend on any axioms":
            if rest.strip():
                fail(f"Trailing text after axiom report for {name}: {rest!r}")
            reports[name] = frozenset()
            i += 1
            continue
        rest = rest.strip()
        if not rest.startswith("["):
            fail(f"Axiom list for {name} does not start with '[': {rest!r}")
        body = rest[1:]
        # F3e: collect continuation lines until the closing bracket. A continuation line is
        # indented; a non-indented line or the end of the log before ']' means truncation.
        while "]" not in body:
            i += 1
            if i >= len(lines):
                fail(f"Truncated axiom report for {name}: no closing ']'")
            continuation = lines[i].rstrip("\r")
            if not continuation[:1].isspace():
                fail(f"Truncated axiom report for {name}: line {i + 1} is not a continuation")
            body += "\n" + continuation
        inside, _, after = body.partition("]")
        if after.strip():
            fail(f"Trailing text after axiom report for {name}: {after!r}")
        if "[" in inside:
            fail(f"Nested '[' in axiom report for {name}")
        entries = [entry.strip() for entry in inside.split(",")]
        if any(entry == "" for entry in entries):
            fail(f"Empty axiom entry in report for {name}: {inside!r}")
        # F3f: exact comparison with the authorised names, never substring extraction.
        unexpected = sorted(set(entries) - AUTHORISED_AXIOMS)
        if unexpected:
            fail(f"Unexpected axioms for {name}: {unexpected}")
        reports[name] = frozenset(entries)
        i += 1
    return reports


def main() -> None:
    source_files = sorted(Path("FoldBounds").glob("*.lean"))
    if len(source_files) < 10:
        fail(f"Missing expected Lean modules: found only {len(source_files)}")
    for path in source_files:
        if FORBIDDEN_SOURCE.search(path.read_text(encoding="utf-8")):
            fail(f"Forbidden proof hole, custom axiom, or implementation override: {path}")
    print(f"Source audit: {len(source_files)} Lean modules, no forbidden proof bypasses")

    audits = re.findall(
        r"^#print axioms (FoldBounds\.[A-Za-z0-9_.]+)\s*$",
        Path("FoldBounds/AxiomAudit.lean").read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    if len(audits) < 34 or len(audits) != len(set(audits)):
        fail("Audit must preserve 34 original declarations and have no duplicates")
    if "FoldBounds.Executable.check_sound" not in audits:
        fail("Executable.check_sound must be included in the axiom audit")

    if len(sys.argv) > 1:
        reports = parse_axiom_reports(Path(sys.argv[1]).read_text(encoding="utf-8"))
        if set(reports) != set(audits):
            missing = sorted(set(audits) - set(reports))
            extras = sorted(set(reports) - set(audits))
            fail(f"Axiom-audit mismatch. Missing: {missing}; extras: {extras}")
        print(f"Axiom audit: {len(audits)} declarations, standard axioms only")


if __name__ == "__main__":
    main()
