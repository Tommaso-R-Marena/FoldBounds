#!/usr/bin/env python3
"""Fail-closed source and axiom audit for the Coupled Ensemble Certificate (CEC) Lean files.

Usage:
    python3 scripts/check_coupled_lean.py                 # source audit only
    python3 scripts/check_coupled_lean.py AXIOM_LOG       # source audit + axiom-output audit

AXIOM_LOG is the output of `lake env lean FoldBounds/AxiomAudit.lean`.

Source audit: every FoldBounds/Coupled*.lean file must be free of proof holes and trust
escapes.  The match is purely lexical and also applies to comments, so it fails closed.
Axiom audit: every audited FoldBounds.Coupled* declaration must be reported, and only the
standard axioms propext, Classical.choice, Quot.sound are allowed; in particular any
dependence on Lean.ofReduceBool (native evaluation) or Lean.trustCompiler fails the check.
"""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent.parent
LIBRARY = ["CoupledRatio", "CoupledSeparator", "CoupledBinding"]
TESTS = "CoupledTests"
STANDARD = {"propext", "Classical.choice", "Quot.sound"}
REQUIRED_AUDITS = {
    "FoldBounds.CoupledRatio.forall_le_Q_iff",
    "FoldBounds.CoupledRatio.forall_Q_le_iff",
    "FoldBounds.CoupledRatio.checkLower_sound",
    "FoldBounds.CoupledRatio.checkUpper_sound",
    "FoldBounds.CoupledRatio.checkEnclosure_sound",
    "FoldBounds.CoupledRatio.le_Q_of_lowerDualCert",
    "FoldBounds.CoupledSeparator.local_enclosure",
    "FoldBounds.CoupledSeparator.sharp_closed",
    "FoldBounds.CoupledSeparator.checkSeparator_sound",
    "FoldBounds.CoupledBinding.selectivity_mem_Icc",
    "FoldBounds.CoupledBinding.nonidentifiable",
}
FORBIDDEN = [
    (r"\bsorry\b", "sorry"),
    (r"\badmit\b", "admit"),
    (r"\baxiom\b", "axiom declaration (the singular word 'axiom' is reserved)"),
    (r"\bnative_decide\b", "native_decide"),
    (r"\bnative\b", "native evaluation (e.g. decide +native)"),
    (r"\bimplemented_by\b", "implemented_by"),
    (r"\bextern\b", "extern"),
    (r"\bunsafe\b", "unsafe"),
    (r"\bopaque\b", "opaque"),
    (r"ofReduceBool", "Lean.ofReduceBool"),
    (r"trustCompiler", "Lean.trustCompiler"),
    (r"skipKernelTC", "kernel type-checking bypass"),
]


def fail(msg: str) -> None:
    raise SystemExit(f"check_coupled_lean: FAIL: {msg}")


def source_audit() -> None:
    files = sorted((ROOT / "FoldBounds").glob("Coupled*.lean"))
    names = {p.stem for p in files}
    for required in LIBRARY + [TESTS]:
        if required not in names:
            fail(f"missing FoldBounds/{required}.lean")
    for path in files:
        text = path.read_text(encoding="utf-8")
        for pattern, label in FORBIDDEN:
            match = re.search(pattern, text)
            if match:
                line = text.count("\n", 0, match.start()) + 1
                fail(f"{path.relative_to(ROOT)}:{line}: forbidden token ({label})")
    root = (ROOT / "FoldBounds.lean").read_text(encoding="utf-8")
    for module in LIBRARY:
        if not re.search(rf"^public import FoldBounds\.{module}\s*$", root, re.MULTILINE):
            fail(f"FoldBounds.lean must publicly import FoldBounds.{module}")
    if re.search(rf"^\s*(public\s+)?import FoldBounds\.{TESTS}\b", root, re.MULTILINE):
        fail("the test module must not be imported by the library root")
    tests = (ROOT / "FoldBounds" / f"{TESTS}.lean").read_text(encoding="utf-8")
    if len(re.findall(r"\bdecide\b", tests)) < 10:
        fail("CoupledTests.lean must contain kernel-evaluated decide checks")
    print(f"Coupled source audit: {len(files)} files, no forbidden tokens")


def audited_declarations() -> list[str]:
    audit = (ROOT / "FoldBounds" / "AxiomAudit.lean").read_text(encoding="utf-8")
    names = re.findall(
        r"^#print axioms (FoldBounds\.Coupled[A-Za-z0-9_]*\.[A-Za-z0-9_.]+)\s*$",
        audit,
        re.MULTILINE,
    )
    if len(names) != len(set(names)):
        fail("duplicate Coupled entries in FoldBounds/AxiomAudit.lean")
    missing = REQUIRED_AUDITS - set(names)
    if missing:
        fail(f"required declarations missing from the axiom audit: {sorted(missing)}")
    if len(names) < 60:
        fail(f"expected at least 60 audited Coupled declarations, found {len(names)}")
    return names


def axiom_audit(log_path: str, names: list[str]) -> None:
    log = Path(log_path).read_text(encoding="utf-8")
    if re.search(r"\berror\b", log):
        fail("the axiom-audit log contains an error")
    # Reports for long names are wrapped over several lines by Lean's pretty printer, so the
    # bracketed axiom list is matched across newlines.  Every report header must parse.
    report = re.compile(
        r"^'(FoldBounds\.Coupled[A-Za-z0-9_]*\.[A-Za-z0-9_.]+)' "
        r"(?:does not depend on any axioms|depends on axioms:\s*\[([^\]]*)\])",
        re.MULTILINE,
    )
    headers = re.findall(r"^'FoldBounds\.Coupled", log, re.MULTILINE)
    reports: dict[str, set[str]] = {}
    for match in report.finditer(log):
        name = match.group(1)
        if name in reports:
            fail(f"duplicate axiom report for {name}")
        reports[name] = set(re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", match.group(2) or ""))
    if len(headers) != sum(1 for _ in report.finditer(log)):
        fail("an axiom report for a Coupled declaration could not be parsed completely")
    for name in names:
        if name not in reports:
            fail(f"no axiom report for {name}")
        axioms = reports[name]
        if "Lean.ofReduceBool" in axioms or "Lean.trustCompiler" in axioms:
            fail(f"{name} depends on native evaluation: {sorted(axioms)}")
        if axioms - STANDARD:
            fail(f"{name} depends on non-standard axioms: {sorted(axioms - STANDARD)}")
    extras = set(reports) - set(names)
    if extras:
        fail(f"unexpected Coupled reports not listed in the audit: {sorted(extras)}")
    print(f"Coupled axiom audit: {len(names)} declarations, standard axioms only")


def main() -> None:
    source_audit()
    names = audited_declarations()
    if len(sys.argv) > 1:
        axiom_audit(sys.argv[1], names)


if __name__ == "__main__":
    main()
