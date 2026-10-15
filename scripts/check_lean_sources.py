#!/usr/bin/env python3
"""Fail-closed Lean source hygiene and declaration-driven axiom-output audit."""
from pathlib import Path
import re
import sys

source_files = sorted(Path("FoldBounds").glob("*.lean"))
if len(source_files) < 10:
    raise SystemExit(f"Missing expected Lean modules: found only {len(source_files)}")

for path in source_files:
    source = path.read_text(encoding="utf-8")
    if re.search(r"\b(?:sorry|admit)\b|\baxiom\s+\w+|@\[\s*implemented_by\b", source):
        raise SystemExit(f"Forbidden proof hole, custom axiom, or implementation override: {path}")
print(f"Source audit: {len(source_files)} Lean modules, no forbidden proof bypasses")

audits = re.findall(
    r"^#print axioms (FoldBounds\.[A-Za-z0-9_.]+)\s*$",
    Path("FoldBounds/AxiomAudit.lean").read_text(encoding="utf-8"),
    re.MULTILINE,
)
if len(audits) < 34 or len(audits) != len(set(audits)):
    raise SystemExit("Audit must preserve 34 original declarations and have no duplicates")
if "FoldBounds.Executable.check_sound" not in audits:
    raise SystemExit("Executable.check_sound must be included in the axiom audit")

if len(sys.argv) > 1:
    log = Path(sys.argv[1]).read_text(encoding="utf-8")
    reports = {}
    for line in log.splitlines():
        match = re.match(
            r"^'(FoldBounds\.[A-Za-z0-9_.]+)' "
            r"(does not depend on any axioms|depends on axioms:?\s*.*)$",
            line.strip(),
        )
        if match:
            name, result = match.groups()
            if name in reports:
                raise SystemExit(f"Duplicate axiom output: {name}")
            reports[name] = result
            if result.startswith("depends on axioms"):
                tail = result.split("depends on axioms", 1)[1]
                axioms = set(re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", tail))
                unexpected = axioms - {"Classical.choice", "Quot.sound", "propext"}
                if unexpected:
                    raise SystemExit(f"Unexpected axioms for {name}: {sorted(unexpected)}")
    if set(reports) != set(audits):
        missing = sorted(set(audits) - set(reports))
        extras = sorted(set(reports) - set(audits))
        raise SystemExit(f"Axiom-audit mismatch. Missing: {missing}; extras: {extras}")
    print(f"Axiom audit: {len(audits)} declarations, standard axioms only")
