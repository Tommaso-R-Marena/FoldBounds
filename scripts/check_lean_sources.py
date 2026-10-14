#!/usr/bin/env python3
"""Static hygiene plus strict 34-theorem Lean axiom-output gate."""
from pathlib import Path
import re
import sys

files=sorted(Path("FoldBounds").glob("*.lean"))
if len(files)!=10:
    raise SystemExit(f"Expected 10 FoldBounds Lean modules, got {len(files)}")
for file in files:
    data=file.read_text()
    if re.search(r"\b(?:sorry|admit)\b|\baxiom\s+\w+|@\[\s*implemented_by\b", data):
        raise SystemExit("Forbidden proof hole/axiom/implemented_by in "+str(file))
print("Source audit: 10 Lean modules, no sorry/admit/new axiom/implemented_by")
if len(sys.argv)>1:
    log=Path(sys.argv[1]).read_text()
    statements=[
        line for line in log.splitlines()
        if "depends on axioms" in line
    ]
    if len(statements)!=34:
        raise SystemExit(f"Expected 34 axiom audit statements, found {len(statements)}")
    for line in statements:
        if not line.startswith("'FoldBounds.") or "depends on axioms" not in line:
            raise SystemExit("Unexpected theorem/axiom output: "+line)
        tail=line.split("depends on axioms",1)[1]
        names=re.findall(r"\b(?:Classical\.choice|Quot\.sound|propext|[A-Za-z_][A-Za-z0-9_.]*)\b",tail)
        unexpected=set(names)-{"Classical.choice","Quot.sound","propext"}
        if unexpected:
            raise SystemExit(f"Additional axioms detected: {sorted(unexpected)}")
    print("Axiom audit: 34 declarations; only standard Lean axioms")
