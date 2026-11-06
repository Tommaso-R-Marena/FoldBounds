"""Reject import-order mistakes in integration source before starting a costly Lean build."""
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PATHS=(ROOT/"FoldBounds/IntegrationImports.lean",ROOT/"FoldBounds/AxiomAudit.lean")
for path in PATHS:
    lines=path.read_text(encoding="utf-8").splitlines()
    stripped=[line.strip() for line in lines if line.strip()]
    if not stripped or not stripped[0].startswith("import "):
        raise SystemExit(f"ERROR: {path}: first nonempty line must be a Lean import")
    saw_other=False
    count=0
    for line in stripped:
        if line.startswith("import "):
            if saw_other:
                raise SystemExit(f"ERROR: {path}: import occurs after a non-import command")
            count+=1
        else:
            saw_other=True
    if count==0:
        raise SystemExit(f"ERROR: {path}: no imports")
    print(f"PASS: {path.name} has {count} leading imports and no late imports")
