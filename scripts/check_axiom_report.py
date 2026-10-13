"""Validate Lean's printed axiom dependencies. Kernel build remains required."""
import re
import sys
from pathlib import Path

text=Path(sys.argv[1]).read_text(encoding='utf-8')
expected=len(re.findall(r'(?m)^#print axioms FoldBounds\.',
                        Path('FoldBounds/AxiomAudit.lean').read_text()))
entries=re.findall(r'(?m)^FoldBounds\.\S+ depends on axioms:\s*\[([^\]]*)\]',text)
if len(entries)<expected:
    raise SystemExit(f'Axiom report incomplete: {len(entries)} entries, expected {expected}')
allowed={'propext','Classical.choice','Quot.sound'}
for entry in entries:
    axioms={item.strip() for item in entry.split(',') if item.strip()}
    if axioms-allowed:
        raise SystemExit(f'Unexpected axioms: {axioms-allowed}')
print(f'Checked {len(entries)} Lean axiom reports.')
