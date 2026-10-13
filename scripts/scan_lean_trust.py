"""Defense-in-depth source scan, not a substitute for Lean's kernel."""
import re
from pathlib import Path

paths=[Path('FoldBounds.lean'), *Path('FoldBounds').glob('*.lean')]
if not paths or not all(path.exists() for path in paths):
    raise SystemExit('Lean sources not imported')
for path in paths:
    data=path.read_text(encoding='utf-8')
    for bad in (r'(?m)^\s*(?:private\s+)?axiom\s+',
                r'(?m)^\s*(?:sorry|admit)\b',
                r'@\[implemented_by\s+'):
        if re.search(bad,data):
            raise SystemExit(f'Forbidden source construct in {path}: {bad}')
print(f'Source scan completed on {len(paths)} Lean source files.')
