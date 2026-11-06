#!/usr/bin/env python3
"""Offline fallback test runner (used when pytest is not installed).

Supports only pytest.raises and pytest.mark.parametrize, which is all the
FoldBounds test-suite uses.  Prefer `python -m pytest` when available.
Run from the repository root:  python3 scripts/run_tests_offline.py
"""
import sys, inspect, importlib.util, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "offline_pytest")); sys.path.insert(0, "src")
ok=fail=skip=0
for f in sorted(pathlib.Path("tests").glob("test_*.py")):
    spec=importlib.util.spec_from_file_location(f.stem, f)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    for n,fn in inspect.getmembers(m, inspect.isfunction):
        if not n.startswith("test_") or fn.__module__!=m.__name__: continue
        cases = getattr(fn, "_params", [{}])
        params = inspect.signature(fn).parameters
        for c in cases:
            if set(params) - set(c): skip+=1; continue
            try: fn(**c); ok+=1
            except Exception as e:
                fail+=1; print("FAIL", f.name, n, c, type(e).__name__, str(e)[:100])
print(f"passed={ok} failed={fail} skipped(fixtures)={skip}")
