"""Strict replay of authenticated GA/GB NMR structural JSON, with only tiny float jitter.

All JSON keys, array lengths/order, source hashes, strings, model counts,
integer values, booleans, classification labels, and null values must match
EXACTLY. Only real JSON floating-point values are compared to 1e-12 relative
and absolute tolerance to accommodate cross-platform math library rounding.
Neither source files nor the frozen expected measurements are rewritten.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

REL_TOL = 1e-12
ABS_TOL = 1e-12


def _reject_constant(value: str):
    raise ValueError(f"Forbidden nonfinite JSON constant: {value}")


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=_reject_constant)


def compare_json(actual, expected):
    """Return (exactly-described errors, tiny-ULP float difference count, max delta)."""
    errors: list[str] = []
    rounded = 0
    max_delta = 0.0

    def visit(a, e, location: str):
        nonlocal rounded, max_delta
        if type(a) is not type(e):
            errors.append(f"{location}: JSON type differs ({type(a).__name__} vs {type(e).__name__})")
            return
        if isinstance(e, dict):
            missing = sorted(set(e) - set(a))
            extra = sorted(set(a) - set(e))
            if missing or extra:
                errors.append(f"{location}: keys differ, missing={missing}, unexpected={extra}")
            for k in sorted(set(a) & set(e)):
                visit(a[k], e[k], f"{location}.{k}")
        elif isinstance(e, list):
            if len(a) != len(e):
                errors.append(f"{location}: array length {len(a)} differs from frozen {len(e)}")
            for i, (x, y) in enumerate(zip(a, e)):
                visit(x, y, f"{location}[{i}]")
        elif isinstance(e, float):
            if not math.isfinite(a) or not math.isfinite(e):
                errors.append(f"{location}: non-finite floating point value")
            else:
                delta = abs(a - e)
                max_delta = max(max_delta, delta)
                if not math.isclose(a, e, rel_tol=REL_TOL, abs_tol=ABS_TOL):
                    errors.append(f"{location}: numeric mismatch: actual={a!r}, frozen={e!r}, abs_delta={delta:.17g}")
                elif a != e:
                    rounded += 1
        elif a != e:
            errors.append(f"{location}: exact mismatch: actual={a!r}, frozen={e!r}")

    visit(actual, expected, "$")
    return errors, rounded, max_delta


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("actual", type=Path)
    parser.add_argument("frozen", type=Path)
    args = parser.parse_args()
    actual = load_json(args.actual)
    frozen = load_json(args.frozen)
    errors, rounded, max_delta = compare_json(actual, frozen)
    if errors:
        for issue in errors[:30]:
            print("FAIL:", issue)
        if len(errors) > 30:
            print(f"... and {len(errors)-30} more mismatches")
        parser.exit(1, f"FAIL: {len(errors)} structural/scientific JSON discrepancies; frozen reference NOT accepted\n")
    print(
        f"PASS: frozen JSON replay; exact structure, labels/counts/provenance;"
        f" {rounded} sub-tolerance float differences; max absolute difference {max_delta:.17g}"
    )


if __name__ == "__main__":
    main()
