"""Verify the archived HP benchmark CSVs, and optionally a fresh regeneration of them.

usage:
  python3 experiments/coupled/verify_hp_archive.py ARCHIVE_DIR [REGEN_DIR]

ARCHIVE_DIR holds the six CSVs retrieved from the archive ref listed in
experiments/coupled/results/DATA_MANIFEST.json. Each file's SHA-256, size and row count must
match the manifest.

If REGEN_DIR is given, it holds the same six files regenerated with the manifest's `generator`
commands. Every row must then agree with the archive in every column except the
non-deterministic wall-clock column(s) listed in the manifest (`t_sec`). The script exits
non-zero on any mismatch.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

MANIFEST = Path(__file__).resolve().parent / "results" / "DATA_MANIFEST.json"


def main() -> int:
    if len(sys.argv) not in (2, 3):
        print(__doc__)
        return 2
    archive = Path(sys.argv[1])
    regen = Path(sys.argv[2]) if len(sys.argv) == 3 else None
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    failures = 0
    for entry in manifest["files"]:
        name = entry["path"].rsplit("/", 1)[-1]
        data = (archive / name).read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        rows = list(csv.reader(data.decode("utf-8").splitlines()))
        ok = digest == entry["sha256"] and len(data) == entry["bytes"] and len(rows) - 1 == entry["data_rows"]
        failures += not ok
        print(f"{'OK ' if ok else 'BAD'} archive  {name}  sha256={digest[:16]}  rows={len(rows) - 1}")
        if regen is None:
            continue
        fresh = list(csv.reader((regen / name).read_text(encoding="utf-8").splitlines()))
        header = rows[0]
        if fresh[0] != header or len(fresh) != len(rows):
            print(f"BAD regen    {name}  header or row count differs ({len(fresh) - 1} vs {len(rows) - 1})")
            failures += 1
            continue
        skip = {header.index(c) for c in entry["nondeterministic_columns"] if c in header}
        diff = sum(
            1
            for a, b in zip(rows[1:], fresh[1:])
            if [x for j, x in enumerate(a) if j not in skip] != [x for j, x in enumerate(b) if j not in skip]
        )
        failures += diff > 0
        print(f"{'OK ' if diff == 0 else 'BAD'} regen    {name}  differing rows (ignoring {sorted(entry['nondeterministic_columns'])}): {diff}")
    print("VERIFIED" if failures == 0 else f"FAILED ({failures} problem(s))")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
