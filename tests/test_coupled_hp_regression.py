"""Regression fixture for the archived HP coupled-certificate benchmark (PR #96, owner decision 2).

The raw benchmark CSVs live in versioned git storage (see
experiments/coupled/results/DATA_MANIFEST.json). This test keeps a representative,
fast-to-regenerate slice in the repository. It is the 1,290 archived `oracle_n10` rows for the
first sequence with fewer than 100 regions, and they must be reproduced exactly (every column
except the wall-clock `t_sec`) by the recorded generator run with a smaller cap.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "coupled_hp" / "oracle_n10_seq0_regions_lt100.csv"
MANIFEST = ROOT / "experiments" / "coupled" / "results" / "DATA_MANIFEST.json"


def _rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_oracle_slice_regenerates_exactly(tmp_path: Path) -> None:
    out = tmp_path / "oracle_n10_small.csv"
    subprocess.run(
        [sys.executable, "experiments/coupled/run_oracle.py", "10", str(out), "100", "1"],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": "src"},
        check=True,
        capture_output=True,
        text=True,
    )
    fixture = _rows(FIXTURE)
    columns = list(fixture[0].keys())
    assert "t_sec" not in columns
    regenerated = [row for row in _rows(out) if int(row["regions"]) < 100]
    key = lambda row: tuple(row[c] for c in columns)  # noqa: E731
    assert len(fixture) == 1290
    assert Counter(map(key, regenerated)) == Counter(map(key, fixture))


def test_fixture_rows_are_sound_and_bounded_by_trivial() -> None:
    # The archived negative result in miniature: every interval is sound, and the coupled
    # width never exceeds the trivial ||V||_inf width.
    for row in _rows(FIXTURE):
        assert row["sound_triv"] == row["sound_sep"] == row["sound_coup"] == "1"
        assert float(row["coup_w"]) <= float(row["triv_w"]) + 1e-9


def test_manifest_lists_every_archived_file_with_checksums() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    names = {entry["path"].rsplit("/", 1)[-1] for entry in manifest["files"]}
    assert names == {
        "oracle_n10.csv", "oracle_n12.csv", "oracle_n14.csv",
        "deep_n12.csv", "large_n30.csv", "large_n50.csv",
    }
    for entry in manifest["files"]:
        assert len(entry["sha256"]) == 64 and int(entry["bytes"]) > 0 and entry["generator"]
    assert manifest["archive"]["commit"].startswith("5303e02")
    assert sum(e["data_rows"] for e in manifest["files"] if "oracle_" in e["path"]) == 71859
