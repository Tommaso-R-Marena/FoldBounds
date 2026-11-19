#!/usr/bin/env python3
"""Long-budget MINWALK on hard S1-* instances; dumps PATH for FoldBounds certs."""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BENCH = json.loads((Path(__file__).parent / "benchmarks.json").read_text())
BIN = ROOT / "third_party" / "minwalk" / "dist" / "hpfold"
OUTDIR = Path("/tmp/minwalk_long")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ids = ["S1-6", "S1-7"]
    for iid in ids:
        inst = next(x for x in BENCH["ground_state_2d"] if x["id"] == iid)
        inp = OUTDIR / f"{iid}.in"
        outp = OUTDIR / f"{iid}.out"
        log = OUTDIR / f"{iid}.log"
        inp.write_text(inst["sequence"] + "\n")
        t0 = time.time()
        print(f"START {iid} n={inst['n']} contacts={inst['contacts']}", flush=True)
        try:
            proc = subprocess.run(
                [str(BIN), str(inp), str(outp)],
                capture_output=True,
                text=True,
                timeout=7200,
                check=False,
            )
            wall = time.time() - t0
            text = (outp.read_text() if outp.exists() else "") + "\n" + (proc.stdout or "")
            log.write_text(f"wall={wall}\nrc={proc.returncode}\n{text}")
            print(f"DONE {iid} wall={wall:.1f}s\n{text[:800]}", flush=True)
        except subprocess.TimeoutExpired:
            print(f"TIMEOUT {iid} after 7200s", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
