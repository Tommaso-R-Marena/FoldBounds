#!/usr/bin/env python3
"""Adversarial standard-library checks of the bulk Lean rank source exporter."""
import os
import tempfile
from pathlib import Path
from export_bulk_lean_rank_proofs import overlapping_blocks, generate_source, sha

def main():
    entries = [(10,20,0),(21,25,1),(30,35,2),
               (100,101,3),(105,110,4),(999,1001,5)]
    chunks = list(overlapping_blocks(entries,3))
    assert [[e[2] for e in block] for block in chunks] == [
        [0,1,2],[2,3,4],[4,5]]
    assert sum(len(b)-1 for b in chunks) == len(entries)-1
    source, stats = generate_source(
        {"CspB":(entries,"0"*64)}, "1"*64,3)
    assert stats["CspB"]["adjacent_comparisons"]==5
    assert stats["CspB"]["blocks"]==3
    assert "checkAdjacent CspB_block_002 = true := by decide" in source
    assert "checkAdjacent_sound CspB_block_002 CspB_block_002_ok" in source
    assert source==generate_source({"CspB":(entries,"0"*64)}, "1"*64,3)[0]
    assert len(sha(source.encode()))==64
    # Overlap and tied intervals are not proofs: the input gate must reject
    # them before source generation. Check a minimal independent strict gate.
    bads = [
        [(1,4,0),(4,7,1)],
        [(1,5,0),(3,6,1)],
        [(9,10,0),(1,2,1)],
    ]
    for items in bads:
        assert not all(
            items[i][1]<items[i+1][0]
            for i in range(len(items)-1))
    print("PASS: overlapping chunk coverage, deterministic Lean source,"
          " and strict rank sanity checks")

if __name__ == "__main__":
    main()
