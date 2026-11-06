# Version-2 frontier certificate format (`FOLDBOUNDS-FRONTIER 2`)

Trusted core: `FoldBounds.FrontierNext.check` (trie form) and
`FoldBounds.FrontierNextFlat.checkFlat2` (flat form).  Soundness:
`FrontierNext.check_sound`, `FrontierNextFlat.checkFlat2_sound`.
Text parsing (`FoldBounds/FrontierNextCLI.lean`) and the `foldcheck2`
executable's I/O are an **unverified boundary**; the parsed parameters are
echoed in the JSON output so the user can confirm them.

```
FOLDBOUNDS-FRONTIER 2
sequence HPHHPHHPHHHH          # H/P string, length >= 2
q 3/2                          # natural or natural/positive-natural; checker rejects q < 1
threshold 4                    # event A = { k >= threshold }
<KIND> <WORD> [<WITNESS> ...]  # one line per frontier leaf
...
end
```

`WORD` is the prefix direction word (digits `0`–`3` = the four lattice
directions of `LatticeModel.advance`: `0` = +x, `1` = −x, `2` = +y,
`3` = −y; `.` = empty word / root).  The word lists the steps *after*
the fixed initial bond `(0,0) → (1,0)`.  Witnesses
are suffix words completing the prefix to a full conformation.

| kind | meaning | what the verifier recomputes |
|---|---|---|
| `E` | exact complete conformation | length = steps, self-avoidance, exact weight |
| `D` | dead prefix | no legal child exists |
| `X` | exhaustive local evaluation (v2) | DFS over the subtree; exact masses |
| `B` | version-1 basic bound | `3^{n-m} q^{K}` caps, witness checks |
| `T` | version-1 transfer bound | prefix-obstacle transfer DP |
| `M<k>` | memory-`k` transfer + geometric caps + families (v2) | memory-`k` DP, `geomCap`, degree cap, family count, witness checks |

Witnesses (only after `B`, `T`, `M<k>`) must be distinct, of the right
length, self-avoiding; duplicates or invalid ones cause rejection.  For
`M<k>` leaves, witnesses outside the certified constructive family are
added to the family mass (disjointness kernel-checked); witnesses inside
the family are ignored for the lower bound.

The flat list must form an exact prefix-free cover: duplicated entries,
ancestor/descendant overlaps, omitted legal branches and entries on
colliding prefixes are rejected (`buildTrie2` + `checkNode`).

A `FOLDBOUNDS-FRONTIER 1` header restricts the kinds to `E D B T`;
for such certificates `foldcheck2` returns exactly the version-1 result
(kernel test `legacy_agrees_with_v1`; 3000 differential certificates,
`experiments/lean_v2_adversarial.py`).

Output (one JSON object per file): echoed `version, sequence, q, threshold,
nodes`, `accepted`, and if accepted the eight rational endpoints
`zLo zHi eLo eHi nLo nHi pLo pHi` and `verify_ms`.

Guarantee (from `checkFlat2_sound`): if accepted, the exact anchored
(no reflection quotient) HP-model values for the parsed `(sequence, q,
threshold)` satisfy `zLo ≤ Z ≤ zHi`, `eLo ≤ Z_A ≤ eHi`,
`nLo ≤ Z_¬A ≤ nHi`, `pLo ≤ Z_A/Z ≤ pHi`.
