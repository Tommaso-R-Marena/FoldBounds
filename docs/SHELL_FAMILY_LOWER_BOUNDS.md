# Fixed-core, disjoint extremal-quadrant lower-mass shells

**Status:** A mathematically motivated, exact-rational **Python research**
algorithm tested against independent complete-geometry enumeration. **Not yet
Lean verified**, and not an assertion about real proteins or scientific
novelty. This module is intentionally developed on its own branch, separate
from the Aristotle frontend-checker integration.

## A non-enumerative lower mass at a fixed core depth

Let a canonical self-avoiding path `p` have `m` vertices, with final vertex
`v`. For each sign pair `(sx,sy)` in `{+1,-1}²`, define
`phi(x,y)=sx*x+sy*y`. A sign pair is *eligible* when
`phi(v)=max_{u in p} phi(u)`. Its quadrant suffix comprises `r=n-m`
unit steps, each either `(sx,0)` or `(0,sy)`.

**Self-avoidance:** Every new step increases `phi` by exactly one. Since all
fixed-prefix sites have `phi<=phi(v)` and each suffix site after j steps
has `phi(v)+j`, no suffix site can equal a prefix site. Distinct suffix
sites have distinct potential. Thus every completed walk is self-avoiding.

**Contact locality:** Two suffix vertices differing by two or more suffix
steps have `phi` difference >=2, whereas neighboring square-grid sites
have absolute `phi` difference 1. So nonconsecutive suffix/suffix
contacts are impossible. Similarly, a suffix site after two or more steps
cannot contact the prefix: every prefix vertex has potential <=phi(v).
All additional H-H contacts are therefore *exactly* those formed by the
first suffix vertex and earlier fixed H residues, excluding the backbone.

For each eligible quadrant and first-step cardinal direction `d`, there
are `2^(r-1)` suffix direction words, and all have exactly the same
contact score `k_fixed+k_cross(d)`. Their exact rational weight mass is
`2^(r-1)*q^(k_fixed+k_cross(d))`, for `q>=1` (indeed for any positive q).

**Disjoint union:** Two distinct signed quadrants intersect only along
their common cardinal ray (when they share one direction). Three
quadrants cannot share a complete nonempty suffix. Thus for every `d`,
subtract exactly one copy of its all-`d` ray when both quadrants containing
`d` are eligible. The implementation uses
`multiplicity(d)*2^(r-1) - 1_{multiplicity(d)=2}`, assigning its entire
weight to the correct event/complement class. The formula also holds for
`r=1`. `r=0` is handled by exact scoring of the core.

Enumerate all canonical self-avoiding cores of the **same** length m
extending a supplied prefix. Their completion sets are disjoint because
every completed direction word has exactly one m-vertex initial segment.
Summing the quadrant-union masses over all cores therefore gives **exact
restricted-family event and complement masses** that can both be used
as lower bounds on their unrestricted counterparts. No double counting
occurs across cores or inside each eligible quadrant union.

These masses can be intersected with independent previously certified
lower bounds using the `max` operation. Event and complement are disjoint
sets, so `max(event lower)+max(complement lower)` is also a sound lower
bound on total partition mass, even if the two chosen families differ.
Do **not** add this shell's event mass to overlapping two-turn, ray,
monotone, or independently enumerated witness masses.

## Complexity and scope

The search enumerates all valid **m-vertex cores**, and only those cores,
not the full n-vertex conformation space unless `m=n` is explicitly chosen.
There are at most `3^(m-2)` distinct cores after the anchored first edge.
For each core, the formula counts exponentially many valid completions
using exact integer multiplicities, without enumerating those suffixes.
The algorithm is exponential in m and can be impractical for large m.
It is not a polynomial-time folding algorithm. Its effectiveness varies
by sequence, q, threshold, and core depth.

## Reproducible verification

    python -m pytest tests/test_shell_family.py
    python experiments/compare_shell_family.py
    python -m foldbounds HPHHPHHPHHHH --q 3/2 --event-at-least 2 \
      --enhanced-lookahead 0 --matching --transfer --transfer-memory 4 \
      --terminal-box --monotone --turn-family --shell-core 10

Independent oracle tests explicitly enumerate complete small self-avoiding
walks, separately compute the eligibility predicate and their contact
scores, and compare **exact Fraction event/complement masses** with the
formula. For `HPHHPHHPHHHH`, `q=3/2`, threshold 2 and m=10,
`6867` is the restricted **event** lower mass from 17,227 counted full
paths across 4,067 seed cores (complement paths are also included). The
previous two-turn-family event lower mass was only `459/8`. This is a
restricted ensemble claim, not an assertion that true event mass is 6867.

Any probability interval improvement depends on the separately checked
upper bounds, particularly the complement upper mass. Exact small-model
oracle verification is a software test, not a Lean proof of the Python
implementation. Literature novelty and competitive performance remain
unverified.

## Monotonicity under core-depth refinement

The selected *full-path* shell sets satisfy `F_m ⊆ F_{m+1}` for
`2 <= m < n`, with identical sequence and starting prefix. To see this,
consider a complete monotone suffix in `F_m`. Its first appended point
lies in the same signed outward quadrant, and has potential strictly
higher than all sites in the previous core. Therefore the newly extended
`(m+1)`-vertex core has an extreme terminal point for that same signed
potential, and its remaining suffix is still quadrant-monotone. This
places the full conformation in `F_{m+1}`. Consequently, **both the event
and complement exact restricted-family masses are nondecreasing in m**,
for any fixed positive q. We test this property independently. It does
not mean increasing m is computationally cheap: core enumeration still
grows exponentially, and choosing `m=n` performs full enumeration.