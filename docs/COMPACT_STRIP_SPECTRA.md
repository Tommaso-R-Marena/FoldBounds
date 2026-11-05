# Compact Hamiltonian-strip certificates for the anchored HP model

**Status:** Python exact-integer/rational research prototype, model-relative. New work is intentionally additive to FoldBounds PR #63; no existing lattice solver, protein-motion code, Lean proof, CI file, or other agent branch is changed. This is **not** a discovery in experimental protein folding or a claim of priority for strip/transfer-matrix polymer combinatorics.

## Construction

For any tile height `h` between 3 and 8, a tile is the entire `4 × h` rectangular set of square-lattice sites. Enumerate its Hamiltonian paths starting at `(0,0)` and ending at `(3,0)`. The first tile is restricted to paths whose second site is `(1,0)`, matching FoldBounds' canonical first edge. Translate each subsequent tile four columns East, and join its entrance to the previous tile's exit by one East backbone edge. Tile site sets are disjoint and fully filled; a final straight East ray accommodates lengths not divisible by `4h` without creating new nonconsecutive contacts. Different selections of the finite tile paths produce different valid anchored SAWs.

The tile library is exhaustively constructed for a **fixed, bounded** tile height. A separate path validator rejects duplicate sites, missing cells, invalid endpoints, or nonadjacent steps. The Python library has no *a priori* claim that larger rectangles can be solved efficiently.

## Exact heterogeneous HP contact-score spectrum without enumerating combinations

For a given HP sequence, each tile variant determines: (a) its internal HH contact score, (b) a bitmask of H occupation along its left boundary, and (c) a bitmask along its right boundary. A dynamic program carries `(last_right_H_mask, score, direction_flags)` and exact integer multiplicities. Joining tiles creates precisely one backbone edge at row `y=0`. For every other row `y=1,...,h-1`, the seam contributes one new nonconsecutive HH contact iff both adjacent boundary sites are H. Thus the extra contact score is `bit_count((previous_right & current_left) >> 1)`.

Every actual tiled conformation has a unique sequence of tile variants, so the recurrence computes the exact contact-count histogram of the *restricted tile family*. Runtime scales with the DP states rather than the exponential number of assembled walks, although fixed-tile enumeration and worst-case score-state growth remain limitations.

For a complete all-H rectangle made of `B` tiles, all variants have the same number of HH contacts:

    n = 4 h B
    k = (h-1)(4B-1)

This follows exactly from counting occupied grid edges and subtracting the `n−1` backbone edges. At height `h=8`, the **finite verified tile library** has 360 Hamiltonian variants, 147 of which start with the canonical East edge. Therefore it constructs

    147 · 360^(B-1) distinct full-rooted folds

of contact score `28 B − 7`, for every integer `B >= 1`. Reflection about the line `y=0` preserves the canonical East first edge, self-avoidance, and all HH contacts. Every tiled conformation occupies positive y sites and its reflection negative y sites, so the two families are disjoint. Enabling `include_reflection=True` therefore doubles the witness multiplicity to **294·360^(B−1)** for each `B>=1`. These conclusions require only the explicit 32-site tile witnesses and their independently validated Cartesian concatenation; they do not assume an unproved general formula for tile counts. For example, when `B=2` there are **52,920 unreflected and 105,840 including reflected distinct conformations** of length 64 with **49 HH contacts each**. The elementary box/perimeter bound `k ≤ 64+1−ceil(2√64) = 49` shows that each is maximally contacting in this 2D all-H model.

For any rational `q >= 1`, this certifies the exact *restricted-family* mass and a rigorous lower bound on the total HP partition function:

    Z_(32B)(q) >= 294 · 360^(B−1) · q^(28B−7)   (including reflections)

Consequently, along `n=32B`, `liminf_B log(Z_(32B)(q))/(32B) >= [log(360) + 28 log(q)]/32`. This is **only** a model-relative analytic lower bound, not a new physical thermodynamic limit theorem.

## Certified disjoint addition to the existing partially directed witnesses

FoldBounds' PR #56 partially directed subset is the disjoint union of:

- East-progressive walks (no West steps);
- North-progressive walks with at least one West step (no South steps);
- South-progressive walks with at least one West step (no North steps).

The new `disjoint_from_partially_directed=True` certificate counts **only** compact paths containing at least one West, one North and one South step. Such a path satisfies none of the three partially directed languages. Therefore these compact paths and the PR #56 paths are provably disjoint, and their integer contact histograms and weighted masses can be **added** without overlapping witnesses. The code checks the direction mask for each finite tile combination through the transfer recurrence, rather than assuming every tile variant is disjoint. Horizontal reflection exchanges North/South but preserves the joint W/N/S predicate, so the mirrored family is also disjoint from PDSAW and from the original compact family. It deliberately does not add the entire unfiltered compact family or multiple competing tile heights.

The `compact_plus_directed_mass` consumer recomputes the existing partially directed spectrum and combines it with the separately checked disjoint compact certificate. It does not modify the preexisting thermodynamic-response calculation; reviewers can add it as an independent trusted input after auditing the exact bridge.

## Controlled model-only benchmark

All-H `n=64`, `q=10`, event threshold `k>=1`:

| Witness family | Exact lower mass |
|---|---:|
| Existing partially directed family | 13,394,219,636,175,161,757,789,276,985,333,749,427,840,163,576,552,435 |
| New disjoint compact family (including reflections) | `105,840 × 10^49` |
| Union (verified disjoint) | prior mass + new mass = **80.01916×** the previous partially directed mass |

At `q=100`, the new compact subset alone is about **6,822×** the partially directed lower mass in the same model. These factors are improvements in *proven lower bounds*, not tighter upper bounds, accuracy, a practical folding-prediction benchmark, or an empirical discovery. The true partition sum may exceed both lower bounds by a large factor. At `q=2` the new subset can have relatively little weight and the previous method may dominate.

## Independent tests and assurance boundary

Tests use an independently coded bitmask Hamiltonian-path enumerator to reconstruct each tile library (`h=3..8`); direct Cartesian enumeration with an independent all-pairs HP contact scorer for small `B`; full rooted SAW enumeration at length 12; error injection into supplied spectra; brute-force direction-language disjointness checks; exact `q` tests; and a nonenumerative long-tile test. The PR #56 partially directed module is a Python reference and is itself **not formally refined to Lean**.

    python -m pytest tests/test_compact_strip_spectra.py -q
    python experiments/benchmark_compact_strip_spectra.py

**Remaining work:** prove the finite tile-library/checker interface in Lean, integrate disjoint subset masses into the richer moment certificate API without double-counting exact-frontier witnesses, compare against published strip-polymer baselines, and most importantly evaluate physically calibrated protein structure/thermodynamics outside the toy HP model. No claim is made about experimental proteins.

## Continuum (all-q-tail) improvement certificate

If the minimum contact exponent in the filtered compact spectrum is **at least** the maximum exponent of the PDSAW spectrum, then the *ratio of the disjoint compact mass to the PDSAW mass* is nondecreasing in `q>=1`. Proof: for each compact exponent `a` and PDSAW exponent `b`, the exponent difference `a-b` is nonnegative, so multiplying the nonnegative finite sums shows `C(q2)P(q1) >= C(q1)P(q2)` whenever `q2>=q1>=1`. Therefore `1+C(q)/P(q)` is nondecreasing. This is a simple finite polynomial order argument, not an original general theorem.

`build_tail_dominance_certificate('H'*64, height=8, q_min=Fraction(10))` recomputes the spectra and confirms this exponent-order condition. Its separately replayed rational certificate then establishes that **for every real `q>=10`, not just sampled integer q**, the PDSAW-plus-compact lower mass is at least `80.01916...` times the original PDSAW lower mass. The factor approaches `1+105840/2=52921` as `q` increases without bound, because the existing PDSAW coefficient at contact score 49 is two. This continuum certificate uses the default mirrored disjoint subset and independent tests check the reflection geometry explicitly. The claim concerns a ratio of lower bounds, not a bound on the unknown true partition-function error.
