# Non-enumerative bounds: research report

The question was: *can a machine-verified algorithm rigorously bound HP
lattice-protein partition functions and probabilities without enumerating
every conformation?*

This report separates five kinds of statement: (A) theorems proved in Lean,
(B) measured computational behaviour, (C) negative results, (D) conjectures
and open directions, and (E) novelty, which has **not** been established.
Numbers in §B come from `experiments/results/lean_frontier_benchmark.log` and
`.json`, produced by `experiments/lean_frontier_benchmark.py`.

## A. What is proved (kernel-checked, standard axioms only)

1. **Compact certificates are sound.** `FrontierCheck.check_sound` and
   `checkFlat_sound`: any certificate the checker accepts yields rational
   intervals that contain the exact `Z`, `Z_A`, `Z_¬A`, and `P(A)` of the
   Lean HP model. A bound leaf stands for its entire completion set, so the
   checker never lists those completions.
2. **Compactness.** `check_root_bound_accepted`: a one-node certificate is
   accepted for every model. `checkNode` recurses structurally over the
   submitted trie only.
3. **Exact partition structure.** `ext_eq_biUnion`, `ext_pairwiseDisjoint`,
   `sum_ext_eq_sum_leaves`, and `leaf_cover_unique`: an accepted frontier
   partitions the conformation space.
4. **Three upper-bound levels with a proved ordering.**
   * basic: `3^r q^Kmax` (`card_ext_le_pow_three`, `contacts_le_kmax`);
   * degree-capped: `3^r q^(known + futureCap)` (`contactsAt_le_slots`,
     `contacts_le_degree_caps`), proved never worse than basic
     (`basicZ_degree_le_basic`);
   * transfer: `wsum` over the merged DP of relaxed walks that avoid the
     prefix and never reverse, with weights `q^min(known + cross + ffCap,
     known + futureCap)` (`wsum_dp`, `relaxed_sound`, `sum_ext_le_dp`),
     proved never worse than degree-capped (`transfer_total_le_basic`).
   The transfer step is a weighted generating-function bound in the sense of
   §6.1 of the assignment. The table is the exact cross-contact histogram of
   the relaxed walk set, and the envelope is applied coefficient by
   coefficient.
5. **Constructive lower bounds without enumeration.** `family_lower` and
   `famRaise_sound`: when the prefix end is extremal for a diagonal potential,
   the `2^r` directed continuations are certified completions, so
   `Z(p) ≥ 2^r q^known`. Explicit witnesses are also supported. They are
   distinct and re-scored by the checker, and each one tightens the upper
   bound by `envelope − exact` (`boundEnc_sound`). Witness and family lower
   bounds are combined with `max`, which is safe under overlap.
6. **Model foundations.** The representation bijection, sum preservation,
   nonemptiness, checkerboard parity, and Lean ≡ Python contact semantics;
   see `docs/ARISTOTLE_LATTICE_VERIFICATION_REPORT.md`.

## B. Measured behaviour (empirical, not proved)

Setup: all-H chains with `q = 2` and threshold `⌊n/4⌋`; the mixed sequences use threshold 2. Depth is the number of
residues fixed at each leaf. The Python producer builds the frontier, TW adds
DFS witnesses with budget 60 per leaf, and the Lean checker is the compiled
`foldcheck` executable. Verification time is measured inside the Lean process
around `checkFlat`. "exact" comes from an independent exhaustive Python
enumeration and is shown only for `n ≤ 16`. Every accepted interval contained
the exact value (`encl=True`) in every case where the exact value was
computed.

| n | depth | leaves | mode | verify | Z interval | exact Z |
|---|---|---|---|---|---|---|
| 14 | 8 | 543 | B | 0.1 s | [0, 1.2e14] | 1.53e6 |
| 14 | 8 | 543 | T | 15 s | [6.6e4, 1.9e7] | 1.53e6 |
| 14 | 8 | 543 | TW | 21 s | [1.8e5, 1.7e7] | 1.53e6 |
| 16 | 8 | 543 | B | 0.1 s | [0, 9.0e18] | 1.65e7 |
| 16 | 8 | 543 | T | 63 s | [2.6e5, 5.9e9] | 1.65e7 |
| 16 | 8 | 543 | TW | 68 s | [4.7e5, 5.8e9] | 1.65e7 |
| 18 | 2 | 1 | T | 1.2 s | [6.6e4, 2.2e15] | — |
| 18 | 8 | 543 | B | 0.2 s | [0, 2.7e24] | — |
| 18 | 8 | 543 | T | 166 s | [1.1e6, 9.0e11] | — |
| 18 | 8 | 543 | TW | 176 s | [1.35e6, 9.0e11] | — |
| HPHPPHHPHPPHPH (14) | 8 | 543 | TW | 17 s | [4.3e4, 1.35e6] | 3.80e5 |
| HPHPPHHPHPPHPHHPPH (18) | 8 | 543 | T | 91 s | [5.9e5, 8.3e8] | — |

The full table, including the mixed sequences, peak memory, certificate
bytes, and producer time, is in the log and JSON files. It covers 78 runs,
and none of them showed a problem. Peak RSS of the checker process stayed
below about 120 MB. Observations:

* The transfer relaxation lowers the upper bound by 2 to 13 orders of
  magnitude compared with the basic `3^r q^Kmax` envelope, and the gain grows
  with chain length. At `n = 14`, depth 8, it is
  within a factor of about 12 of the exact `Z`.
* The certified family gives nonzero `Z` lower bounds from a single node.
  For H18 the root node alone gives `Z ≥ 65536`, about 190 times the old
  Python lower bound of 343 reported in the assignment. At depth 8 the bound
  is `1.1e6`.
* Verification cost: basic mode is essentially free. Transfer mode costs
  about 0.03–0.3 s per leaf in the compiled checker, and grows with the
  remaining length. It is always independent of the number of completions.
* Comparison with the existing Python stack (H18, `q = 2`, a shallow
  frontier, `343 ≤ Z ≤ 2.98e10`): our verified transfer mode at depth 8 gives
  `1.1e6 ≤ Z ≤ 9.0e11`. The lower bound is much better and the upper bound is
  about 30 times *worse*, because the Python matching, isoperimetric, and
  bounding-box caps are not formalised. The two configurations differ, so
  this is **not** a controlled comparison.

## C. Negative results and limitations (honest assessment)

* **The probability intervals are still almost uninformative** for
  moderate thresholds at these depths (`P ∈ [≈0, ≈1]`). A leaf contributes to
  the event or complement *lower* bound only when every completion is
  certainly on one side. Shallow leaves rarely satisfy this, so `eLo` and
  `nLo` stay tiny. Informative probabilities would need far deeper
  frontiers, exact leaves, or new lower-bound families whose contact scores
  are controlled.
* The `Z` intervals still have ratios of `10^2`–`10^6`. The relaxed walks
  still allow self-intersection among the new sites, and the family lower
  bound grows like `2^r` against roughly `2.64^r` for all completions, so the
  gap grows exponentially with the number of unfolded steps.
* The transfer DP is polynomial per leaf but is the dominant cost
  (`O(r³·n·C)` states and work, see the certificate specification §7).
  There is no formal complexity theorem.
* An earlier version of the checker used Lean's core `List.mergeSort`, which
  the kernel cannot evaluate in `decide +kernel`. It was replaced by a
  fuel-based structural sort (`msortF`) with a proved permutation lemma.
  This changes speed only, not soundness.

## D. Conjectures and open directions (unproved)

1. **Monotone memory hierarchy.** Adding finite memory to the relaxed walk,
   for example forbidding returns to the last `k` new sites, should give
   relaxations `R_k ⊇ R_{k+1} ⊇ … ⊇ C(p)` whose weighted sums are monotone in
   `k`. The generic `wsum_dp` theorem applies to any such state machine.
   What remains to prove is the inclusion of each genuine completion into
   `R_k`, together with a simulation between consecutive levels. Not done.
2. **Partially directed families.** Families of walks using `{+x, +y, −y}`
   without vertical reversal grow like `(1+√2)^r` and are self-avoiding. A
   contact-weighted DP over them would give much stronger lower bounds,
   including event-classified ones. Not done.
3. **Witness/frontier optimisation.** No theorem relates certificate size to
   the achievable probability resolution. We have no proven trade-off.
4. Formalising the Python matching and isoperimetric caps would combine them
   with the transfer bound through `min`. Since every cap is sound, the `min`
   is too. Not done.

## E. Novelty assessment

No literature search was possible in this environment, which had no network
access. Therefore:

* Representation, parity, counting, and contact results: these are
  established, elementary facts, **newly formalised** here in Lean.
* Transfer-matrix DP over walks that avoid a fixed obstacle set, and
  branch-and-bound over SAW prefix tries: these are established techniques
  in lattice-polymer enumeration. Our contribution is a **formally verified**
  combination: a proof-carrying checker with a kernel-checked end-to-end
  soundness theorem, plus a proved hierarchy (transfer ≤ degree ≤ basic).
  Whether such a verified checker exists in the literature is
  **unverified**.
* The diagonal-family lower bound is elementary (directed walks are
  self-avoiding). Its use as a certified, non-enumerative lower bound inside
  a verified checker may be new, but **novelty is not established**.
* We do not claim a performance improvement over established exact
  enumeration or transfer-matrix methods. Exact enumeration of `n ≤ 18` HP
  chains is routine, and our bounds are much looser than exact values. The
  value of this work is assurance: every number the checker outputs is
  backed by a theorem. Speed or precision is not where it competes.
