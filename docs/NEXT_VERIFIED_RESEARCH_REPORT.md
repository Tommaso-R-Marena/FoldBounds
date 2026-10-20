# FoldBounds — version-2 verified bounds: research report

This report covers the version-2 extension of the FoldBounds Lean checker
built on top of the integration snapshot available in this environment
(the `FoldBounds-main/` tree, which contains the previous Aristotle
frontier work: `Frontier*.lean`, `LatticeEquivalence.lean`, `Transfer*.lean`,
117-declaration audit).  The newer GitHub commit
`7ab6349542ec4b2dd0cda47daa029228f304aa16` could **not** be inspected
(no GitHub access was attempted, as instructed).  Python files that may
have changed upstream were therefore *not* modified; all Python additions
are new files.

Everything labelled **Theorem** below is a Lean declaration that compiles
under Lean 4.28.0 / Mathlib v4.28.0, contains no `sorry`/`admit`, and is in
the axiom audit with only `propext`, `Classical.choice`, `Quot.sound`
(see `VERIFICATION_REPORT.md`).  Everything labelled **Measured** is an
experiment, not a proof.

---

## 1. The mathematical problem

Model (unchanged from version 1; see `docs/LEAN_LATTICE_MODEL.md`): an
H/P sequence `S` of length `n = steps + 2`; conformations are anchored
square-lattice self-avoiding walks starting `(0,0) → (1,0)`, encoded by
direction words of length `steps` (`LatticeModel.validWords`); **no**
rotation/reflection quotient.  `k(c)` is the number of H–H pairs `i < j`,
`j ≥ i + 3`, `j − i` odd, that are lattice-adjacent
(`parityContactCount`; proved equal to the direct definition in the
previous run).  For rational `q ≥ 1`, weight `q^{k(c)}`.  For a threshold
`t`, event `A_t = {k ≥ t}`.  We want exact-rational enclosures of

  `Z = Σ_c q^{k(c)}`,  `Z_A`, `Z_{¬A}`,  `P(A) = Z_A / Z`

from an *untrusted* certificate, without enumerating every conformation.

## 2. Precise theorem statements and hypotheses

### 2.1 Flagship (version 2)

`FoldBounds.FrontierNext.check_sound`
```
theorem check_sound (steps : ℕ) (seq : Fin (steps + 2) → Bool) (q : ℚ) (thr : ℕ)
    (t : NTree) (r : FrontierResult) (h : check (List.ofFn seq) q thr t = some r) :
    massQ (weightQ steps q seq) univ ∈ Set.Icc r.zLo r.zHi ∧
    massQ (weightQ steps q seq) (contactEvent steps seq thr) ∈ Set.Icc r.eLo r.eHi ∧
    massQ (weightQ steps q seq) (contactEvent steps seq thr)ᶜ ∈ Set.Icc r.nLo r.nHi ∧
    massQ (weightQ steps q seq) (contactEvent steps seq thr) /
      massQ (weightQ steps q seq) univ ∈ Set.Icc r.pLo r.pHi
```
Hypotheses: none beyond acceptance.  The checker itself rejects `q < 1`
and sequences of length `< 2`.  `t : NTree` is arbitrary untrusted data.
`check_coupled_probability_sound` restates the probability claim over ℝ in
the form `Z_A/(Z_A+Z_¬A)`, matching `LatticeModel.checkHP_coupled_sound`.
`FrontierNextFlat.checkFlat2_sound` gives the same conclusion for the flat
(list of `(prefix, kind)` lines) certificate used by the CLI.

Probability endpoints use the previously proved coupled bound
`L_A/(L_A+U_¬A) ≤ P ≤ U_A/(U_A+L_¬A)` with the version-1 zero-denominator
policy (endpoint `0` resp. `1` when the corresponding denominator is `0`).

### 2.2 Geometric contact caps (`FoldBounds/VerifiedGeometry.lean`)

* **Theorem `Geometry.contacts_le_rows_cols`** (perimeter inequality).
  For any list `W` of distinct lattice sites forming a nearest-neighbour
  chain (`adjacentChain W = true`) with `W.length ≥ 1`:
  `#contactPairs W + #rows W + #cols W ≤ W.length + 1`,
  where `contactPairs` are index pairs `i < j`, `j ≠ i+1`, at lattice
  distance 1, and `rows`/`cols` are the sets of occupied y/x coordinates.
  Proof: count lattice edges inside the occupied set two ways — at most
  `N − #rows` horizontal and `N − #cols` vertical edges (each row/column
  segment loses one), and at least `(N−1) + #contacts` (backbone + contacts).
* **Theorem `Geometry.rows_ge_span` / `cols_ge_span`** (connectivity):
  occupied rows form an interval, so `#rows ≥ y_b − y_a + 1` for any two
  occupied sites; proof by a discrete intermediate-value theorem
  (`discrete_ivt`).  Hence `#rows`, `#cols` *are* the bounding-box height
  and width, and the user's candidate `K ≤ max(0, n + 1 − w − h)` (w,h in
  lattice sites) is proved exactly, without extra hypotheses.
* **Theorem `Geometry.length_le_rows_mul_cols`**: `N ≤ #rows · #cols`.
* **`geomCap N K := N + 1 − minPerim N (xSpan K) (ySpan K)`**, where
  `minPerim N w h = min { c + max(h, ⌈N/c⌉) : w ≤ c ≤ N }`.
  **Theorem `Geometry.contactPairs_le_geomCap`**: if every site of `K`
  lies on the chain, `#contactPairs W ≤ geomCap W.length K`.
  **Theorem `Geometry.parityContacts_le_geomCap`**: every valid complete
  word whose walk contains `K` has `k ≤ geomCap (steps+2) K`.
  Example: `geomCap 18 [] = 10` (all-H 18-mer has at most 10 contacts
  by geometry alone; the true maximum is smaller, not claimed).

### 2.3 Finite-memory transfer relaxation (`FoldBounds/VerifiedMemoryTransfer.lean`)

State `MState = (cur, dir, hist, cross, ffLo, ffHi, idx)`: current site,
last direction, the last `k` placed suffix sites with their residue
indices (`hist`), the exact number of suffix-to-prefix H–H contacts
observed (`cross`), and a certified lower/upper count of suffix-to-suffix
contacts (`ffLo`, `ffHi`).  Transition `memStep`: forbids prefix sites and
the remembered `k` sites (`memAllowed`); for the new residue `j`,
`loAt` counts contacts with *remembered* earlier H residues that are
actually adjacent (observed ⇒ valid lower bound), `upAt = loAt + outAt`
where `outAt` counts parity-compatible earlier H residues that have fallen
out of memory (possible but unobserved contacts ⇒ upper bound).  The
two directions are never interchanged: `ffLo` is only ever increased by
observed contacts; `ffHi` additionally pays for every unobservable one.

* **Theorem `MemTransfer.memState_inv`** (simulation invariant): along any
  genuine self-avoiding completion, the relaxed state `gstate` reached is
  well-formed, `hist` is a genuine suffix of the walk, `cross` is the exact
  cross-contact count, and `ffLo ≤ ff ≤ ffHi` for the genuine
  suffix–suffix count `ff`.
* **Theorem `MemTransfer.memCover`** (coverage + weighted domination):
  for every nonnegative state function `G`,
  `Σ_{l ∈ ext steps p} G(gstate S k p l) ≤ wsum G (memTable steps S k p)`.
  Every genuine completion is a relaxed path; the relaxed table may also
  contain invalid paths, which only increase the right side.
* **Theorem `MemTransfer.wsum_dpM`** (DP correctness): the executable,
  stack-safe dynamic program (merge-sort based state merging, `combineTR`)
  computes exactly the weighted sum over all relaxed paths
  (iterated `memNext`), for any weight function.
* **Theorem `MemTransfer.memContacts_bounds`**: for every genuine
  completion `l` with final relaxed state `st`:
  `known + cross + ffLo ≤ k(l) ≤ known + cross + ffHi` and
  `k(l) ≤ geomCap (steps+2) (knownSites p st)`.
* The checker uses `mExp = min(known+cross+ffHi, kcN, geomCap(knownSites))`
  (upper exponent) and `mLow = known+cross+ffLo` (lower exponent):
  event mass bound counts a state only if `thr ≤ mExp`; complement mass
  bound counts a state only if `mLow < thr`.

### 2.4 Monotone refinement (`FoldBounds/CertifiedRefinement.lean`)

* **Theorem `Refinement.memUppers_mono`**: for `1 ≤ q` and `k ≤ k'`,
  all three memory-`k'` leaf upper bounds (Z, Z_A, Z_¬A) are `≤` the
  memory-`k` ones.  Proof via a simulation relation `Rel` between
  memory-`k'` and memory-`k` states (`step_rel`, `spec_mono`): more memory
  forbids more moves, observes more contacts (`loAt_mono`), and pays for
  fewer unobserved ones (`lo_add_out_anti`, `upAt_anti`); geometric caps
  are antitone in known sites (`geomCap_anti`, `minPerim_mono`).
* **Theorem `Refinement.memUppers_le_basic`**: for every memory `k`, the
  memory leaf uppers are `≤` the version-1 basic and degree-capped
  counting bounds `3^{n−m} q^{K}`.
  *Not proved:* comparison with the version-1 prefix-obstacle *transfer*
  bound (`T` leaves).  Experimentally memory leaves were never worse
  (all rows in §7), but this is measured only.

### 2.5 Constructive families (`FoldBounds/VerifiedFamilies.lean`, `FrontierNext.lean`)

* Partially directed family: if the prefix endpoint is the unique strict
  maximum of the prefix along axis direction `m` (`pdOK m p`), all words
  generated by `pdWords m r 0` (steps in `m`, `up`, `down`, never `up`
  immediately followed by `down` or vice versa) extend `p` to genuine
  self-avoiding conformations.
  **Theorem `Families.pdFamily_subset_ext`**, **`card_pdFamily`**
  (`#family = pdCount r`, computed by a 3-term linear recurrence
  `pdTriple`), **`two_pow_le_pdCount`** (`2^r ≤ pdCount r`; growth
  `≈ (1+√2)^r` versus `2^r` for the previous diagonal family).
* Every family member has `k ≥ known prefix contacts`; so its mass
  `pdCount(r) · q^{known}` is a certified lower bound for `Z`, for `Z_A`
  if `thr ≤ known`, and for `Z_¬A` if `kcN < thr`.
* **Disjoint family + witness combination** (`FrontierNext.fam_wit_lower`,
  `memEnc_lower_sound`): explicit witness completions that fail the
  executable membership test `famTest` (proved complete for the family:
  `mem_pdFamily_test`, `inFam_of_mem_family`) are disjoint from the family
  and their exact weights are *added*.  Witnesses inside the family are
  dropped, never double counted.  Kernel example
  `FrontierNextTests.accept_root_mem3_witness`: HHHHHH, q=2: family `41`
  plus disjoint witness `4` gives `Z ≥ 45 > max(41, 4)`.
* Witness upper-bound correction: each witness's envelope contribution is
  replaced by its exact weight (`memEnc_upper_sound`).

### 2.6 Exhaustive local leaves

`NTree.enum` (`X` in text) asks the verifier to evaluate a subtree by
depth-first search (`enumEnc_sound`, exact).  This makes the cost/precision
trade-off explicit inside one certificate: the proposer chooses which
subtrees are bounded (`M<k>`) and which are evaluated (`X`).

## 3. The algorithm and why it is sound

The checker walks the untrusted trie; at each node it re-validates
geometry (`(walkFrom p).Nodup`, lengths), and at each leaf computes an
enclosure using only verified functions:

| leaf | lower endpoints | upper endpoints | lemma |
|---|---|---|---|
| `E` | exact | exact | `exactEnc` |
| `D` | 0 | 0 (checked no legal child) | `deadOK` |
| `X` | exact (DFS) | exact | `enumEnc_sound` |
| `B`,`T` | version-1 | version-1 | `boundEnc`, `famRaise` |
| `M<k>` | family ⊎ outside witnesses | min(count×cap, memory-k DP envelope) − witness corrections | `memEnc_sound` |

Internal nodes add children; the children of a nonterminal prefix
partition its completion set (version-1 theorem `ext_partition`-style
lemmas reused), so sums of child enclosures enclose the parent.  At the
root the completion set is the whole conformation space.

## 4. Independently verified refinement

`memUppers_mono` (§2.4) is the monotone hierarchy: memory 0 ⊒ memory 1 ⊒
… ; at memory `≥ steps` every relaxed path is genuinely self-avoiding and
the observed contacts are exact, which is why measured M4 leaves on 4–8
remaining steps hit the exact Z.  (That memory `≥ remaining steps` gives
exactness is observed, **not** proved as a theorem.)

## 5. What is new relative to the previous run

1. Geometric perimeter/bounding-box contact cap, fully proved (new module).
2. Finite-memory transfer relaxation with observed-vs-possible contact
   separation, coverage, DP correctness and monotone refinement (new).
3. Partially directed constructive family (`≈2.414^r` vs `2^r`), plus a
   kernel-checked disjointness rule allowing family + witness *addition*.
4. Version-2 checker integrating all of the above with the old leaf kinds,
   flagship soundness theorem, flat/CLI front end, compiled `foldcheck2`.
5. 46 new audited declarations (117 → 163).

## 6. What remains unverified / not done

* **Matching / max-flow / dual certificates (§3.3 of the brief): not
  implemented in Lean.**  Only the existing version-1 degree budget
  (`kcap`) is used.  This is the most natural next step.
* Endpoint-conditioned terminal-box bound (per terminal location): not
  formalized; geometric caps use spans of known sites only.
* No proof that memory leaves are ≤ version-1 *transfer* leaves.
* No kernel-checked complexity theorem; the cost analysis below is
  informal.  Memory-`k` table size is bounded by
  `#states ≤ (#reachable sites) · 4 · (#site-sequences of length k) · (contact counters)`;
  measured: practical only for ≲ 10 remaining steps (root M2 for n = 18
  took ≈50 s; root M4 for n = 18 failed, apparently out of memory).
* The text parser (`FrontierNextCLI.lean`) and the JSON output of
  `foldcheck2` are an **unverified boundary**; the parsed parameters are
  echoed so the user can confirm they match the intended model.
* Python exporter `lean_export_v2.py` is untrusted; correctness of what it
  proposes is irrelevant to soundness (wrong certificates are rejected or
  give wider-but-sound intervals).
* The reals-valued physical Boltzmann model is not addressed (exact
  rational `q` only).

## 7. Benchmark improvements and failures (Measured)

Settings identical across modes: same frontier (Python `build_frontier`
at the stated depth), same sequence, `q`, threshold, same witness budget
(`propose_dfs_witnesses`, 60 nodes per leaf) for `*W` modes; Lean-compiled
verifier `foldcheck`(v1)/`foldcheck2`(v2) on the same machine.  `v1-B`,
`v1-T`, `v1-TW` = previous basic / transfer / transfer+witness leaves;
`M<k>` = version-2 memory leaves, `M<k>W` with witnesses.  Every accepted
row enclosed the exact value (oracle: independent Python enumeration,
n ≤ 16).  Full logs: `experiments/results/v2/lean_v2_benchmark_allH_partial.log`,
`experiments/results/v2/lean_v2_benchmark_mixed.log` (JSON alongside).

**All-H, q = 2 (Z upper bound / exact; verify time).** Bold = equal to the exact Z to the 4 significant digits printed in the log (this all-H run was stopped during n = 16 and only its log is kept; no JSON):

| n | thr | depth | nodes | v1-T Zhi | M0 Zhi | M2 Zhi | M4 Zhi | exact Z | M4 ms | v1-T ms |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 3 | 4 | 9 | 1.656e5 | 3.672e4 | 2.321e4 | **1.332e4** | 13323 | 385 | 164 |
| 10 | 3 | 6 | 71 | 2.309e4 | 1.768e4 | 1.332e4 | 1.332e4 | 13323 | 446 | 248 |
| 12 | 4 | 4 | 9 | 5.039e7 | 9.433e5 | 6.73e5 | 2.412e5 | 1.428e5 | 3278 | 676 |
| 12 | 4 | 6 | 71 | 1.787e6 | 5.335e5 | 2.784e5 | **1.428e5** | 1.428e5 | 4047 | 1701 |
| 14 | 4 | 6 | 71 | 5.474e8 | 1.269e7 | 8.618e6 | 2.713e6 | 1.533e6 | 31184 | 6695 |
| 14 | 4 | 8 | 543 | 1.927e7 | 6.548e6 | 3.163e6 | **1.533e6** | 1.533e6 | 38895 | 16100 |
| 16 | 5 | 6 | 71 | 8.34e10 | 3.23e8 | 2.296e8 | 8.423e7 | 1.646e7 | 207274 | 19934 |

Z lower bounds (partially directed family vs v1): e.g. n=14 d=8:
v1-T 6.611e4, M0 1.688e5; with witnesses v1-TW 1.837e5, M2W/M4W 3.214e5.

**Probability intervals** (all-H, q = 2, P(k ≥ thr)):

| n | thr | d | v1-TW | M4W | exact |
|---|---|---|---|---|---|
| 10 | 3 | 6 | [0.291, 0.640] | [0.349, 0.549] | 0.4431 |
| 12 | 4 | 6 | [0.0058, 0.996] | [0.0756, 0.896] | 0.4068 |
| 14 | 4 | 8 | [0.0133, 0.997] | [0.168, 0.948] | 0.5806 |

Without witnesses both v1-T and M<k> give `P ∈ [0,1]` (event lower mass 0)
in the all-H rows above: the constructive families only certify `k ≥ known`,
which is below threshold for most leaves.  **This is the main negative
result: the total-mass upper bound becomes exact or near-exact, but the
probability interval is limited by event lower mass.**

**Hybrid certificates** (`experiments/lean_v2_hybrid.py`; M4W leaves, but
leaves whose prefix already has ≥ κ contacts are evaluated exactly by `X`):

| case | κ | X leaves | verify ms | P interval | exact |
|---|---|---|---|---|---|
| H12 q2 thr4 d6 | none | 0/71 | 4713 | [0.0756, 0.8962] | 0.4068 |
| | 1 | 30/71 | 3165 | [0.3357, 0.5712] | |
| | 0 | 71/71 | 535 | exact | |
| HPHHPHHPHHHH q2 thr4 d6 | none | 0/71 | 4783 | [0.0053, 0.6287] | 0.1196 |
| | 1 | 24/71 | 3540 | [0.0904, 0.2310] | |
| H14 q2 thr4 d8 | none | 0/543 | 46842 | [0.1679, 0.9475] | 0.5806 |
| | 2 | 124/543 | 38806 | [0.4713, 0.8261] | |
| | 1 | 308/543 | 25195 | [0.5514, 0.6820] | |
| | 0 | 543/543 | 5311 | exact | |
| HPHPPHHPHPPHPH q2 thr2 d8 | none | 0/543 | 44952 | [0.0413, 0.8274] | 0.2908 |
| | 1 | 86/543 | 39027 | [0.1915, 0.6308] | |

**Honest negative result:** at these sizes (n ≤ 16) full exhaustive
evaluation inside the verifier (`κ = 0`) is *cheaper* than the memory-4
relaxation (e.g. 5.3 s vs 46.8 s for H14 d8).  The memory DP is a
precision tool, not a speed tool, at n ≤ 16; its benefit is that its cost
grows with the memory and remaining length rather than with the number of
completions, but in our implementation the constant factors dominate in
the tested range.  The all-H run was stopped after the n = 16, d = 6, M4
row (M4 took 207 s there); all-H n = 16 at d = 8 and n = 18 were **not run**,
and all-H n = 8 was not in the benchmark list.  An n = 18 mixed sequence
is reported in §7b.

### 7b. Mixed sequences (complete run, 120 rows)

Full table: `experiments/results/v2/lean_v2_benchmark_mixed_table.md`
(generated by `experiments/lean_v2_tables.py` from
`lean_v2_benchmark_mixed.json`).  Selected rows (P = P(k ≥ thr)):

| sequence | q | thr | d | v1-TW Z | M4W Z | v1-TW P | M4W P | exact P | M4W ms | v1-TW ms |
|---|---|---|---|---|---|---|---|---|---|---|
| HPHHPHHPHHHH | 1 | 4 | 6 | [4432, 3.326e4] | [1.44e4, 3.007e4] | [0.0007, 0.864] | [0.0008, 0.172] | 0.0174 | 4738 | 1944 |
| HPHHPHHPHHHH | 3/2 | 4 | 6 | [5515, 1.542e5] | [1.711e4, 5.0e4] | [0.0008, 0.963] | [0.0026, 0.427] | 0.0575 | 4714 | 2011 |
| HPHHPHHPHHHH | 2 | 2 | 6 | [7500, 5.226e5] | [2.06e4, 8.189e4] | [0.0126, 0.995] | [0.120, 0.948] | 0.6070 | 4392 | 1857 |
| HPHHPHHPHHHH | 2 | 4 | 6 | [7500, 5.226e5] | [2.06e4, 8.189e4] | [0.0008, 0.986] | [0.0053, 0.629] | 0.1196 | 4755 | 1978 |
| HPHPPHHPHPPHPH | 2 | 2 | 8 | [4.252e4, 1.348e6] | [1.275e5, 3.801e5] | [0.0080, 0.983] | [0.0413, 0.827] | 0.2908 | 45893 | 16887 |

At d = 6 the M4 upper bound for `HPHHPHHPHHHH` equals the exact Z
*exactly* (rationals 30073, 50005, 81891 for q = 1, 3/2, 2); v1-T gave
3.326e4, 1.622e5, 5.56e5.

Across all 120 mixed rows and all comparable mode pairs (v1-T vs M0/M2/M4,
v1-TW vs M2W/M4W), **no version-2 endpoint was worse than the version-1
endpoint** (checked by script over all eight endpoints); this is measured,
not proved (only `memUppers_le_basic` is proved).

**Failure / scaling limit (n = 18, `HPHPPHHPHPPHPHHPPH`, d = 10, 4073
leaves, 400 s per-run timeout):** v1-T accepted in 324 s
(Z ∈ [1.161e6, 1.024e8]); M0 accepted in 371 s
(Z ∈ [4.654e6, 1.010e8], i.e. 4× better lower bound, similar upper);
v1-TW, M2, M4, M2W, M4W **timed out**.  The witness-carrying certificates are
≈1.3 MB.  One row is flagged `v1_v2_compatible = false`: there the v1
executable timed out while the v2 executable's result was not recorded;
this is a timeout artefact, not an observed value mismatch (all 3000
completed differential v1/v2 comparisons agreed exactly).

## 8. Scientific significance (assessment)

* The main contribution is engineering-plus-proof: a sound, versioned,
  proof-carrying checker whose total-mass upper bounds are dramatically
  tighter than the previous verified baseline (often by 1–3 orders of
  magnitude at equal frontier, sometimes exact) without visiting every
  conformation, with a kernel-checked monotone refinement theorem.
* It does **not** yet deliver narrow probability intervals for moderate
  n without substantial exhaustive evaluation; for n ≤ 16 exhaustive
  verification remains competitive or faster.
* No claim of polynomial-time protein folding is made or implied.

## 9. Literature and originality (uncertain)

* The HP lattice model and exact enumeration of short lattice proteins are
  classical (Lau & Dill; many exact-enumeration studies).  Exact
  enumeration of self-avoiding walks by transfer-matrix / finite-lattice
  methods is far more powerful than anything here.
* The perimeter/edge-counting contact bound is a standard argument for
  compact lattice polymers; we believe the *statement* is known and only
  its formalization is new to this project.
* Relaxing self-avoidance to a finite memory window (walks avoiding only
  short loops) is a classical technique for bounding connective constants;
  our weighted variant with separated observed/possible contact counters
  inside a proof-carrying checker is, to our knowledge, not in the
  formal-methods literature, but **originality has not been verified** by
  a literature search in this environment.
* Partially directed walks and their `(1+√2)^n` growth are classical.
* We do not claim novelty for any mathematical statement; the plausible
  contribution is the combination: a kernel-verified certificate checker
  for statistical-mechanical bounds with these relaxations.

## 10. Strongest next directions

1. **Proof-carrying matching/dual bounds:** have Python emit a vertex
   cover / fractional dual for the residual H–H contact graph; verify in
   Lean only that it is a valid cover (cheap) and prove
   `#future contacts ≤ cover weight`.  This avoids verifying max-flow.
2. **Event lower mass:** families with *guaranteed* contacts (e.g.
   hairpin/zig-zag motifs with a counted family generating function), or
   an exact-contact transfer DP restricted to a strip where
   self-avoidance is decidable by finite memory (strip width ≤ memory),
   giving lower bounds with exact contact spectra.
3. Prove `memory ≥ remaining steps ⇒ exact` and `M ≤ T`.
4. Coefficientwise spectra in the memory DP (track `(ffLo, ffHi)`
   histograms rather than maxima) for joint event/complement bounds.
5. Performance: hash-consing states, avoiding re-sorting, compiled
   `Std.HashMap` with a proved refinement to the list-based DP.
