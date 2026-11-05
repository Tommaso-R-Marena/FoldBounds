# FoldBounds — adaptive proof-carrying contact topology, November 5, 2026

## Registered biological outcome (unchanged)

The original RNase A 7RSA -> 2AAS challenge remains a **FAIL**. Weighted ANM Spearman was 0.421357, below the frozen 0.50 threshold, and trailed the best original GNM-resistance baseline (Spearman 0.511009) by 0.089652. None of the following results changes the registration or indicates a novel law of protein dynamics. The 32 deposited NMR structures are **not** independent samples of the protein's physical equilibrium or kinetic trajectory.

## Exact certification problem

The 7RSA structure has 124 Cα atoms and 2,614 eligible pair responses. Each Cα Cartesian coordinate is allowed to change independently by at most **0.05 Å**. Exact rational interval geometry of the *archived IEEE-754 coordinates* classifies the 7,626 unordered residue pairs into **972 certain contacts**, **108 optional contacts**, and **6,546 always-absent contacts**, using the 10 Å cutoff inclusive. The certain-contact network is connected.

A naive uncertainty certificate evaluates (i) the graph with all 108 optional edges and (ii) the graph with none. Rayleigh monotonicity encloses every contact realization's pair response between the two extremes. But applying these separate score intervals to the previously revealed 2AAS target ranks gave a signed directional concordance lower bound of **−15117/1138397 ≈ −0.0132792**: unable to certify positive agreement.

## Adaptive case decomposition

We select six of the optional contacts using **only the 7RSA reference structure and its previously frozen eligible pair list**, ranking each optional edge by summed squared effective-resistance sensitivity at the certain-contact graph. This sensitivity uses the established Sherman–Morrison rank-one identity; its floating inverse is used **only to select cases**, never accepted as a mathematical witness. The selected edges, zero-indexed Cα indices, are:

1. (51,115) — PDB residues (52,116)
2. (1,9) — residues (2,10)
3. (6,32) — residues (7,33)
4. (43,121) — residues (44,122)
5. (1,7) — residues (2,8)
6. (52,76) — residues (53,77)

Any physical contact-network assignment has **exactly one** 6-bit selected-edge mask. For each of the 64 cases, the minimum response is bounded by the network with **all other 102 optional edges ON** and the maximum response by the network with **all other 102 OFF**. This is a *superset* of the networks jointly realizable by coordinate perturbations, so the certificate remains valid even if optional edge switches are correlated through shared coordinates. No enumeration of the remaining 2^102 assignments occurs.

The independent Python verifier reconstructs every graph from an **integer-only dyadic coordinate decoder**, validates the fixed spanning tree (anchor 45, path constant B = 262), and checks all 128 untrusted integer inverse proposals by grounded-Laplacian residual arithmetic. The 2,614 score intervals from each endpoint are reconstructed from integer operations with fixed denominator Q², Q = 2^55, before computing 3,415,191 distinct signed comparisons with the observed 2AAS variance ranks. Stage 1 and Stage 2 archive hashes are pinned. Zero target-tied comparisons occurred for this dataset.

### Main proven output (conditional on defined graph model and archived inputs)

For every edge assignment allowed by the ±0.05 Å coordinate box (indeed for the larger set of independent optional-edge assignments), the signed model-target directional concordance C satisfies:

**33131/3415191 ≤ C ≤ 2201509/3415191**,

approximately **+0.00970107 ≤ C ≤ +0.64462251**.

The lower bound is **strictly positive**, unlike the previous two-extreme-graph interval. The lower bound is attained by the certificate's conservative case calculation at selected-contact bitmask 36; it is **not** claimed to be attainable by a real protein structure. The certified interval is **not** a statistical confidence interval, Spearman rho, or Kendall tau-b. It is a conditional combinatorial statement about the ranking sign of one explicitly defined coarse-grained GNM contact network against fixed archived NMR structure variance ranks.

Exact 4-edge intermediate case bound: **11495/3415191 ≈ +0.00336584**, using 16 cases. This shows refinement improves over the inconclusive unsplit interval without requiring all 64 cases.

### Meaning of “exact” and trust limits

1. The model here is **scalar, unweighted Gaussian Network Model effective resistance** with inclusive 10 Å contact cutoff. Directional ANM/Hessian variations within a fixed graph are **not** covered by this GNM theorem.
2. Exact rational coordinate/contact status: all arithmetic uses the dyadic rationals represented by uploaded binary64 coordinates and a rational 0.05 Å box. This is not laboratory uncertainty calibration.
3. Exact integer response proof: proposed approximate inverse matrix N is untrusted. The independent checker evaluates residual R = QI−AN exactly, with a grounded-graph spanning-tree Poincaré constant B. For each pair, lower = (Q fᵀu + uᵀr)/Q², upper = lower + B‖r‖²/Q². All 128 endpoint witnesses and 334,592 individual pair intervals are replayed.
4. The universal branch covering and Rayleigh-response implication are additionally formalized in `FoldBounds/MotionGNMAdaptiveCases.lean` building on previous Lean modules. **The complete uploaded NPZ-to-64-case arithmetic itself is not Lean-kernel re-executed**, so this is a layered proof-carrying workflow, not an end-to-end Lean certification of a protein dataset.
5. Every target statistic is **retrospective**. Selection of the six split edges is **reference-only**, but the *research question* and uncertainty radius were investigated after seeing the prior failed RNase result. Generalization requires the already sealed thioredoxin challenge and additional genuinely independent proteins.
6. Prior art: Rayleigh monotonicity, Sherman–Morrison effective-resistance updates, and exhaustive Boolean case decomposition are **established mathematics**. This project combines them into a protein-specific exact certificate; no claim that these mathematical ingredients are newly discovered.

## Reproduction

Required original archived ZIPs, unchanged:

- `FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip`: SHA256 `ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff`
- `FoldBounds_RNaseA_STAGE2_20261105_02_SEND_BACK.zip`: SHA256 `4f1ea752ada45453a44ea9128b6e7a43f4070e6824062b796a27144c340dfc6c`

Commands from the included scripts folder, with NumPy and local files available:

```bash
export PYTHONPATH=scripts
python scripts/gnm_adaptive_partition.py STAGE1.zip STAGE2.zip reproduced.json --splits 6
python scripts/gnm_adaptive_independent_replay.py STAGE1.zip STAGE2.zip reproduced.json
pytest -q scripts/test_gnm_adaptive_partition.py
```

To audit the delivered full certificate without rerunning an untrusted numerical inverse:

```bash
python scripts/gnm_adaptive_independent_replay.py STAGE1.zip STAGE2.zip research/adaptive_0050_six.json
```

Reference-only file-provenance controls use the original fixed SHA256 digests; all scripts stop on any byte or target-coordinate mismatch. The delivered certificate's independently checked SHA256 is `01409446a15a7d49ad7ed38c3142a88523608c5d3fb954d32c133d961e8bf9fc`.

## Literature context

The rank-one resistance-update formula is standard (e.g. *The resistance perturbation distance: A metric for the analysis of dynamic networks*, 2017, https://www.sciencedirect.com/science/article/pii/S0166218X17304626). See also *Evolution of Robustness in Growing Random Networks* (Entropy 2023), https://pmc.ncbi.nlm.nih.gov/articles/PMC10528878/ . No priority claim is made for contact-set partitioning or effective-resistance monotonicity.

## Follow-up implementation: low-rank witness reuse

The optional `--inverse-policy woodbury` in `gnm_adaptive_partition.py` replaces **128 separate floating-point dense inverse proposals** by **two** source-graph inverses plus low-rank updates for the same 64 cases. Every resulting *untrusted* approximate inverse is independently rechecked with the **same exact integer residual gate**. For the full 64-case RNase dataset, the low-rank approach independently replayed to the **identical exact concordance interval** `[33131/3415191, 2201509/3415191]`; its own valid certificate digest is `e6ce76fdda4575126aa5b94c5fef603f0b1e7694d1f60112f7fe159dc7636baa`. This is a **different valid witness file**: the floating proposals can differ in their least significant bits and therefore the untrusted inverse matrices need not be identical.

Performance on one local CPU session with `OPENBLAS_NUM_THREADS=1`, 124 residues and 128 graph endpoints:

| Task | Direct inverses | Low-rank reuse |
|---|---:|---:|
| Proposal production alone | 0.697 s | 0.381 s |
| Proposal **and exact integer residual verification** | 3.272 s | 3.345 s |

Thus source-only numerical proposal production is ~1.83× faster in this one run, but **there is no demonstrated end-to-end speedup at 124 residues**, because exact integer verification dominates. Larger-protein performance remains to be measured. No low-rank numerical output is trusted until the integer checks pass. Sherman–Morrison–Woodbury is established mathematics, not a new inverse law.

Five deliberate attacks on the real certificate were tested: missing case, fabricated case concordance, strongly corrupted **relevant** inverse-witness entry, changed contact identity and invalid spanning-tree constant. All five were rejected. The independent checker does not require unused numerical proposal entries to match a canonical value: a different approximation can remain valid if the *resulting certificate statements* still check. Byte-level integrity of the distributed reference certificate is separately provided by the published SHA-256 and archive manifest.

## Small exact mathematical witness of the improvement

A five-vertex unit-resistor path with edges (0,1),(1,2),(2,3),(3,4) and optional edge (0,2) has two possible networks. Let A denote the effective resistance between (0,2) and B between (0,3):

| Contact (0,2) | A = R(0,2) | B = R(0,3) | Strict A < B? |
|---|---:|---:|---|
| OFF | 2 | 3 | Yes |
| ON | 2/3 | 5/3 | Yes |

The independent marginal intervals are A ∈ [2/3,2] and B ∈ [5/3,3]. They **overlap**, so an all-on/all-off *per-score* hull cannot prove A < B. Exhaustive splitting on the ONE selected edge proves A < B for both possible networks, with exact rational Laplacian inverses. This constructive example establishes that case-level correlations can provide a strictly stronger ranking certificate than marginal response bounds. An executable exact SymPy test reproduces it; the corresponding universal case-cover principle has a separate Lean theorem. This is not itself a novel graph-theoretic discovery.