# FoldBounds: Variable-depth, proof-carrying GNM contact certificates

**Date:** 2026-11-05  
**Status:** Implemented, independently replayed on the sealed 7RSA/2AAS historical data, Lean CI pending on new branch.  
**Scientific status:** *Retrospective mathematical/model-certification improvement, not a new blind biological success.*

## Registered biological result remains unchanged

The predeclared 7RSA (RNase A X-ray) -> 2AAS (32 deposited NMR models) weighted ANM test failed: Spearman 0.421357, short of the 0.50 threshold and below the GNM-resistance baseline of 0.511009. Deposited NMR conformers are not equilibrium or kinetic samples. No new results below alter that failure.

## The problem addressed

For a hypothetical independent ±0.05 Å Cartesian uncertainty box around 124 Cα reference coordinates, exact dyadic geometry identifies 972 contacts always present, 108 potentially switching, and 6,546 always absent at 10 Å. The certain-contact graph is connected, so the scalar GNM effective-resistance response is defined for every optional-contact assignment. Standard electrical-network Rayleigh monotonicity bounds the response of each uncertain graph between its two corresponding endpoint networks.

Previous FoldBounds PR #50 partitions the first six *reference-only selected* optional contacts into all 64 Boolean masks, then verifies two exact-rational residual witnesses per mask. On the already-revealed 2AAS target ranks, the 64-case scheme guarantees signed directional concordance **C ≥ 33131/3415191 ≈ 0.00970106796** for every contact network in the box. The 64-case certificate has 128 integer-matrix endpoint witnesses (~33.0 MB uncompressed JSON).

## New implementation: heterogeneous-depth prefix trees

A full binary prefix decision tree on the six selected contacts **does not require equal-length prefixes**. Every full six-bit contact mask has exactly one terminal prefix. At a depth-d leaf, the first d optional contacts are fixed according to that prefix; *every other* optional contact is still universally quantified. The lower endpoint graph includes all not-forced-absent contacts, and the upper endpoint includes only forced-present and certain contacts. The exact original grounded-Laplacian residual enclosures remain valid. The independent verifier checks the prefix tree against every one of 64 complete masks, so missing or overlapping assignments fail closed.

A leaf is terminal only if its target-signed lower bound satisfies the explicitly requested rational certificate threshold (or the maximum depth is reached). **The stopping criterion uses already revealed NMR labels.** The six contacts and their order are frozen from the 7RSA reference, but the tree structure is a **retrospective, target-dependent compression** and must never be represented as a prospective predictor.

**New contribution:** graph endpoint inverse witnesses are *shared* between parent and child prefixes. If one contact splits a parent branch, the OFF child's empty-remaining endpoint equals the parent's, and the ON child's all-remaining endpoint equals the parent's. Hence only two fresh endpoint graphs are required per internal split. A full binary prefix tree with L terminal leaves has L−1 internal splits and **exactly 2L endpoint graph witnesses** (no need to store redundant parent copies). For these inputs the producer maintains an exact graph-key cache and the verifier checks graph IDs are neither misused nor duplicated. The mathematical tree-shape identity is formalized in Lean.

## Independently verified results

| Retrospective required lower C | Leaves | Distinct graph witnesses | Verified actual lower bound | JSON SHA-256 |
|---|---:|---:|---:|---|
| ≥1/10000 (+0.0001), hence strictly positive | **13** | **26** | **497/3415191 ≈ +0.0001455263** | `c3be0f1c79711eba552c48c46c65c13a0c4fde7edb3805795da59d059a7103a6` |
| ≥1/200 (+0.005) | **22** | **44** | **263/50973 ≈ +0.0051595943** | `0c609166d0f9a31647169a4a4b9f538a8367477a1ef0bd8287ac10b3b39fa237` |
| ≥97/10000 (+0.0097) | **54** | **108** | **33131/3415191 ≈ +0.0097010680** | `3017cdb7526d34f3d67cac8792ebff50fc5a3dcafb9e8125549614a2afca1120` |
| Fixed-depth six selected contacts (PR #50) | 64 | 128 | 33131/3415191 ≈ +0.0097010680 | `01409446a15a7d49ad7ed38c3142a88523608c5d3fb954d32c133d961e8bf9fc` |

The 13-leaf certificate uses **79.6875% fewer endpoint witnesses and approximately 79.7% less JSON storage** than the fixed-depth certificate, though its certified positive margin is much thinner. The 54-leaf certificate retains the **same lower bound** with **15.625% fewer witnesses**. All three compressed files passed independent exact integer replay. A 13-leaf producer rerun generated the exact same SHA-256. A single local run of the 13-leaf producer took **3.83 seconds** with OPENBLAS_NUM_THREADS=1; independent replay took **3.70 seconds** on the same environment. These timings are illustrative, not hardware-general speedup claims.

### Provenance and exact-integer trust boundary

* Original reference ZIP SHA-256: `ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff`.
* Original NMR evaluation ZIP SHA-256: `4f1ea752ada45453a44ea9128b6e7a43f4070e6824062b796a27144c340dfc6c`.
* Checked 64 possible selected-contact masks; all 108 potentially switching contacts are still covered by Rayleigh sandwiching, without assuming they vary independently in any actual coordinates.
* All score numerators are exact integers with common denominator (2^55)^2; every integer witness is rejected unless it satisfies the grounded Laplacian residual bound, verified with a source-derived spanning tree constant B=262.
* The exact rational certificate checks concern **the defined unweighted GNM effective-resistance model**. They neither certify real physical motion nor validate the choice of coordinate box as an experimental distribution.
* The byte-level dataset and 124-node residual arithmetic are replayed by **independent Python**; Lean 4.28 checks only the general heterogeneous-depth covering/concordance theorem, refinement monotonicity, and a small arithmetic witness of the *reported* leaf counts. Lean does **not** replay the 26–108 full-protein matrices. Avoid calling this end-to-end Lean validation.
* Independent malicious-input tests reject missing/overlapping leaves, altered split-order metadata, fabricated leaf counts, altered tree certificates, corrupted graph witness references, oversized bogus inverses, and false global lower bounds.

## Proofs and incremental guarantees

`FoldBounds/MotionGNMPrefixTree.lean` proves: (1) every Boolean world uniquely routes to a terminal leaf of any finite, potentially unbalanced tree; (2) given valid score envelopes in every leaf, a common leaf lower-concordance bound applies to every full contact world; (3) safely intersecting score enclosures preserves truth and certified signed lower bounds cannot decrease when score intervals narrow; (4) the structural identity leafCount = splitCount + 1. The existing Lean PR #50 Rayleigh theorem provides the link from contact endpoint graphs to the realized model response. The numerical premises come from exact integer Python replay, **not** kernel loading of the archived PDB/NPZ files.

These are general mathematical results about verified Boolean case trees; decision trees, Rayleigh monotonicity, effective resistance, and interval refinement are classical ideas. No novel theorem-priority claim is made merely because the protein application uses these ingredients.

## Independent replication commands

From the ZIP's `scripts/` directory, with Python 3.11+, NumPy and pytest installed and the *original* two archived Stage 1/2 RNase experiment ZIPs on disk:

```bash
export PYTHONPATH=scripts
python scripts/gnm_prefix_tree_producer.py STAGE1.zip STAGE2.zip reproduced.json --threshold 1/10000
python scripts/gnm_prefix_tree_independent_replay.py STAGE1.zip STAGE2.zip reproduced.json
pytest -q scripts/test_gnm_prefix_tree.py
```

Optional real-source attack tests (their paths are supplied explicitly):

```bash
export FOLDBOUNDS_RNASE_STAGE1=/path/to/STAGE1.zip
export FOLDBOUNDS_RNASE_STAGE2=/path/to/STAGE2.zip
export FOLDBOUNDS_PREFIX_TREE_CERT=/path/to/research/prefix_tree_positive.json
pytest -q scripts/test_gnm_prefix_tree_real.py
```

Both the producer and the verifier stop upon SHA-256 mismatch or invalid input. Do not modify the archived original files to make a failing check pass.

## Disconfirming investigation

We also investigated whether the case partition creates much larger sets of *globally invariant, source-only score orderings* independent of the NMR labels. A 13-leaf exploratory check found 2,274,578/3,415,191 ≈ **66.60%** globally certified strict score comparisons at ±0.05 Å, only a modest change from the previously known ~66.47% unconditional ranking baseline. This did **not** demonstrate a predictive breakthrough and is not the principal result of this contribution. No prospective thioredoxin target 1XOA was fetched or evaluated.

## Next scientific gates

1. Complete and review Lean CI on the stacked draft PR, including its entire original axiom inventory. No merge until independently reviewed.
2. Extend dataset-to-Lean executable checking for **all** 26 endpoint integer matrices; the current independent Python checker is a residual trust boundary.
3. Obtain the *unseen* thioredoxin Stage 1 PyMOL ZIP, audit its source-only model freeze, and only then open its sealed Stage 2 target. Any target-informed model changes after that freeze must be explicitly labeled exploratory.