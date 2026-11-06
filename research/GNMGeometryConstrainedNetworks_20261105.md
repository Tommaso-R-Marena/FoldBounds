# FoldBounds — exact geometry-aware contact constraints

**Research date:** November 5, 2026. **Status:** Exact-rational original-reference tests independently replayed; Lean general soundness separately CI-gated. **Scope:** structural uncertainty certificates for a specified scalar GNM, not a new physical protein-folding law or a new blind experimental prediction.

## Scientific reason for the new approach

Earlier FoldBounds PRs #43, #50, #55, #58, #62, #72, and #78 bounded Gaussian Network Model (GNM) effective resistances across all **independent** on/off assignments of 108 near-threshold contacts. However, these contacts depend on **shared residue coordinates**. The original uncertainty set allowed every Cartesian coordinate of every reference Cα position to vary independently within ±0.05 Å. A contact is **on** if its squared Cα distance is ≤100 Å², and **off** otherwise. Not every abstract assignment of these contacts can come from a single set of perturbed coordinates.

For edge \(e=(i,j)\), write \(r=x_i-x_j\), \(\Delta=\delta_i-\delta_j\), with \(|\delta_{i,k}|\le\varepsilon=1/20\) Å. Then exactly

\[
 \|r+\Delta\|^2=\|r\|^2+2r\cdot\Delta+\|\Delta\|^2,
 \qquad 0\le\|\Delta\|^2\le12\varepsilon^2.
\]

Every **on** assignment therefore implies the necessary rational linear inequality \(2r\cdot\Delta\le100-\|r\|^2\). Every **off** assignment implies \(-2r\cdot\Delta\le \|r\|^2-100+12\varepsilon^2\). These are necessary, not sufficient, conditions. Any contradiction among them excludes an actual contact assignment, without ever making an invalid positive-realizability claim.

For a candidate set of 2–3 contact assignments, choose nonnegative rational weights \(\lambda_e\). The producer uses an **untrusted** two-constraint breakpoint search or SciPy LP to propose them. The independent verifier does not import SciPy or the producer. It checks, with exact rational arithmetic, that

\[
 \sum_e \lambda_e b_e
 \; <\; -\varepsilon\sum_{i,k}\left|\sum_e\lambda_e a_{e,i,k}\right|.
\]

The right side is a rigorous lower bound on the combined linear left side for *every* coordinate perturbation in the box. The strict inequality is a finite Farkas-style impossibility certificate. No numerical tolerance, Cholesky, or floating-point decision enters its acceptance check. The original binary64 coordinates are promoted to exact rational numbers; the independent verifier uses a common integer denominator.

## Results from three already frozen references

| Protein, reference-only | Residues | Potentially switching contacts | New 2-edge no-goods | Additional 3-edge no-goods | Fraction of abstract networks ruled out |
|---|---:|---:|---:|---:|---:|
| RNase A (7RSA) | 124 | 108 | 8 | 12 | **95.3990727663%** |
| CspB (previously frozen) | 67 | 52 | 7 | 2 | **87.3828887939%** |
| Crambin (previously frozen) | 46 | 36 | 14 | 1 | **95.0733184814%** |

RNase A exact source archive SHA-256: `ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff`. CspB and crambin original prior prospective archive SHA-256: `666a5d734fbd61be0f182c954474841bb552b71b3a1e9156d6c21b8c2f095f95`. For CspB and crambin **only `*_FROZEN_BEFORE_TARGET.npz` was read**; the NMR evaluation files present in the same historical archive were never opened for this analysis.

For RNase A, the 20 exact clauses induce connected components of at most **7 Boolean contacts**. Independently enumerating each component leaves at most

\[
 \mathbf{14,930,862,513,666,869,914,792,120,811,520}
 \quad\text{of}\quad 2^{108}=\mathbf{324,518,553,658,426,726,783,156,020,576,256}
\]

abstract contact assignments, hence **at least** `32010617/33554432` (~95.39907%) are geometrically impossible. The remaining networks are only *not yet excluded*, **not proven geometrically realizable**. The exact source-only JSON no-good certificate is 12,238 bytes with SHA-256 `5543f17fa7c57d956a6bf0516247ad23d6d59ba2d9c32dfea8a7527c2630c97d`. All 20 rational margins are strictly positive; the smallest is approximately 0.0020243 in the units of its normalized separation, but the checker compares exact fractions, not rounded summaries.

The combinatorial proportions are **not probabilities** for protein motion or physical conformations. Coordinate perturbation box width is hypothetical and not empirically calibrated. The no-good clauses are only a subset of all possible geometric contradictions; higher-order constraints may exclude more networks.

## Improved global exact GNM model-rank certificate

To show a computational consequence, we independently intersected seven **monotone** geometric no-goods with the prior 64-case exact-resistance and coupled-ranking certificates. For any all-off no-good \(\{e_1,\ldots,e_k\}\), at least one contact must be on. If a case's sparse graph is \(L\), then for every real coordinate perturbation the true resistance satisfies

\[
 R_{\mathrm{world}}(p)\le\max_{e\in\{e_1,\ldots,e_k\}}R_{L+e}(p).
\]

Similarly an all-on no-good guarantees at least one edge absent from a case's dense graph, yielding

\[
 R_{\mathrm{world}}(p)\ge\min_{e\in\{e_1,\ldots,e_k\}}R_{U-e}(p).
\]

Every endpoint is generated by the preexisting exact integer Sherman–Morrison update identity, starting from base graphs whose matrices satisfy `A*N=D*I`. The seven clause bounds are intersected per pair, then the old **849,112 independently certified coupled comparisons** are added only when not already covered by the new interval bound (no double counting). Every one of the 64 case rows was independently reproduced with a separate checker.

| Guaranteed signed concordance with **nominal GNM**, not experiment | Exact guaranteed lower numerator over 3,415,191 | Decimal |
|---|---:|---:|
| Original independent contact interval | 1,252,145 | 0.366640 |
| Previous coupled-rank PR #78 | 1,278,611 | 0.374389 |
| **Geometry-aware coupled certificate (this research)** | **1,302,351** | **0.381341** |

The new source-only guarantee improves by **0.0069513 signed-concordance units** over PR #78 and certifies **986,410 additional case-specific pair-order comparisons** summed across all 64 cases after removing overlap. The minimum occurs in case 55. This calculation never opens 2AAS NMR data or 1XOA. It does not change the original RNase A prospective weighted-ANM **FAIL**.

## Formal assurance and honest trust boundaries

`FoldBounds/MotionGNMGeometricNogoods.lean` formalizes the box linear-separation inequality, on/off contact linearization, nonnegative weighted halfspace soundness, and rejection of an impossible contact pattern. The Lean kernel verifies these **general statements**, without knowing an optimizer or trusting its output. The concrete deposited Cα bytes, SHA-256 checks, 20 exact λ vectors, exact forbidden-pattern count, and source-only rank result are independently checked in Python with exact rational/integer arithmetic. **The 124-residue source data have not been parsed and checked within Lean itself.**

The producer is intentionally **untrusted**; the independent replay imports neither producer nor SciPy and rejects missing/forged patterns, changed coordinates, nonpositive separation margins, false exact network counts, and modified source ZIPs. Real-data archive tests are **opt-in** in hosted CI because repository CI does not vendor the user's frozen protein archives. The portable research package includes unchanged frozen original archives, exact certificate outputs, a reproducibility script, and SHA-256 manifest.

## Replication

From the extracted release directory, with Python 3.11+, NumPy, SciPy, SymPy, and pytest:

```bash
OPENBLAS_NUM_THREADS=1 python scripts/gnm_geometry_producer.py data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip data/recomputed_nogoods.json
python scripts/gnm_geometry_independent_replay.py data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip data/geometry_nogoods.json
python scripts/gnm_geometry_crossprotein_independent_replay.py data/FoldBounds_Prospective_SEND_BACK.zip data/cspb_geometric_nogoods.json
python scripts/gnm_geometry_crossprotein_independent_replay.py data/FoldBounds_Prospective_SEND_BACK.zip data/crambin_geometric_nogoods.json
```

For the 64-case report, verify the existing PR #78 witness first:

```bash
python scripts/gnm_coupled_independent_replay.py data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip data/gnm_certain.json data/gnm_possible.json data/gnm_nominal.json data/gnm_source_exact_64cases.json data/gnm_coupled_witnesses.json
```

Then use the geometry-aware independent checker (whole run or four verified 16-case batches):

```bash
python scripts/gnm_geometry_rank_independent_replay.py data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip data data/geometry_nogoods.json data/gnm_geometry_rank_verified.json --skip-existing-audit
```

`--skip-existing-audit` is only appropriate **after** explicitly verifying the previous coupled certificate as shown. The default verifier replays prior work itself. All exact 64 rows were checked in the original analysis in four 16-case batches using the independently validated predecessor.

This is a model-relative verification improvement, not a new biological experiment, a proof of a novel universal physical law, or a license to alter the original preregistered threshold.