# Exact Cartesian contact-box geometry improves real 1D3Z GNM rank certificates

**Result:** An exact per-edge coordinate-box analysis eliminates three spurious optional contacts in the previous reference-only FoldBounds GNM uncertainty envelope. A complete, independently checked integer reduced-matrix computation over **4,096 graphs** then certifies **2,096/2,346 invariant** residue-flexibility rankings (plus **250 opposite-sign witness pairs**) for the public ubiquitin 1D3Z reference structure at hypothetical ±0.05 Å Cartesian uncertainty. No NMR motion labels or blind targets were used.

This is a rigorous **model-bound improvement**, not a demonstrated new protein-folding mechanism or experimentally validated fluctuation predictor.

## Scientific protocol and source

- Same previously source-verified actual PDB 1D3Z ubiquitin model 1, 69 Cα coordinates, residues 2–70, exact 0.001 Å, as stacked PRs #76/#79/#86. Public ENCORE Git blob `a438d36816e15647177e32d2eafa6aa9a0782747`.
- The same Gaussian-network model: unit-weight mandatory neighboring-chain edges, unit-weight nonlocal edges for sequence separation ≥3 and Cα distance **strictly less than 8 Å**.
- The independent Cartesian box ±0.05 Å per atomic coordinate is **hypothetical**, not a measured or calibrated uncertainty model. Graph ranking refers only to diagonal entries of the GNM pseudoinverse.
- This PR is fully additive and based on the green, still-unmerged draft #86. No changes to original RNase, sealed protein labels, Lean, NMR, other agents' branches, or main.

## Mathematical improvement: exact componentwise distance extrema

For a nominal 3D difference vector \(r=x_i-x_j\) (exact integer milliångström source coordinates) and independent Cartesian errors \(\delta_i,\delta_j\in[-\varepsilon,\varepsilon]^3\), the difference in each coordinate ranges over the **attainable** interval \([r_k-2\varepsilon,r_k+2\varepsilon]\). Therefore the squared distance has **exact** per-edge extrema

\[
D_{\min}^2=\sum_{k=1}^3 \max(0,|r_k|-2\varepsilon)^2,
\qquad
D_{\max}^2=\sum_{k=1}^3 (|r_k|+2\varepsilon)^2.
\]

The exact GNM contact rules at strict radius \(R\) are:

- **Guaranteed contact:** \(D_{\max}^2<R^2\).
- **Impossible contact:** \(D_{\min}^2\ge R^2\).
- **Optional:** \(D_{\min}^2<R^2\le D_{\max}^2\).

Unlike the old generic triangle-distance slack \(2\sqrt3\varepsilon<(7/2)\varepsilon\), these componentwise minima/maxima are **sharp for each individual residue pair**, with the exact source integer coordinates and strict inequality handled by rational arithmetic. As before, optional contacts sharing residues are correlated, so the family of independently chosen edge statuses remains an **overapproximation** of all jointly realizable perturbed 3D coordinates. The proof does not claim each optional-edge combination is geometrically feasible.

## Critical negative discovery in prior source-only split selection

Three edges were spuriously treated as individually optional under the earlier conservative 7ε/2 bound. They are **impossible to form anywhere within the stated coordinate boxes**:

| PDB residue pair | Exact minimum squared CA distance (mÅ²) | Contact threshold squared (mÅ²) |
|---|---:|---:|
| **25–52** | 64,514,030 | 64,000,000 |
| **43–51** | 64,585,089 | 64,000,000 |
| **50–59** | 64,223,700 | 64,000,000 |

In particular, PDB contact **43–51** was the first selected edge in #79's source-only greedy search, and **50–59** was used by the earlier alternate trace heuristic. These selection choices **were mathematically valid for the older, looser graph-envelope certification objective**, but cannot be interpreted as high-impact uncertainties under the stated Cartesian error boxes; neither contact can occur. This PR documents that correction explicitly instead of hiding the limitation of the earlier objective.

The new exact per-edge box classifier retains the same **152** certainly-on nonlocal contacts and reduces possibly-on nonlocal contacts from **167 to 164**, leaving **12** individually optional contacts. The resulting independent-edge graph overapproximation has **2¹²=4,096** worlds, instead of 2¹⁵=32,768.

## Complete exact integer result, not a float extrapolation

The new computation enumerates all 4,096 *reduced* network witnesses. Every graph is backed by a verified integer Woodbury update of dimension ≤12, using the single previously verified 69×69 base inverse. **Every one of 69 diagonal model-flexibility scores is calculated with exact integers** in each world; comparisons do not depend on floating-point matrix inverses, numerical tolerances, or sampling. For each of the 2,346 unordered residue pairs, the checker records both possible strict signs and any exact equalities across **all** networks, and supplies opposite-sign full-network masks when a reversal occurs.

| Method | Networks in conservative contact envelope | Strictly invariant residue pairs | Pairs with opposite-sign network witnesses |
|---|---:|---:|---:|
| Original independent diagonal intervals (#76) | 32,768 | 1,884 certified (not known complete) | Not classified |
| Improved joint covariance and 3-edge case certificate (#79) | 32,768 | 1,947 certified (not known complete) | Not classified |
| Complete exact rank-two graph certificate (#86) | 32,768 | **1,980 complete** | **366** |
| **New exact pairwise geometry + full reduced-graph audit** | **4,096** | **2,096 complete** | **250** |

This certifies **116 additional strictly invariant pairs compared with the complete old envelope**, or **212 additional pairs compared with the original diagonal certificate**, without altering the nominal GNM model or using NMR motion labels. The remaining 250 reversing pairs have exact graph witnesses; a witness graph might still be infeasible jointly under Cartesian coordinates. It would be incorrect to call the number 250 a measured frequency of physical conformational changes.

The uncertainty-width controls give identical expected counts where the older geometry already was tight:

| Cartesian ε | Optional edges | Exact graph worlds | Invariant | Reversing |
|---:|---:|---:|---:|---:|
| 0 Å | 0 | 1 | 2,346 | 0 |
| 0.01 Å | 4 | 16 | 2,268 | 78 |
| 0.025 Å | 7 | 128 | 2,220 | 126 |
| **0.05 Å** | **12** | **4,096** | **2,096** | **250** |

These are not experimental confidence intervals.

## Reproduction and trust boundary

Run from the **stacked** branch (requires #76's pinned public first-model Cα fixture, #86's reduced exact integer inverse checker):

```sh
python -m pip install -e '.[test]' sympy
python -m pytest -q tests/test_gnm_exact_coordinate_boxes.py
python experiments/real_1d3z_exact_coordinate_boxes.py --epsilon 0.05
```

Tests verify the exact per-edge strict cutoff at equality; integer-only minimum/maximum distance formula; exactly the three spuriously optional edges; nested-conservative envelope property; all 4,096 independently reduced matrix inversions and 2,346 signed pair classifications; zero/0.01/0.025 Å controls; malformed graph and budget rejection. A dedicated GitHub workflow tests Python 3.11, 3.12 and 3.13. No GPU, molecular dynamics or target label required.

**Mathematical trust:** These calculations rely on the standard Woodbury matrix identity, Laplacian inverse Loewner order, and exact integer matrix multiplication in Python; there is **no kernel-checked Lean refinement** of these new Python/Cα input paths. Pairwise geometry extrema are elementary exact box optimization, not an original mathematical theorem.

**Scientific trust:** The deposited NMR model provides real atomic coordinates but the GNM is still a simplified toy elastic model. These results prove properties of its *predicted flexibility rankings* for precisely defined coordinate errors; they neither establish real protein folding kinetics nor prove better agreement with independently measured motions. Independent prospective experimental validation is still required.
