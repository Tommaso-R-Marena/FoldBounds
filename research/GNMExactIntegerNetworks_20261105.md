# FoldBounds — exact integer GNM contact-network certificates

**Research date:** November 5, 2026. **Status:** Local numerical and independent integer replay complete. Lean integration/CI status is tracked separately on draft PR (do not assume a green build until it has finished). This is a mathematical verification advance in a *specified scalar contact-network model* — it is not a new biological motion law or a positive blind validation.

## 1. Exact results, no floating inversion or residual enclosures in the verifier

Given an integer grounded graph Laplacian \(A\) with dimension 123, a completely **untrusted proposal** supplies an integer matrix \(N\) and strictly positive integer \(D\). The separate checker recomputes the integer identity

\[
A N = D I,\qquad D>0.
\]

Every one of the 15,129 entries is checked with arbitrary-precision integer arithmetic; the inverse is exactly \(A^{-1}=N/D\). Every frozen eligible pair \(i,j\) has **exact rational** effective resistance

\[
R_{ij}=\frac{N_{ii}+N_{jj}-2N_{ij}}{D},
\]

where rows/columns for the grounded vertex are interpreted as zero. There is **no numerical residual interval** in these final score certificates and no dependency on the producer's float inversion. The producer uses SymPy's fraction-free \`DomainMatrix.inv_den()\` to propose \(N,D\); the checker does not import SymPy or the producer, and rejects malformed input, erroneous contact lists, altered source ZIP content, mismatched claims, or one incorrect integer inverse entry.

| Frozen 7RSA reference graph | Contacts | Eligible pair scores | Distinct exact values | Exact ties | Certificate SHA-256 |
|---|---:|---:|---:|---:|---|
| Original nominal 10 Å graph | 1,020 | 2,614 | **2,614** | **0** | `6b413040176838eaf90cac7ac907ac95bbedff2de4197fb351b4348e62cd0706` |
| Definitely-on graph, coordinate box ±0.05 Å | 972 | 2,614 | 2,591 | 23 | `8b0a26737dbcbbe55052148459585da05d4798ba1bc7970c87e0ef8fe7b806c5` |
| Possibly-on graph, coordinate box ±0.05 Å | 1,080 | 2,614 | **2,614** | **0** | `4ca5a4046f89b1a69e7e17d34b89033d344182cf44e5cc602e47891d3a8633e1` |

The **nominal** original model graph has 2,613 strictly ordered adjacent comparisons out of 2,613. All exact scores agree with the original frozen floating GNM predictor to maximum absolute difference \(3.33\times10^{-16}\); the complete rank permutations match. The separate reference-only producer generated the 123×123 nominal integer inverse in approximately **0.86 s** on this environment. The exact solver's performance is not a general scaling guarantee.

The original blind RNase A **weighted ANM** acceptance test remains a failure (Spearman +0.421357, versus GNM +0.511009). Nothing here retroactively satisfies its registration.

## 2. Exact contact changes, with no matrix reinversion

For an edge incidence \(b=e_i-e_j\), start from a *verified* scaled inverse \(A N=D I\). Set \(v=N b\) and \(\alpha=b^T N b\). Adding contact \(+bb^T\) yields the exact integer update

\[
N'=(D+\alpha)N-vv^T,\qquad D'=D(D+\alpha).
\]

Removing a contact gives the \(-bb^T\) version

\[
N'=(D-\alpha)N+vv^T,\qquad D'=D(D-\alpha),
\]

provided \(D-\alpha>0\). A common integer gcd is removed after each step; every requested destination network is checked **again** by directly recomputing \(A'N'=D'I\), so a flawed proposer cannot create a valid false certificate.

On the original reference coordinates, 48 of the 108 potentially switching contacts are nominally present. Starting from the exact 972-edge definitely-on inverse, a source-only run added those **48 edges** algebraically, completed in about **1.43 seconds**, and reproduced all 15,129 entries of the independent exact nominal inverse after denominator normalization. This is a constructive test of exact rank-one network updates at real-protein scale.

The rank-one identity is established algebra; the contribution is integrating exact source, updater and checker into FoldBounds, not a novel matrix theorem.

## 3. Exhaustive six-contact cases over all 108 uncertain contacts

We retained the same **six reference-only optional contacts** used in FoldBounds' prior case-partition experiment:

\((51,115),(1,9),(6,32),(43,121),(1,7),(52,76)\), using zero-based Cα indices.

Every possible contact network consistent with the uncertainty box belongs to exactly one of 64 selected-contact cases. In a given case, all remaining 102 optional contacts are bounded by two endpoint graphs; adding an edge can only **decrease effective resistance** (Rayleigh monotonicity). The verified exact scaled inverses of the all-on and all-off reference graphs, together with exact rank-one additions/removals, generate the **128 exact endpoints**. Every endpoint has its full integer inverse identity checked. No floating inverse or residual interval is used at these endpoints.

All 3,415,191 comparisons between the 2,614 reference pairs were assessed against their **exact nominal GNM ordering** using exact \(\mathbb Q\) comparisons, for each case. The producer and the independent verifier use different ranking algorithms (bisect insertion versus an offline Fenwick sweep). Both obtain the exact source-only guaranteed signed pair-order concordance

\[
\boxed{\frac{1\,252\,145}{3\,415\,191}\le C_{\rm nominal}\le1.}
\]

This applies to **every contact state in all 64 cases** (including a conservative superset of coordinate-realizable graphs), conditioned on the exact specified resistance model and ±0.05 Å coordinate boxes. The universal lower bound is approximately **+0.366640**. It is a worst-case **model-internal rank stability** assertion, **not** a correlation with NMR observations, prediction accuracy, physical protein dynamics, statistical confidence interval, or independent prospective success.

The tiny **12,058-byte** case-summary JSON (SHA256 `2036f073f6f12ab30bf74370eb8ba237ac9fe8226f3357a1ce204b428ca990b0`) is sufficient when paired with the three full base inverse witnesses and their original Stage 1 archive: all 128 endpoint matrices are recomputed on demand, and the independent verifier confirms their exact identities. The full 64-case producer completed in approximately **17.1 seconds**; the separate replay in approximately **21.6 seconds** locally.

## 4. Lean kernel boundary and limitations

`FoldBounds/MotionGNMIntegerInverse.lean` introduces an executable integer acceptance gate and proves that, if the full scaled inverse is accepted, the rational matrix is genuinely a right inverse and supplies an exact rational solution to **every** force vector. Its correctness builds on the existing kernel-checked finite inverse theorem in `MotionVerifiedInverse`. A separate example checks a chain graph, an added triangle contact, false inverse, and zero denominator.

The generic Lean theorem and source audit **are not a claim that the real 7RSA NPZ ZIP, exact 123×123 integer matrix, or 64-case JSON are parsed and evaluated within Lean**. Those original real-protein bytes are independently checked by the Python integer verifier. The existing FoldBounds graph monotonicity and case partition soundness theorems apply to the already checked model intervals, conditional on faithful data translation.

**Scientific boundaries:** This is scalar GNM's *model-defined effective resistance*, not a proven protein-motion law. The uncertainty radius is stipulated, not empirically calibrated. The deposited 2AAS NMR ensemble cannot be assumed to be an equilibrium/kinetic trajectory. The previously revealed RNase A and CspB/crambin results cannot validate a new predictor. **Thioredoxin 1XOA remains unopened** pending a returned, independently checked Stage 1 freeze.

## 5. Reproduction

The release ZIP contains the scripts, SHA256-pinned original Stage 1 archive, three exact base witnesses, the 64-case case summary, regression tests, and SHA256SUMS manifest. With Python 3.11+, `numpy`, `sympy`, and `pytest`, from the extracted directory:

```bash
OPENBLAS_NUM_THREADS=1 python scripts/gnm_exact_inverse_independent_replay.py \
  data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip data/gnm_nominal.json

OPENBLAS_NUM_THREADS=1 python scripts/gnm_exact_case_independent_replay.py \
  data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip \
  data/gnm_certain.json data/gnm_possible.json data/gnm_nominal.json \
  data/gnm_source_exact_64cases.json

OPENBLAS_NUM_THREADS=1 FOLDBOUNDS_RNASE_STAGE1=data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip \
  FOLDBOUNDS_EXACT_CERTAIN=data/gnm_certain.json \
  FOLDBOUNDS_EXACT_POSSIBLE=data/gnm_possible.json \
  FOLDBOUNDS_EXACT_NOMINAL=data/gnm_nominal.json \
  FOLDBOUNDS_EXACT_64CASES=data/gnm_source_exact_64cases.json \
  python -m pytest -q scripts/test_gnm_exact_inverse.py scripts/test_gnm_exact_cases.py
```

These commands access **only the reference-side source archive**. The ZIP itself is the immutable frozen-protein input; the separate scripts do not silently retrieve target PDBs or use a network connection. All source digests and model acceptance checks are fail-closed.

## Frozen CspB and crambin portability — completed source-only audit

The same exact integer matrix certificate was independently tested on the earlier, already revealed CspB and crambin reference structures. Input: the unchanged `FoldBounds_Prospective_SEND_BACK.zip` (SHA-256 `666a5d734fbd61be0f182c954474841bb552b71b3a1e9156d6c21b8c2f095f95`). The verifier authenticates the complete ZIP by its SHA-256 but opens **only** the two `FROZEN_BEFORE_TARGET.npz` members, never target NMR evaluation members. Each frozen reference also has its own pinned digest. The numerical inverse is proposed with SymPy; acceptance independently checks every integer entry of `A*N=D*I`.

| Reference | CA residues | Exact 10 Å contacts | Frozen pairs | Exact tie comparisons | Strict order violations vs frozen float GNM | Max abs float discrepancy |
|---|---:|---:|---:|---:|---:|---:|
| CspB | 67 | 543 | 1,170 | **9** | **0 / 683,856** | 4.996e-16 |
| Crambin | 46 | 334 | 547 | **0** | **0 / 149,331** | 4.441e-16 |

**Important correction:** CspB has nine exact tied effective resistance pairs. Its frozen binary64 approximations can break those ties differently. Therefore a naive equality check between total sorted index arrays is not valid in CspB. The audit checks every **mathematically strict** pairwise ordering with a floating tolerance of `1e-12` and counts exact ties separately. This is not an integer-inverse failure or a new biological finding.

Reproduce using the original signed prospective ZIP:

```bash
OPENBLAS_NUM_THREADS=1 python scripts/gnm_portability_audit.py /path/to/FoldBounds_Prospective_SEND_BACK.zip
FOLDBOUNDS_PROSPECTIVE_ZIP=/path/to/FoldBounds_Prospective_SEND_BACK.zip python -m pytest -q scripts/test_gnm_portability.py
```

The hermetic GitHub CI test covers integer graph construction, exact threshold handling, and forged-source rejection without downloading or opening private experiment files. The real-reference test runs only with an explicitly supplied frozen ZIP. No NMR observations were used.

### Overall trust boundary
These tests establish mathematical correctness of the declared scalar GNM graph computation, **not** accuracy of biological motion predictions. The generic scaled-inverse acceptance theorem is checked by the Lean kernel; the complete real-protein data are checked by independent Python unbounded-integer arithmetic rather than automatically parsed by Lean. The RNase A original registered blind result remains a failure; the prospective thioredoxin target remains unopened.
