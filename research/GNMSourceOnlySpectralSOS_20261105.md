# FoldBounds — Exact source-only spectral SOS certification

**Date:** November 5, 2026. **Status:** Independently reproduced on the frozen 7RSA reference archive; Lean generic SOS integration gated by GitHub CI. **Scope:** scalar GNM mathematical conditioning, *not* a new biological prediction.

## 44.493102× improvement in a global residual bound

For the **972 certainly-on contacts** in the frozen RNase A reference structure (7RSA) under an independently specified ±0.05 Å Cartesian coordinate box, form the **123×123 grounded unweighted Kirchhoff/Laplacian matrix** \(A\) with anchor residue index 45. Previously, FoldBounds used a source-only spanning-tree constant \(B=262\) in its rigorous residual response enclosure. An independently checked rational electrical-flow witness improved this to \(B=72886878439/4294967296\approx16.970299\).

The new certificate improves this to the exact **spectral inequality**

\[
\forall u\in\mathbb Q^{123},\qquad
\frac{169821}{1000000}\,\|u\|_2^2\le u^\mathsf{T} A u.
\]

Consequently a rigorous response upper bound can use

\[
\boxed{B=\frac{1000000}{169821}\approx5.8885532414},
\]

**44.493102× smaller than 262**, and **2.881913× smaller than the previous checked electrical-flow bound**. This coefficient applies to **every contact-network supergraph** of the 972 certain-on contacts, because additional nonnegative springs only increase the quadratic energy. This does not rely on all 108 uncertain contacts being independently realizable, or on the target NMR ensemble.

## Fully exact integer SOS witness

The source-only producer proposes a numeric upper-triangular Cholesky matrix, then rounds each entry to integers in \(W\), scaled by \(Q=2^{30}\). The proposed matrix and the numerical Cholesky solver are **not trusted**. Write \(p=169821\) and \(q=1000000\), and recompute the **entire** residual matrix

\[
S\;=\;Q^2(qA-pI)-qW^\mathsf{T}W
\]

using arbitrary-precision integers. Its 123 diagonal dominance inequalities all passed:

\[
\forall i,\qquad S_{ii}\ge\sum_{j\ne i}|S_{ij}|.
\]

In particular, the minimum exact integer surplus \(S_{ii}-\sum_{j\ne i}|S_{ij}|\) equals **210342143431688704** (positive). This gives an explicit nonnegative sum-of-squares decomposition:

\[
S=\sum_i d_i e_i e_i^\mathsf{T}+
\sum_{i<j}|S_{ij}|\,(e_i+\mathrm{sgn}(S_{ij})e_j)
(e_i+\mathrm{sgn}(S_{ij})e_j)^\mathsf{T},
\]

where \(d_i=S_{ii}-\sum_{j\ne i}|S_{ij}|\ge0\). Thus \(S\) is positive semidefinite, \(qW^TW\) is positive semidefinite, and \(qA-pI\) is positive semidefinite. The independent checker explicitly reconstructs the diagonal and all off-diagonal entries from these square factors, not just the minimum surplus. There are 3,944 nonzero off-diagonal signed pair squares and 123 positive diagonal squares in the residual, as well as the 123 rows of the proposed Cholesky factor.

The source-only certificate contains **15,129 integer matrix entries** and is **71,264 bytes** before ZIP compression. Full exact JSON certificate SHA-256: `f2d316f81635ba01b2bbd1331957f98f2759a469b8772b1288a1e4f54871aa62`.

This result is numerically close to the reciprocal of the reference graph's smallest *floating-point diagnostic* eigenvalue (~0.1698221), but the reported conclusion relies **only** on exact integer identities and positive SOS coefficients.

## Experiment integrity

- Source archive: exactly `ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff` (original Stage 1 frozen 7RSA).
- All ZIP-internal file SHA-256 checks and exact binary64 coordinate-to-contact classifications are independently reconstructed.
- 124 reference residues; 972 certainly-on contacts and 108 optional contacts.
- No target data in the producer or checker. **The 2AAS RNase A test outcome remains a registered failure** for the predeclared weighted ANM predictor. The thioredoxin 1XOA target remains sealed.
- Exact integer replay uses no Cholesky solver, spectral eigensolver, or NumPy integer arithmetic in the trusted proof step. Numerical floating-point routines appear **only** in the untrusted producer.
- Hermetic tests: **11 passed**. Separate original-source positive and tamper tests: **7 passed**. Malicious source hashes, contact lists, matrix entries, nonpositive bounds, false reported surpluses, and malformed matrices are rejected.

## Formal proof status and trust boundary

`FoldBounds/MotionGNMSpectralSOS.lean` supplies a general executable rational **SOS identity checker**: a proposed matrix \(A\), coefficient \(\tau>0\), vectors \(b_e\), and nonnegative weights \(w_e\) are accepted only if every entry satisfies \(A=\tau I+\sum_e w_eb_eb_e^\mathsf{T}\). Existing Lean-checked Hessian and exact solvability results then prove a unique exact response and a residual interval with coefficient \(1/\tau\), without assumed coercivity or trusted inverse. The finite-numeric protein certificate is checked by a **separate** exact-integer Python replay, which derives the corresponding nonnegative SOS vectors from the Cholesky-plus-diagonally-dominant residual. Full original NPZ/ZIP parsing and all 15,129 real protein matrix values are **not yet executed by the Lean kernel**; do not conflate the two verification layers.

There is no new claim about protein dynamics, biological agreement, a kinetic law, or blind generalization. This is an improved **mathematical assurance mechanism** for an already-defined GNM model. It does not change the frozen 2AAS findings or the older ANM predictor.

## Reproduction

With the original, unchanged user-provided `FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip`:

```bash
OPENBLAS_NUM_THREADS=1 python scripts/gnm_spectral_sos_producer.py STAGE1.zip spectral.json
OPENBLAS_NUM_THREADS=1 python scripts/gnm_spectral_sos_independent_replay.py STAGE1.zip spectral.json
python -m pytest -q scripts/test_gnm_spectral_sos.py
FOLDBOUNDS_RNASE_STAGE1=STAGE1.zip FOLDBOUNDS_SPECTRAL_CERT=spectral.json python -m pytest -q scripts/test_gnm_spectral_real.py
```

The signed original Stage 1 archive and exact certificate SHA digests must match. If they do not, the verifier rejects rather than adjusting the acceptance threshold.

## Source-only precision ablation: concrete computational impact

To distinguish a merely smaller algebraic coefficient from an improvement in practical certificate strength, I regenerated **the same** grounded-GNM approximate inverse witness with several dyadic denominator sizes, then checked every residual with the existing exact-integer checker. I compared the original tree constant \(B=262\) with the newly proven \(B=1000000/169821\) on all 2,614 frozen **reference-only** pair scores. The numerical proposal is identical at each precision for both bounds, and neither the test nor the producer opens the NMR target.

| Dyadic witness denominator | Certified adjacent score orders with B=262 | With exact spectral SOS B | Maximum certified response interval width (old → spectral) |
|---|---:|---:|---|
| \(2^{16}\) | **102** | **2,138** | 0.000529128 → 0.0000118924 |
| \(2^{18}\) | 1,586 | 2,559 | 0.0000315493 → 0.0000007091 |
| \(2^{20}\) | 2,506 | 2,586 | 0.00000209526 → 0.00000004709 |
| \(2^{24}\) | 2,589 | 2,590 | 8.7794×10⁻⁹ → 1.9732×10⁻¹⁰ |
| \(2^{55}\) | 2,590 | 2,590 | 3.8886×10⁻²⁷ → 8.7398×10⁻²⁹ |

There are 2,613 adjacent comparisons across the 2,614 pairs; some nominal or near-tied comparisons remain uncertifiable with either coefficient. The ordering in this table is the lower-bound ordering of the corresponding exact intervals at each precision, **not an NMR-observed biological ranking**. In particular, at \(2^{16}\) the exact spectral certificate establishes 2,138 adjacent GNM score orders compared with 102 under the same old tree bound. This is a substantial source-only computational assurance gain, not improved biological prediction accuracy. Full generated data are in `research/source_only_precision_ablation.json`; rerun using `scripts/gnm_precision_ablation.py` with the sealed Stage1 archive and the pre-existing FoldBounds GNM checker modules.