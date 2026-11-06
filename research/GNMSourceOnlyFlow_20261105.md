# FoldBounds: Source-only exact electrical-flow coercivity (GNM)

**Date:** 5 November 2026. **Status:** Implemented; independent integer replay of the 7RSA reference passed. Full Lean proof compilation is independently gated in GitHub CI. No changes to `main`; no access to thioredoxin target 1XOA or NMR outcome during model development.

## The gap closed

A previously accepted FoldBounds GNM certificate required the coercivity premise

\[\|u\|_2^2 \le B\,u^\top L u\]

for every grounded vector on the *certain-on* contact graph. An external BFS tree proof supplies `B = 262`, but the Lean GNM theorem had an **assumed** coercivity premise, with the 124-vertex tree argument verified separately in Python.

Here the *untrusted* source-only proposal includes an exact rational **flow matrix** `P` with one row per grounded vertex and one column per certain-on contact. For the grounded edge-incidence matrix `D`, the entire 119,556-entry certificate is checked using **integer** arithmetic:

\[ P D = I. \]

This makes all grounded vertices reconstructible as sums of edge differences. By finite Cauchy–Schwarz,

\[
\sum_v u_v^2=\sum_v(P_vDu)^2
\le \left(\sum_v\sum_e P_{ve}^2\right)\left(\sum_e(Du)_e^2\right).
\]

Thus **the effective error constant is certified to be**

\[B_{\mathrm{flow}} = \sum_{v,e}P_{ve}^2.\]

On every supergraph formed by adding any subset of the 108 uncertain contacts, spring-energy monotonicity gives

\[\|u\|_2^2\le B_{\rm flow}\sum_{e\in E_{\rm certain}}(Du)_e^2\le B_{\rm flow}\sum_{e\in E_{\rm actual}}(Du)_e^2.\]

The proof requires no independence assumption on contact switches, no target/NMR data, no trusted eigensolver and no calculation of an exact rational inverse.

## Real protein result (7RSA RNase A, reference-only)

- Exact original Stage 1 reference archive SHA256: `ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff`.
- Reference residues: 124; `±0.05 Å` independent Cartesian coordinate box; certainly-on edges: 972; uncertain edges: 108.
- Ground vertex: 45 (zero-based); grounded dimension: 123; flow columns: 972.
- Dyadic denominator `Q=65536`; certificate has **119,556 integer flow entries**.
- Previous source BFS constant: `262`.
- **Exact new coefficient:** `72886878439/4294967296 = 16.9702988208737...`.
- **15.4387381604× tighter uniform residual upper-bound coefficient.**
- Independent source-only replay of all graph edges and every conservation equation passed. SHA256 of exact 477,382-byte JSON certificate: `86cecf5afa9ad5d10ba9808ba6b4e03395ef5a98115aff4dfb31c8d17b9dc059`.
- A separate high-precision numerical comparison gives `trace(L^{-1}) ≈ 16.9702892223`; this suggests the exact flow approaches the best possible *sum-of-row-norms* coefficient, **but the numerical trace is diagnostic only and not a certified lower bound**.

## How the source-only witness is made

1. Parse the original SHA-pinned 7RSA Stage 1 archive and its file manifest; never read 2AAS, 1XOA or any NMR result.
2. Recompute exact binary64-as-dyadic coordinate box classifications and the certain-contact graph.
3. Pick a deterministic source-only root and shortest-path spanning tree.
4. Numerically propose `L^{-1} Dᵀ`, round every entry to the dyadic denominator `Q` (untrusted).
5. Calculate the exact integer conservation defect `R = Q I − W D`.
6. Repair the proposal as `W′ = W + R T`, where `T` is a source-only integer tree-path flow with `T D = I`.
7. **Independently check** `W′D = QI`, compute `∑W′²` using arbitrary-precision integers, and verify the exact source and contact SHA metadata. The proposal and repair are not accepted merely because the producer says they work.

## Formal and numerical scope

The accompanying Lean module proves the generic executable flow identity gate, finite Cauchy–Schwarz coercivity, and a uniquely existing exact GNM response with rigorous residual bounds (via existing `MotionExactSolvability` and `MotionGNMExact` theorems). The general theorem is parameterized by finite rational matrices. Its **full-protein data input is not parsed by the Lean kernel**, and is independently checked by a separate Python program. The graph supergraph statement relies on the already Lean-checked positive-edge monotonicity theorem.

All reported quantities are about the idealized, unweighted **scalar GNM**, not experimentally verified dynamics, equilibrium sampling, or a new physical law. The original weighted-ANM blind RNase trial **failed**, and this certificate does not change that. No new prediction has been evaluated on held-out thioredoxin data.

## Independent replication

```bash
OPENBLAS_NUM_THREADS=1 python scripts/gnm_flow_cert_producer.py FROZEN_STAGE1.zip flow.json --q-bits 16
OPENBLAS_NUM_THREADS=1 python scripts/gnm_flow_independent_replay.py FROZEN_STAGE1.zip flow.json
python -m pytest -q scripts/test_gnm_flow_cert.py
FOLDBOUNDS_RNASE_STAGE1=FROZEN_STAGE1.zip FOLDBOUNDS_FLOW_CERT=flow.json python -m pytest -q scripts/test_gnm_flow_real.py
```

The independent checker does not import the producer or use floating-point linear algebra. Small-graph exhaustive rational tests and six real-archive tamper/replay tests passed locally. Any hash, conservation, graph or claimed coefficient mismatch fails closed.

## Further scientific opportunities

The improved single shared coercivity constant can replace `B=262` across **every uncertain contact case** in the existing 13-leaf/26-graph tree, potentially allowing smaller dyadic response proposals or tighter rank intervals. Such downstream performance improvements require explicit regenerated certificates and independent verification; the new source-only `B` alone does not prove a better biological predictor or faster end-to-end runtime.