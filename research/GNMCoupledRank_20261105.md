# FoldBounds: source-only coupled GNM rank certificates

## Scope and scientific status

A **new rigorous certification technique** for the previously defined scalar Gaussian Network Model (GNM), not a newly discovered physical law, not a newly validated biological predictor, and not a success on the RNase A registered weighted-ANM blind test. The dataset is the original frozen PDB 7RSA **reference-only** archive. **No 2AAS target NMR data or sealed thioredoxin 1XOA target has been accessed**.

The previous 64-case certificate treats uncertainties in effective resistances independently. For each case it has two exact, grounded graph Laplacian inverses, with sparse and dense graph endpoints. Writing their inverse matrices as \(G_-\) and \(G_+\), we know for every allowed graph \(G\) that \(G_+\preceq G^{-1}\preceq G_-\) in positive semidefinite (PSD) order. Let \(M=G_- -G_+\succeq0\) and \(K=G_- -G^{-1}\), so \(0\preceq K\preceq M\).

For residue-pair forces \(p\) and \(q\), put \(u=p+q\), \(v=q-p\). Then the difference of resistance responses is \(D=q^TG^{-1}q-p^TG^{-1}p=u^TG^{-1}v\). Because \(K\) and \(M-K\) are PSD, the **paired cone contraction** gives

\[
|u^T(2K-M)v|^2 \leq (u^TMu)(v^TMv).
\]

This is a standard consequence of PSD order and Cauchy-Schwarz/positive-quadratic discriminants, applied here to the **shared contact uncertainty**. It is not valid to assume each resistance can vary independently. Define endpoint decreases \(d_p=p^TMp\), \(d_q=q^TMq\), and transfer \(c=p^TMq\). If \(s\) is the sum of the sparse- and dense-endpoint values of \(D\), the exact sufficient condition

\[
\boxed{s>0\quad\mathrm{and}\quad s^2+4c^2>(d_p+d_q)^2}
\]

certifies \(D>0\) for **every** intermediate contact network. All terms are obtained from exact integer scaled inverse witnesses. The Lean proof establishes the conditional rational cone/discriminant implication; the 123×123 real protein inputs are replayed by a separate exact-integer checker, not parsed in Lean.

The new certifier starts with the **already independently verified** 64-case source-only model-concordance result from draft PR #72 and adds disjoint witnesses for pair orders that independent intervals could not prove. A floating-point prefilter merely proposes promising comparisons; **every emitted comparison is rechecked with integers**, using a common positive denominator \(D_\mathrm{sparse}D_\mathrm{dense}\). The independent checker never imports the producer, uses no floating-point numerics, revalidates all 128 exact inverse identities, reconstructs the old exact Fenwick counts, and verifies each newly listed strict cone inequality and its non-overlap with the old count.

## Results: original frozen 7RSA source, 10 Å contacts, ±0.05 Å coordinate boxes

| Metric | Previous independent intervals | New coupled cone |
|---|---:|---:|
| Model-internal universal signed rank-concordance bound | 1,252,145 / 3,415,191 ≈ 0.366640 | **1,278,611 / 3,415,191 ≈ 0.374389** |
| Additional verified source-only comparisons across all 64 cases | — | **849,112** |
| Extra proven comparisons in worst-case (mask 55) | — | **13,233** |
| Distinct endpoint exact-integer graph certificates | 128 | 128 (reused; no new approximate inverse) |

The 64 cases cover every assignment of six chosen source-only contacts; the other 102 optional contacts are bounded by the sparse and dense graph endpoints. Together this covers a superset of all networks realizable within the stated coordinate boxes. The 64-case statement is **not** a physical uncertainty probability, an NMR correlation, a confidence interval, or prospective validation.

The full witness JSON (~9.5 MB) is distributed separately because it contains 849,112 exact integer-verified index pairs, not data that should be silently embedded in a Lean proof. SHA-256: `830801857b4e6255090bc3903405c9931bc80cc20cd625e9a295a87a160e821b`. It is replayable against the immutable Stage 1, the three pinned exact integer graph inverse certificates, and the older 64-case summary report.

## Research interpretation and limits

The improvement is real **within the GNM mathematical model**: the checker establishes previously unknown stable pair rankings through shared-contact correlations. It is quantitatively modest in aggregate (+0.007749 signed concordance), but a reusable technique for uncertainty propagation, symbolic graph mechanics, and robust electrical networks.

The 7RSA structure and nominal GNM ranking are in the same source; this is a **source-only model-self-consistency exercise**, not prediction of unseen conformational ensembles. In particular, the RNase A weighted ANM preregistered test failed, the previously exposed CspB/crambin data cannot provide a fresh blind test, and the 1XOA thioredoxin target remains sealed pending a reference-only Stage 1 archive.

The generic Lean theorem assumes two nonnegative quadratic forms for the PSD inverse-difference pieces, as justified by the graph Rayleigh theorem. While original protein data and 128 exact integer inverses are independently replayed, **the original NPZ/JSON bytes are not consumed inside Lean**. A future development should close that last trusted-decoding/byte-boundary gap, and separate predictive validity from theorem validity.

## Reproduction

With original, unchanged archive/cert files in `data/`, execute:

```bash
OPENBLAS_NUM_THREADS=1 python scripts/gnm_coupled_producer.py \
 data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip \
 data/gnm_certain.json data/gnm_possible.json data/gnm_nominal.json \
 data/gnm_source_exact_64cases.json data/gnm_coupled_witnesses.json

OPENBLAS_NUM_THREADS=1 python scripts/gnm_coupled_independent_replay.py \
 data/FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip \
 data/gnm_certain.json data/gnm_possible.json data/gnm_nominal.json \
 data/gnm_source_exact_64cases.json data/gnm_coupled_witnesses.json
```

No network access or NMR structures are required. The SHA256SUMS manifest in the release ZIP authenticates every file.