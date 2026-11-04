# FoldBounds structural-uncertainty response pilot (2026-11-04)

## Purpose and method

The 1,717 archived weighted-ANM model scores on crystal reference coordinates for
CspB (1CSP) and crambin (1CRN) had previously been independently checked
numerically at the fixed structures. This study goes further: every Cartesian
C-alpha coordinate is allowed to vary independently in ±epsilon, including
discrete changes in the 10 Å contact network. It is an **uncertainty analysis
of a chosen physical model**, not a measurement of physical protein dynamics.

For each candidate residue-pair score q(x), an untrusted fixed rational solve
proposal v gives
  L(x) = 2 z(x)^T v - v^T H(x) v,
  R(x) = ||z(x) - H(x) v||²,
  D(x) = ||x_a - x_b||².
For H(x) >= tau I (tau=1/5) and exact uniform interval enclosures
L in [Llo,Lhi], R<=Rhi, 0<Dlo<=D(x)<=Dhi, every allowed x obeys
  max(0,Llo)/Dhi <= q(x) <= (Lhi + Rhi/tau)/Dlo.
The contact graph may change; each world must have a complete, nonduplicated
edge enumeration. The Lean modules MotionRobustResponse.lean and
MotionRobustOrdering.lean prove the uniform enclosure and strict ordering
when their envelope premises are supplied. These are generic kernel-level
theorems, not dataset-specific proofs of the Python-produced intervals.

The independent numerical implementation uses outward-rounded signed integers
at denominator 2^96; no floating-point arithmetic is used in model
certification. All eligible reference contacts are regenerated for each
epsilon and included as definitely on, definitely off, or optionally on/off.
An alternative certificate splits 2^k discrete states for the k ambiguous
edges and takes their union.

## All-pair results

| Epsilon (Å, per Cartesian coordinate) | CspB guaranteed pairwise orders | Crambin guaranteed pairwise orders | CspB guaranteed top 10% | Crambin guaranteed top 10% |
|---|---:|---:|---:|---:|
| 0.0005 | 98.44% | 98.09% | 115/117 | 50/55 |
| 0.001 | 95.28% | 96.06% | 110/117 | 46/55 |
| 0.005 | 74.36% | 63.85% | 79/117 | 13/55 |
| 0.010 | 43.81% | 29.22% | 40/117 | 2/55 |
| 0.050 | 0% | 0% | 0/117 | 0/55 |

For epsilon=0.005, 3 CspB contacts and 4 crambin contacts are ambiguous.
Exact 2^k contact-state splitting improves certification modestly:
  CspB: 512,435 / 683,865 = 74.93%; top-decile 80/117.
  Crambin: 97,789 / 149,331 = 65.48%; top-decile 14/55.

The stronger 0.005 Å results are **guarantees about ordering of the exact
defined model's response scores across all allowed coordinates**. They
do not assert a similar level of biological prediction accuracy. At 0.05 Å
the conservative bounds correctly abstain from all pairwise comparisons.

## Validation and remaining trust boundaries

* Every archived fixed-structure exact rational score was enclosed for all
  five epsilon levels.
* Twenty cases (CspB and crambin, contact-stable and contact-switching)
  were independently recomputed using fractions.Fraction and agreed with
  the outward integer bounds.
* The numeric uncertainty implementation and output JSON records are
  distributed separately in the conversation's reproducibility archive
  FoldBounds_Structural_Uncertainty_Breakthrough_20261104.zip.
* Canonical byte-to-rational decoding and replay of each uncertainty
  interval within Lean are **not yet completed**.
* These coordinate boxes are not correlated structural ensembles.
  No mutation, binding free-energy, kinetics, or time-resolved dynamics
  claim follows from the model score.
* Exploratory comparison with earlier revealed NMR ensembles found that
  selecting the most robust predicted pair scores does **not** improve
  empirical top-decile precision. The original strict prospective
  experimental criteria were not achieved and have not been relaxed.

## Independent follow-up

A strictly reference-only RNase A challenge is prepared using RCSB 7RSA
(X-ray reference) with RCSB 2AAS solution NMR ensemble held sealed.
The 7RSA-only PyMOL freeze records predicted scores, benchmarks and a
deterministic 256-pair uncertainty sub-study; SHA-256 hashes must be
returned before evaluating target 2AAS. Protocol constants remain fixed.
No new target results have been examined.
