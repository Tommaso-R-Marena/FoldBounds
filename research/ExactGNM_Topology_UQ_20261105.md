# FoldBounds — exact GNM and structural-topology uncertainty research

**November 5, 2026 · Draft study, independent mathematical/model audit.**

## Result in brief

We implemented an exact-dyadic-input, integer-residual effective-resistance
certificate for the **standard Gaussian Network Model (GNM)** and extended it
using **Rayleigh monotonicity of effective resistance** to bound changes in
predicted scores over **all** coordinate-box-induced contact networks. This
extends the FoldBounds research prototypes to cover the *GNM baseline that
outperformed weighted ANM* in the frozen RNase A 7RSA→2AAS study.

GNM and the Gaussian/anisotropic interpolation model are longstanding physical
approximations; the work does **not** introduce a new force law and does not
validate conformer-ensemble variability as actual protein dynamics. Its useful
engineering contribution is exact, auditable, source-hash-bound score and
order certificates including contact insertions/deletions.

## Methods

Let the undirected, unweighted protein Cα contact graph contain all pairs with
exact squared separation 0 < d² ≤ (10 Å)². Treat each archived binary64
coordinate as its **exact dyadic rational** value. This freezes the meaning of
the numerical model but does not bound uncertainty in the underlying physical
coordinate measurement.

For arbitrary anchor vertex a, delete row and column a from graph Laplacian L.
The grounded Laplacian A is invertible for a connected graph, and each GNM
response is the scalar resistance

    R(i,j) = (e_i - e_j)' A^{-1} (e_i - e_j),

where anchored force components vanish. This equals the standard ungrounded
pseudoinverse resistance and is independent of the chosen anchor.

Select a BFS spanning tree rooted at the vertex minimizing the sum of
distances from all vertices (reference-only), with tree path lengths l_i.
Set **B=Σ_i l_i**. For any grounded displacement u, telescoping along the
tree paths and Cauchy–Schwarz give

    ||u||² ≤ B Σ_(edge in tree) (u_s-u_t)² ≤ B u' A u.

Generate any untrusted integer approximate inverse matrix N and integer
scale Q (here Q=2^55). For force b=e_i-e_j, let u=N*b and r=Q*b−A*u, all
as arbitrary-precision integers. The exact GNM response satisfies:

    (Q*b'u + u'r) / Q² ≤ R(i,j)
      ≤ (Q*b'u + u'r + B*r'r) / Q².

No floating-point solve is a premise of this certificate; the inverse is
only a proposal. The validated algorithm independently rebuilds the graph
using exact rationals, confirms connectedness via BFS and recalculates every
integer product. The independent verifier uses a **different** integer
coordinate decoder, independently checks tree paths, calculates spring
energy directly, and reproduces interval numerators and all certified ranks.

### Changing contact topology without 2^k enumeration

For every coordinate in a box x_{i,k} ± ε independently, each pair difference
has component-wise absolute range

    [max(0, |x_i−x_j|_k − 2ε), |x_i−x_j|_k + 2ε].

The sum of squared lower endpoints and sum of squared upper endpoints classify
all pair contacts as always-on, always-off, or optional. Every admissible
physical graph lies between an **all-on** and **all-off** network.
By electrical-network Rayleigh monotonicity, adding a contact can only reduce
R(i,j). Hence the lower score for the all-on graph and the upper score for the
all-off graph enclose **every admissible graph**, even when arbitrarily many
contacts may switch. Both extreme graph responses are certified by the above
integer residual method. No independence of contact states is assumed; all
optional states are included even if joint geometric realization is impossible.

For the fully disconnected all-off graph, the checker abstains. For possible
coordinate collisions, it abstains.

## RNase A (7RSA reference; 2AAS previously revealed)

Immutable input Stage 1 archive SHA256:
`ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff`.

- 124 reference Cα residues; 1020 contacts; 2614 eligible pairs.
- Optimal anchor index **44**, tree distance sum **B=254**.
- All **2614** GNM exact score intervals enclosed; **2613/2613** adjacent
  rankings strictly separated, independently reproduced.
- Maximum absolute discrepancy of certificate midpoint vs frozen 7RSA GNM
  scores: **3.331e−16**, a floating-point comparison diagnostic.

| Independent Cartesian box radius ε | Optional contacts | Certified pairwise orders | Guaranteed top-10%-rank membership |
| --- | ---: | ---: | ---: |
| 0.0005 Å | 0 | 3,415,191/3,415,191 (100%) | 262/262 |
| 0.005 Å | 5 | 3,361,308/3,415,191 (98.42%) | 258/262 |
| 0.05 Å | 108 | 2,270,083/3,415,191 (66.47%) | 98/262 |

The guaranteed top-decile is the subset of nominal top-decile-ranked model
pairs whose rank remains within the top decile for **every** admitted contact
network and every coordinate in the box. It is NOT the fraction of experimentally
observed pairs predicted correctly.

### CspB and crambin transfer checks (existing revealed-reference data)

The original publicly reproducible exact-rational coordinate reference file,
`FoldBounds_Frozen_Exact_Rational_Inputs_CORRECTED_V2.json`, from the archived
`FoldBounds_Bulk_Rational_Audit_COMPLETE_20261104.zip` (container-side SHA256
`400dd9a242c1f34c71781f50079ad87f8a394eb03b3fc309aa71736aa2595923`)
was reused without accessing new NMR targets.

| Protein | GNM contact count | Evaluated pairs | Nominal adjacent order certs | ε=.005 Å certified pairwise orders | Guaranteed top 10% |
| --- | ---: | ---: | ---: | ---: | ---: |
| CspB | 543 | 1170 | 1160/1169 | 673,132/683,865 (98.43%) | 116/117 |
| Crambin | 334 | 547 | 546/546 | 143,811/149,331 (96.30%) | 50/55 |

CspB had some true or interval-overlapping ties; they are NOT silently
converted into strict ranks.

## Biology and model development

The **prespecified RNase A weighted ANM model failed**: Spearman ρ=0.421357,
strongest baseline GNM resistance ρ=0.511009; margin -0.089652, against a
required positive margin +0.15. This failed test is not replaced or rescored.

A source-only generalized ANM model was implemented with

    K_e = (7Å/d_e)^2 * [(1−α)u_e u_e' + α I_3], α = 0.5,
    H = Σ_e [[K_e,-K_e],[-K_e,K_e]] + 0.2 I.

It reproduces the original weighted ANM exactly in the α=0 limit. The
isotropic α=1 limit is NOT the original unweighted, zero-ridge GNM. In
already-revealed RNase A, α=0.5 gives exploratory ρ≈0.517420; a secondary
75% GNM rank +25% weighted ANM rank fusion gives ρ≈0.522767. These results
were observed **after** RNase A disclosure and must never be described as a
blind successful study. These small gains do not meet the original +0.15 gate.

G-ANM model families predate this work: Zheng (2008), *A Unification of the
Elastic Network Model and the Gaussian Network Model*, Biophysical Journal,
DOI:10.1529/biophysj.107.125831. Original GNM: Haliloglu/Bahar/Erman (1997),
PRL 79:3090, DOI:10.1103/PhysRevLett.79.3090.

## Truly untouched prospective test: thioredoxin

- Reference **2TRX**, E. coli oxidized thioredoxin crystallographic chain A
  (two copies in the crystal asymmetric unit; chain A only, state 1).
- Target **1XOA**, same E. coli oxidized thioredoxin NMR **20 conformers**.
- These identifiers, organism, redox states and experimental methods were
  checked against RCSB structural metadata; **target coordinates were NOT
  downloaded or accessed** for model fitting here.
- Stage 1 script `scripts/fb_thioredoxin_stage1.py` fetches only 2TRX, produces
  all 8 predictor arrays and protocol, records actual source code and hashes,
  and emits a Stage1 ZIP for independent audit. Mock PyMOL tests demonstrate
  2TRX-only fetch behavior. Actual PyMOL run awaits the user's installation.
- Stage 2 script `scripts/fb_thioredoxin_stage2.py` **requires an externally
  supplied 64-character SHA256** of the actual Stage1 ZIP and verifies it
  BEFORE any NMR fetch, then compares on the exact matching residue intersection.
  Both the Stage1/Stage2 scripts are implemented and tested on mocked data,
  but the real 1XOA result remains unknown.
- **Primary** generalized ANM α=0.5 must achieve Spearman ρ≥0.50 and a margin
  ≥+0.15 versus the strongest of *all six* original baselines; a secondary
  75:25 GNM/weighted-ANM rank-fusion comparison is exploratory.
- At least 85% of reference eligible residue pairs must survive identity
  matching, and at least 10 NMR conformers must be present, or abstain.
- Distinguish NMR conformer refinement spread, experimental redox states,
  crystal lattice conditions, and actual dynamic populations.

## Formal trust boundary and outstanding work

Lean module `MotionGNMExact.lean` checks an **exact rational inverse witness**
of the scalar grounded Laplacian, along with the complete reference graph.
This establishes a completely checked exact score **given a supplied inverse**.
The `gnm_residual_enclosure` theorem rigorously derives rational score bounds
from a symmetric matrix and a stated global coercivity premise.

`MotionGNMTopologyRobust.lean` kernel-checks the monotonicity of scalar
spring energy and the global Rayleigh response inequalities over arbitrary
variable-contact assignments. However:

- The **full n=124 grounded Poincaré/BFS lemma is not instantiated into Lean**;
  it is proved mathematically above and numerically applied with independent
  exact Python checks. Full coordinate bytes→Lean graph→integer witness→Lean
  score rank replay remains a future engineering gate.
- The numerical certificate source is exact and independently replayable in
  Python, but its Python implementation/runtime and binary64 decoder are not
  in Lean's trusted computing base. No Lean theorem about these exact NPZ
  archive bytes is claimed.
- The explicit coordinate boxes are mathematical thought-experiments, not
  empirical confidence intervals for protein structural error.
- Thousands of shared-residue pair comparisons are dependent observations;
  no independent-sample p-values or guaranteed biological correctness follow.

## Reproduction (local Python 3.11+ with NumPy, SymPy, pytest)

From the repository branch, with the original Stage1 ZIP at `STAGE1.zip`:

```
python scripts/gnm_exact_certificates.py STAGE1.zip gnm_cert.json
python scripts/gnm_exact_certificates.py STAGE1.zip gnm_cert.json --verify
python scripts/gnm_independent_replay.py STAGE1.zip gnm_cert.json
python scripts/gnm_topology_boxes.py STAGE1.zip gnm_box.json
python scripts/gnm_independent_box_replay.py STAGE1.zip gnm_box.json
python -m pytest -q scripts/test_gnm_exact_certificates.py scripts/test_gnm_topology_boxes.py scripts/test_hybrid_and_blind_stage.py scripts/test_thio_stage2_mock.py
```

All Stage1 model results are reference-only; no source data or predictive model
are silently updated after seeing a target.

## Retrospective target-ranked concordance interval (new independent bridge)

With the separately revealed and SHA-256 verified RNase A 7RSA->2AAS
32-conformer NMR dataset (Stage2 digest
`4f1ea752ada45453a44ea9128b6e7a43f4070e6824062b796a27144c340dfc6c`),
we applied the **reference-only** exact-integer GNM structural uncertainty
intervals to all target-nontied pair comparisons. None of the predictor
parameters or original blind acceptance criteria were altered.

For target pair distance variance values `t_i`, define pair comparison
contribution `sign((q_i-q_j)*(t_i-t_j))` (0 for predicted tie), and average
over distinct `t_i != t_j`. This is Kendall-style *signed concordance*,
NOT Kendall tau-b, Spearman rho, a confidence interval, a p-value, or an
independent trial. For every world whose GNM contact network is induced by
coordinates in the stipulated box, disjoint certified response intervals
fix the contribution; ambiguous intervals contribute conservatively
`[-1,+1]`. We DO NOT assume independent rank comparisons or attainable
simultaneity of every uncertain sign.

| Reference Cartesian coordinate box (each component) | Certified signed concordance bounds | Certified comparison count |
|---|---:|---:|
| ±0.0005 Å | `[404305/1138397, 404305/1138397]` ~`0.355153` | `3415191/3415191` |
| ±0.005 Å | `[1161499/3415191, 1269265/3415191]` ~`[0.340098,0.371653]` | `3361308/3415191` |
| ±0.05 Å | `[-15117/1138397, 2244865/3415191]` ~`[-0.013279,0.657318]` | `2270083/3415191` |

At ±0.005 Å the **defined model** retains positive directional agreement
with the NMR structural ensemble labels even for every admissible contact
switch. This is a conditional proof about an unweighted GNM graph model and
fixed deposited target labels, NOT an empirical guarantee of dynamic
prediction quality or a new blind success.

Reproduce with `scripts/gnm_target_concordance.py` and the source Stage1,
source Stage2, plus `gnm_structural_robustness.json`. The command
**first invokes** an independently written integer-coordinate source
replay (`gnm_independent_box_replay.py`) of all endpoint response envelopes;
then checks source freeze, target 32-conformer arrays, label regeneration,
and computes exact integer interval-comparison counts using a sorted-target
algorithm. `scripts/test_gnm_target_concordance.py` (exhaustive brute-force
comparison against independent O(n²) logic on thousands of small adversarial
box/tie cases) passed. Exact rational endpoint bounds are checked in Python;
the formal proof `MotionRobustConcordance.sum_alignedSign_interval` is in
Lean, but the 3,415,191 dataset-specific comparisons have **not** been
kernel-replayed in Lean.

The large ±0.05 Å enclosure includes negative concordance; the system
correctly abstains from guaranteeing useful agreement at that radius.

## Cross-protein falsification of fixed topology/anisotropy mixtures

The physically motivated 50% isotropic/50% directional hybrid showed a small
post-reveal gain on RNase A (rho 0.5174 vs 0.5110 GNM). However, the exact
same model on the **previously revealed** CspB and crambin data is inferior to
the original weighted ANM, using identical reference-only pair masks:

| Protein | GNM rho | Weighted ANM rho | 50:50 G-ANM rho | GNM 75%/ANM 25% rank-fusion rho |
|---|---:|---:|---:|---:|
| RNase A | 0.5110 | 0.4214 | 0.5174 | 0.5228 |
| CspB | 0.1856 | 0.4104 | 0.3427 | 0.2610 |
| Crambin | 0.1175 | 0.4154 | 0.2902 | 0.1979 |

The observed variation makes a universal fixed hybrid unsupported. It does NOT
invalidate exact GNM model certificates. The thioredoxin Stage1 model is **not
silently retuned** to fit this exploratory observation; its future target is a
falsification test, with the original stringent gate preserved.

A second exploratory **selective model-consensus** diagnostic examines only
pair order comparisons where frozen GNM and weighted ANM agree. This cannot
be called uncertainty-certified dual-model inference because the complete
ANM coordinate-box certificate is unavailable for these 2,614 RNase pairs.
It is also a post-hoc target-revealed diagnostic (all pair comparisons are
statistically dependent):

| Protein | GNM/ANM pair-order consensus coverage | Target-correct among agreed pairs |
|---|---:|---:|
| RNase A | 73.70% | 71.88% |
| CspB | 69.69% | 64.59% |
| Crambin | 76.16% | 61.88% |

The consensus signed agreement is 0.438 (RNase), 0.292 (CspB), and 0.238
(crambin). In the crambin case the selective consensus is **worse than the
weighted ANM's overall signed agreement** (~0.283), so consensus is not
universally superior either. This is an exploratory coverage/accuracy
trade-off, NOT a prospective win. Scripts are
`scripts/audit_exposed_prior_proteins.py` and
`scripts/model_consensus_audit.py`; the input user archives are pinned by
SHA-256 and tests may be rerun locally with NumPy/SciPy. See the companion
machine-readable JSON audits.