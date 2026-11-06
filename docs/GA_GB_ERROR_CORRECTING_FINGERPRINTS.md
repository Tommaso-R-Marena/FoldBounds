# Error-correcting geometric fingerprints from real GA95/GB95 NMR coordinates

**Experimental input:** RCSB 2KDL (GA95) and 2KDM (GB95), authentic 20-model 56-residue solution NMR structural ensembles; data downloaded afresh in GitHub CI with SHA-256 byte verification. The parent research branch \`research/real-ga-gb-experimental-benchmark-20261020\` provides strict chain A and sequence checks. Their 56-residue sequences differ only at positions 20, 30 and 45. This is a **retrospective mathematical classification of deposited models**; it is neither a mutation-to-fold predictor nor a measurement of the solution-fold equilibrium populations.

## Motivation and non-overlap

The existing FoldBounds GA95/GB95 work found a very strong **single** distance reporter Q11–V21 (~16.37 Å minimum band gap). However, a single assay becomes inconclusive if its label is impaired, its fluorescence is ambiguous, or the measurement is corrupted. The new problem is a **fault-tolerant set of reporter distances** in which:

- Every reporter uses only **unchanged residue positions**, with sequence separation ≥4.
- All observed C-alpha pair distances across ALL 20 GA95 and 20 GB95 deposited models must lie within a conservative proxy readout window of **12–40 Å**.
- The two complete NMR distance intervals must be disjoint for EACH reporter.
- Selected reporter pairs must use **mutually distinct residues** (no reporter site shared with any other).
- A single common pair-distance error tolerance must work across ALL models and all selected channels.
- For K=5 or K=7 channels, majority classification must survive at most 2 or 3 **arbitrarily wrong channel decisions**, respectively.

This is a **structural readout proposal**, not a physical demonstration that labeling each of 10 or 14 sites preserves the native fold. Label chemistry is an extremely substantial limitation and must be tested separately.

## A verifiable optimum with integer arithmetic, not an opaque solver

The coordinate PDB columns have three decimal places in Å. Convert each 3D C-alpha coordinate into exact **integer milli-Å** and compute the squared distance between residues i,j entirely with integers. The floor/ceiling of the exact integer square root conservatively enclose the Euclidean distance with a maximum 1 milli-Å outward error at each endpoint. Consequently every reported inter-fold band gap g (integer milli-Å) is a **certified lower bound** on the true worst-case separation of all submitted model coordinates, not a floating-point approximation that can round away a threshold failure.

To make the matching optimum independently checkable, restrict candidate readouts to N-terminal indices 1–14 paired with indices 15–56 (still excluding 20,30,45). This **explicitly defined bipartite candidate graph** permits a compact combinatorial proof. For every candidate residue-pair edge e, its weight g_e is the conservative guaranteed inter-fold distance-band gap.

Find the largest integer margin g* for which there are K **vertex-disjoint** edges with g_e ≥ g*. The checker constructs K such edges using exact augmenting-path matching.

To prove no stronger K-reporter solution exists in the **same declared domain**, separately construct a vertex cover of the graph containing edges with g_e > g*. Its verified cardinality is **strictly less than K**; since every matching must use distinct vertices and every edge touches a cover vertex, a matching with K edges is impossible at a stricter margin. The receiver independently reconstructs ALL candidate edges from raw authenticated coordinates, verifies the chosen K disjoint edges and their gaps, and checks every stronger edge is covered by the displayed vertex-cover witnesses. This is a simple **primal-feasible matching + dual-like upper obstruction**, relying on classical matching theory (König's theorem for why an optimal small cover can be found).

It is not a claim of the global optimum over all 56 choose 2 residue pairs, all labeling chemistries, or all assays. The specific N-terminal division and proxy distance range are intentional auditable design constraints.

## Coding theorem and physical interpretation

For a chosen reporter with GA distance band [a_lo,a_hi] and GB band [b_lo,b_hi], suppose without loss of generality a_hi < b_lo. The fixed midpoint of the two bounds separates the two **finite model sets** by at least g/2 on each side. Any uncorrupted pair-distance measurement deviating by strictly less than g*/2 from its target deposited-model distance will be classified correctly. With K=2f+1 reporters, at most f arbitrary erroneous channel votes leave f+1 correct votes, guaranteeing a correct binary majority.

- For K=5, this corrects **up to two** entirely wrong channels.
- For K=7, this corrects **up to three** entirely wrong channels.
- If each endpoint C-alpha position is changed by at most r in Euclidean distance, then by the reverse triangle inequality a pair distance changes by at most **2r**. Thus a sufficient per-endpoint coordinate error budget is r < g*/4. This is a **mathematical perturbation bound**, not a calibrated physical dye/linker precision.

No distributional or independent-error assumption enters the *conditional* majority theorem. But all good reporters must be valid at the supplied unknown structural state, and the deposited NMR models are not an exhaustive or representative sampling guarantee for arbitrary solution ensembles. Correlated label perturbations affecting more than f reporter decisions defeat the guarantee.

## Independent tests and real-data acceptance gate

\`tests/test_ga_gb_error_codes.py\` uses a completely synthetic 20+20 coordinate fixture and an independent exhaustive oracle for small bipartite matchings. Tests check conservative exact integer square-root intervals, the cardinality matching and cover certificate, every pattern of up to f arbitrarily corrupted channel outcomes on 40 synthetic models, and rejection of mismatched sequences, nonexact coordinates, forged margins and insufficient channels. It does not substitute these fixtures for real PDB data.

The separate GitHub Actions \`independent-real-nmr\` workflow:
1. Downloads public **official** 2KDL and 2KDM from RCSB.
2. Fails closed on either published-coordinate SHA-256 mismatch.
3. Strictly reparses the 20×56-residue chain-A C-alpha coordinates and authenticates the expected 3-residue sequence difference.
4. Independently solves the 5- and 7-channel N-terminal/C-terminal reporter selection problems on ALL 40 deposited conformers.
5. Replays each model's exact CA distance observations through the majority classifier and enumerates **every** corruption subset of sizes up to 2 or 3.
6. Uploads the results as a machine-readable JSON run artifact and prints it in the GitHub Actions log.

The final mathematical bounds and specific residues must be read from a **successful real-data run**; no data-dependent performance is claimed before that run completes.

## Limitations and scientific acceptance criteria

- **Prior art:** distance reporter optimization, bipartite matching, vertex covers and error-correcting majority codes are established techniques. This is a new FoldBounds application and exact assurance workflow, not an asserted original mathematical theorem or experimental breakthrough.
- **Not blind:** all NMR models are used to optimize the reporter set. Testing it on the same models establishes internal consistency, not prediction accuracy on previously unseen structures. A serious next gate is held-out or prospective mutants and experimentally measured labeled distances.
- **Not an actual FRET measurement:** C-alpha distances are not fluorescent dye–dye distances. Ten or fourteen actual labels would be unrealistic for many single-chain assays; the mathematical channels could instead be measured across separately prepared constructs, which introduces comparability and perturbation controls.
- **Not a fold-population estimate:** NMR models are selected conformers, not equilibrium samples. The error-correcting code makes NO claim about GA95↔GB95 equilibrium ratios or folding rates.
- **Not thermodynamically novel:** it certifies **structural separation** within the submitted models. The model does not infer temperature-dependent free energy or a new physical folding mechanism.

Reproduce the real-data audit:

    python -m pip install -e '.[test]'
    python experiments/ga_gb_error_code_demo.py --ga /path/2KDL.pdb --gb /path/2KDM.pdb

Inspect \`.github/workflows/ga-gb-error-codes.yml\` for exact download, provenance pins, independent tests, and verification report. The code remains draft pending CI, expert scrutiny, and real experimental validation.


## Verified real-data results, including the unsuccessful holdout tests

The official RCSB coordinate files downloaded by GitHub Actions passed
the exact pinned SHA-256 checks, full 20×56-residue model completeness,
and expected sequence/mutation checks. At the latest successful
independent-real-NMR run, the exact conservative bipartite candidate
graph had **170** eligible positive-gap edges. The following results
were computed from the **actual 2KDL/2KDM C-alpha coordinates**, not
from the synthetic unit-test fixtures:

| Number of disjoint reporter pairs K | Corrupted readout decisions correctable f | Optimal minimum NMR distance-band gap (Å) | Strict pair-distance observation error allowance (Å) | Required distinct residue sites |
|---|---:|---:|---:|---:|
| 1 | 0 | 16.368 | <8.184 | 2 |
| **3** | **1** | **12.948** | **<6.474** | **6** |
| **5** | **2** | **11.580** | **<5.790** | **10** |
| **7** | **3** | **10.302** | **<5.151** | **14** |
| 9 | 4 | 9.081 | <4.5405 | 18 |

Every full-data code correctly classified **all 40 deposited NMR
conformers**, and passed every re-evaluated adversarial bit-flip subset
up to its listed budget. The full-data codes are NOT prospective tests:
these same 40 models were used to choose the reporters.

The six-site, three-reporter result is the least demanding option that
tolerates one wholly incorrect channel; the 10-site/5-reporter option
is a useful mathematical compromise. These are *distinct amino-acid
sites* that would require labeling, rather than simply five or seven
readout values from one verified assay. The chemical feasibility and
absence of fold perturbations remain untested.

### Fixed within-deposition 10/10 holdout: important negative controls

To evaluate transfer to unseen structures *within the same deposited
NMR entries*, a separate protocol fitted the reporter choices on
models 1,3,5,...,19 and tested on 2,4,6,...,20, then reversed the
split. There was no outcome-dependent model resampling. Note that
individual deposited models are still **not independent biological
replicates** and can have refinement-related correlations.

Every heldout half had **20/20 correct clean binary fold
classifications**, regardless of K. When the stated f channels were
deliberately forced to vote for the wrong fold, however, the reversed
split sometimes failed:

| K | Allowed wrong channels f | First-half-fit adversarial successes | Reversed-fit adversarial successes |
|---|---:|---:|---:|
| 1 | 0 | 20/20 | 20/20 |
| 3 | 1 | 80/80 | **78/80** |
| 5 | 2 | 320/320 | **314/320** |
| 7 | 3 | 1,280/1,280 | **1,240/1,280** |
| 9 | 4 | 5,120/5,120 | **4,980/5,120** |

The stricter reversed-fit evaluation had one or more channels with
**negative signed classification margins** in the omitted structural
models, as low as -6.305 Å. A code trained on that half of the data
therefore **does not** carry its training-derived deterministic
worst-case guarantee over the heldout models, even though ordinary
uncorrupted majority prediction happened to be perfect. This
distinction is central to the scientific interpretation.

An earlier exploratory run perturbed corrupted channels using the
same numerical extremum irrespective of each reporter's fold
orientation. That is NOT the strongest adverse binary error.
The later data above use the corrected test that deliberately forces
each specified channel to vote for the OPPOSITE fold, so the smaller
reversed-split counts here supersede those early exploratory numbers.

These results can support a practical future design protocol only
after independent labeled-construct feasibility, state-specific
dye-distance calibration, genuinely new conformers/mutants, and
measurement-condition controls. This is a proof-carrying
**structural-geometry benchmark**, not evidence of a novel folding
mechanism or measured mutation-dependent populations.
