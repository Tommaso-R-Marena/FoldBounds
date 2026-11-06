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
