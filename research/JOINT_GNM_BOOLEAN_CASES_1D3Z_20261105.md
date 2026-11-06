# Reference-only, exhaustive contact-case refinements of exact GNM rank certificates

**Status:** computational verification advance in a known Gaussian network model (GNM); **not** new protein-motion prediction, a new biophysical law, or a new rank-one matrix theorem. This contribution is a new isolated PR stacked solely on the exact joint-covariance PR #76, itself directly based on `main`.

## Scientific context and input

The biological source is the public deposited ubiquitin solution NMR ensemble **1D3Z**, using the publicly archived `encore-similarity/encore` filtered first-model Cα coordinates, residues 2–70 (69 residues). The source/model 1 fixture is independently checked against the pinned Git blob in PR #76. Coordinates are exact to 0.001 Å as deposited. **No NMR model-ensemble variability, target motions, or biological labels are consumed by the split selector or by the certificate generator.**

Define a binary contact network: all consecutive chain bonds are present, other edges require |i-j| >=3 and <8 Å. The structural uncertainty is an illustrative **independent ±0.05 Å Cartesian coordinate box**, not experimentally calibrated. Using the conservative distance slack 7ε/2, the certain graph has 152 nonlocal contacts and the possible graph 167, so 15 nonlocal contacts are optional. All contact subsets are physically treated as a conservative superset: some assignments might not correspond to any single coordinate realization, but all realizable graphs are included.

The exact Cα GNM covariance is the graph-Laplacian pseudoinverse. A rank ordering is a sign of C_ii-C_jj, not a measurement of true protein mobility. The existing PR #76 computes a **rank-two semidefinite interval** for each such comparison, including the off-diagonal covariance of endpoint differences, using exact integer inverse certificates.

## A universal Boolean partition of contact uncertainty

Select s distinct optional contacts e1,...,es. For each of the **2^s** assignments to their on/off statuses, form a separate interval: add selected-on contacts to both endpoints; remove selected-off contacts from both. All other optional contacts remain uncertain. These intervals cover **every** admissible contact graph, and a strict ordering holds globally precisely when **every** case independently certifies the same strict sign.

A case certificate is accepted only after an **independent exact integer verifier** reconstructs each endpoint's `B=nL+11^T`, checks its integer matrix-inverse witness `B*N=dI` across **all n² entries**, verifies symmetric and positive denominators, then checks the rational rank-two signs. Final signed comparisons are *the intersection of equal-sign conclusions across every Boolean case*, not a union of case findings. This gives a deterministic, fail-closed, model-wide sound implication; incomplete case coverage or contradicted signs produce no global ranking.

**Known math:** graph Laplacian Rayleigh/Loewner monotonicity, rank-two SDP interval optimization, Boolean partition coverage, and Sherman–Morrison rank-one updates are classical. No mathematical-priority claim is made.

## Exact real-data output (no biological label tuning)

The first two selected optional contacts were found by a source-only greedy search maximizing the number of **provable GNM rankings**, not observed NMR agreement. The fixed, independently replayed selections are PDB residues **43–51** and **39–42** (zero-based indices 41–49 and 37–40). They partition all 2^15 optional-edge assignments into **four** mutually exclusive and exhaustive classes.

| Exact method | Certified strict residue-pair orders / 2346 | Additional over prior |
|---|---:|---:|
| Separate diagonal endpoint intervals | 1884 | baseline |
| Joint rank-two covariance endpoint interval | 1919 | +35 |
| Joint rank-two with 2 selected contacts, 4 cases | **1941** | **+57** over separate / **+22** over joint |
| **Completed depth-three greedy, 3 contacts, 8 cases** | **1947** | **+63** over separate / **+28** over joint |

For the two-contact design, every case verifies two 69×69 inverse witnesses: eight final endpoint identities, plus a separately verified unsplit baseline. The final signed list contains **1,941** rankings. At depth three, the eight exhaustive cases verify **16** endpoint inverse identities and certify **1,947/2,346** rankings (approximately **83.0%**). These are exact acceptance counts, not sampled estimates.

A simpler source-only trace-susceptibility selector, which ranks all uncertain contacts by the exact single-edge reduction in the trace of the certain-network covariance, takes PDB contacts (39,42), (46,67), (50,59). These three edges yield **eight** cases, **16** checked endpoint matrices, and **1,939** globally certified orders. The source-only greedy 2-edge method therefore gives more coverage at lower final-certificate cost. This is an **exploratory** optimization of model-certificate coverage on this particular reference structure, not a learned physical model nor a general runtime scaling guarantee.

## Faster untrusted proposer with unchanged trusted checks

To avoid reinverting a large exact matrix for every optional-edge set, `gnm_integer_updates.py` proposes a classical exact rank-one update. With an accepted parent inverse `B N=dI` and new contact incidence b=e_i-e_j (the scaled matrix changes by `n bbᵀ`), let v=Nb, α=bᵀNb. Then propose

    N' = (d+nα) N - n v vᵀ,
    d' = d(d+nα).

Remove the common integer GCD from `(N',d')`. **Neither producer nor cache is trusted**: accept only if the independent integer verification of the resulting edge set gives `B' N' = d' I`. The higher-level rank certificate independently rechecks each endpoint anyway, so the cached proposal cannot bypass the checker. A separate exact rational Fraction-oracle and direct SymPy inversion comparison test the proposal on dozens of seeded small connected graphs and poisoned-parent cases.

The **previously unfinished depth-three greedy search has now completed** using the existing exactly verified integer-inverse cache. Starting from the selected contacts PDB (43,51) and (39,42), the third greedy contact was **PDB (7,13)** (0-based indices (5,11)), improving certified rankings from 1,941 to **1,947**. The independently verified eight-case result reproduces all 1,947 accepted signs over all 15 optional contacts; the source-only producer never uses NMR variability or target labels. The depth-three choice is **greedy conditional on the already selected first two contacts, not a proven globally optimal 3-contact subset**. A few next-best candidate contacts produced 1,946 rankings. The trace-heuristic eight-case method remains a separate baseline at 1,939. The optional greedy search is much slower than a fixed-certificate replay and is deliberately not run inside hosted CI; instead, CI directly verifies the frozen source-only three-contact certificate.

## Reproduction

From this stacked branch in a Python 3.11+ environment:

```bash
python -m pip install -e '.[test]' sympy
python -m pytest -q tests/test_gnm_integer_updates.py tests/test_gnm_joint_cases.py
python experiments/real_1d3z_joint_cases.py --epsilon 0.05 --selected-count 2 --fixed-two
python experiments/real_1d3z_joint_cases.py --epsilon 0.05 --selected-count 3
# Completed depth-three greedy selection, exactly checked without NMR labels:
python experiments/real_1d3z_joint_cases.py --epsilon 0.05 --selected-count 3 --fixed-three
# Optional slower replay of the source-only greedy candidate search:
python experiments/greedy_1d3z_cases.py --depth 3
```

The 7 new tests pass locally. Combined with PR #76's 10 tests, **17/17** relevant tests pass. New dedicated CI checks the small exhaustive Boolean graph oracle, real pinned-coordinate four-case result, update-identity tampering, source-only heuristic result and both reproducible demos on Python 3.11, 3.12 and 3.13. The optional greedy search is not required in hosted CI due to cost.

## What is NOT established

- This is a **certificate-strength improvement**, not a new GNM scoring rule; biological prediction quality is unchanged. No new held-out or blind accuracy assessment has been conducted.
- No result on the sealed thioredoxin 1XOA target or the previously failed RNase A 7RSA→2AAS weighted-ANM blind experiment. They are untouched.
- NMR ensemble spread is not a kinetic trajectory, independent solution-state fluctuation measurement, or proof of equilibrium dynamics.
- ±0.05 Å boxes and cutoff choice are illustrative; no empirical coordinate-error calibration is claimed.
- Independent exact-integer Python checking is not an end-to-end Lean kernel refinement. Passing an unrelated Lean workflow does not formally certify these new source-specific arrays.
- The case-split method has exponential complexity in the *number of selected contacts*; it is not a general resolution of very-large-network computational complexity.

The next scientifically consequential milestone is a preregistered **held-out experiment** using actual motion-related observations from more than one independent protein, including negative baselines, plus extension of certified geometry/energy bounds to more realistic physically calibrated models.

## Completion and reproducibility correction

The initial hosted workflow failed despite the then-15 passing unit tests, because executing `python experiments/real_1d3z_joint_cases.py` did not place the parent directory in the module search path. The direct/package import was repaired in both experiment entry points on PR #79; the CI on the fixed prior head then passed all Python and Lean workflows. An explicit subprocess regression test now compares direct CLI and `python -m` JSON outputs and exercises `greedy_1d3z_cases.py --help` in both modes. The new three-contact result is independently reproduced by the exact eight-case checker in both a unit test and a direct CLI invocation in hosted CI. The previous two-contact and trace-heuristic figures remain unchanged.
