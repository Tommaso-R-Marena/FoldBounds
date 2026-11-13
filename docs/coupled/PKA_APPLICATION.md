# Multi-site protonation model with separator certificates (application of CEC)

Status: results of a pre-registered run (2026-11-13). No novelty is claimed: the model is a
standard Coulomb/Ising titration model, the certificates apply the separator and T1 enclosures of
`CEC_THEORY.md` (sections 1, 2). The certificate concerns the MODEL, never experiment.

## 1. Model and provenance
Frozen in `experiments/coupled_pka/PROTOCOL.md` (commit fbaf401, before any experimental comparison)
plus `PROTOCOL_ADDENDUM_1.md` (commit d6a8057, also before any comparison). Summary: binary
protonation states; model pKa from PROPKA 3.5.1 `propka.cfg` (sdist SHA-256 d8083b58...2257);
site-site Coulomb with eps(r)=4r, hard cutoff 12 A, T=298.15 K, no desolvation or screening; couplings
quantised to 1/1000 pK; every 10^(e/1000) replaced by a 15-significant-digit rational, so the model is an exact rational
antiferromagnetic MRF (coefficient of sigma_i sigma_j is +w_ij > 0 for all pairs). Structures: RCSB PDB files;
experiment: PKAD-2 page (https://compbio.clemson.edu/pkad-2/, retrieved 2026-11-13, SHA-256 a71dfa05...8123).
All URLs and hashes: `data/coupled_pka/MANIFEST.json`. Model specs (integers and decimal strings): `data/coupled_pka/specs/`.

## 2. Method
For site i and pH, Q_i = P(sigma_i=1)/P(sigma_i=0). Interior S = i plus up to 8 nearest sites within R
(R in 0,4,6,8; 10 only in the small exact tests); the boundary is the set of outside sites coupled to S. The
boundary cube tree (static order, strongest coupling first) partitions {0,1}^boundary; on each cube the regional A-masses and
regional ratios are bracketed (separable in the boundary), and the T1 test (Thm 1.1/1.2) yields a bound on Q. The tree
is expanded best-first up to a leaf budget; fully expanded it is exact enumeration of the boundary (the hard
separator enclosure of Thm 2.1); R=0 is the closed-form separable case.
`prop.py` (numpy floats) proposes; `checker.py` (Fractions only, own weight derivation) verifies the order, the tree and
the exact T1 test on every leaf. Claims are backed off 1e-9 relative.
**Soundness argument.** For every boundary state tau in a cube, d_xi(tau) lies in [L,U] and r_xi(tau) in
[rmin,rmax] (products of per-boundary-site minima/maxima, exact); Q(tau)=sum r d / sum d; Thm 1.1 shows g_-(lambda)>=0
implies Q>=lambda on the admissible box, which contains the true values. A partition of the boundary states
into cubes then gives the bound for all tau, and, since no factor joins interior and exterior (hard cutoff),
the exterior only enters through tau (Lemma 2.0, Thm 2.1). The checker re-establishes each step exactly.
Per-site monotonicity in pH is NOT assumed (see 4.1); the pKa bracket (integer grid 1..13, 3 bisections, certified
decisions at Q=1) asserts only that a crossing exists inside it.

## 3. Exact validation (<= 20 sites, brute force)
BPTI (4PTI, 20 sites) and four 20-site subsets (nearest sites to site 0 of 2LZT, 3RN3, 1STN, 2RN2).
Exact marginals by integer depth-first enumeration of 2^20 states. 5 systems x 5 pH (3,5,7,9,11) x 20 sites x 5 radii x 4
budgets (1, 8, 64, 16384) = 10,000 certificate pairs: **0 violations (exact Q inside every certified interval), 0
checker rejections** (`results/validate_exact.json`). Mean theta width (budget 1 / exhaustive):
R=0 0.336/0.336, R=4 0.322/0.273, R=6 0.254/0.176, R=8 0.199/0.108, R=10 0.152/0.053.
Width decreases with R and budget. A float64 variable-elimination reference
for 4HHB (174 sites, min-fill width 17; `results/mid_exact_4HHB.json`): 3480 certificates (5 pH x 4 radii), 0
violations, 338 bracket endpoint checks, 0 violations; mean theta width at pH 3-11 by R: 0.466, 0.361, 0.224, 0.156; 2913 certificates
checked exactly, 0 rejected.

## 4. Results
### 4.1 Non-monotone titration (negative for a simplification)
The total protonation is monotone in pH, but single sites are not: float scans of the exact systems (pH 0-14, step 0.25)
find increasing theta_i in 16 of 100 sites, using an increase threshold of 1e-12. The bracket procedure therefore
does not assume monotonicity and flags several crossings (`multi`).
*Correction after adversarial review (P2-m1).* The BPTI increases are physically negligible. The independently
reproduced maximum increases are GLU7 9.7e-11, TYR10 2.9e-5, ARG20 3.9e-6 and TYR35 2.4e-6. The largest increase,
0.013, occurs in the 3RN3/2RN2 subsets, not in BPTI. Non-monotonicity is real in the model, but "16 of 100" counts
negligible increases. Non-monotone single-site titration in coupled systems is a known phenomenon in the titration
literature (e.g. Onufriev, Case and Ullmann 2001; citation not re-verified here).

### 4.2 Large system: GroEL 1OEL (7 chains, 1064 sites, 4035 couplings)
Greedy min-fill treewidth bound 31: tables of 2^32 doubles (32 GB) exceed the 4 GB budget; exact inference not run.
pH 7, 1024-leaf budget (all certificates for R=8 verified exactly: 19,482 accepted, 0 rejected; CPU
2238 s in 16 chunks, most of it in the checker):

| R | mean width | median | 90th pct | frac width<=0.01 | <=0.1 | >0.5 |
|---|---|---|---|---|---|---|
| 0 | 0.434 | 0.249 | 0.999 | 0.258 | 0.427 | 0.431 |
| 4 | 0.305 | 0.031 | 0.992 | 0.437 | 0.572 | 0.293 |
| 6 | 0.165 | 0.0010 | 0.808 | 0.627 | 0.720 | 0.161 |
| 8 | 0.099 | 0.0001 | 0.323 | 0.680 | 0.807 | 0.082 |

Non-converging sites (theta width > 0.1 at R=8): 205 of 1064 (19%): GLU 120/322, ASP 37/231, LYS 11/273,
TYR 14/49, CYS 14/21, HIS 7/7, N-term 2/7, ARG 0/147, C-term 0/7. They form 57 connected clusters in the coupling graph (largest 7 sites),
so non-convergence is concentrated in small frustrated clusters of acids, not in one giant component. Failure rate does not increase
monotonically with a site's total coupling strength (by quartile: 0.21, 0.33, 0.11, 0.11), so the simple strength proxy does not explain it.
We did not test whether this is non-uniqueness (weak-spatial-mixing failure, Thm 2.8) or only budget/interior-cap limits; the
width keeps shrinking with R and budget in all tests, which is consistent with either.
pKa brackets at R=8: single bracket 364, none-above 392 (theta>1/2 certified across pH 1-13, mostly Lys/Arg), none-below
308; widths of single brackets: mean 3.9 pH units, only 9.6% <= 1 (bisection stops at the first undecided midpoint).
Monte Carlo (Gibbs, 3 seeds x 64 chains, 1000 burn-in + 4000 sweeps, pH 3,5,7,9,11): the MC mean lies inside the certified
interval for 53-59% of (site, pH) pairs.
**Pre-registered criterion (MC mean inside the interval widened by 3 standard errors): 55.8 / 56.5 / 58.4 / 65.3 / 54.9 %
inside at pH 3 / 5 / 7 / 9 / 11, with 369-480 of 1064 sites outside at each pH.** This figure was reproduced by the
adversarial review (P2-M1).
An earlier version of this document reported "100% (0 outside)". That figure used an additional absolute slack of +1e-3
which the protocol did not register. It is a **post-hoc deviation** (§6) and must not be quoted as the registered result.
Likely benign cause, not proven: many sites are frozen (fully protonated or deprotonated) in every MC sample, so their MC
standard error is exactly 0, while their certified interval is a tiny interval that excludes exactly 0 or 1. Under the
post-hoc +1e-3 slack, 100% of MC means are inside, and the largest strict overshoot is 0.003. MC agreement is a consistency
check, not a proof, and it does not affect certificate soundness.

### 4.3 Pre-registered comparison with experiment (PKAD-2, 9 proteins)
134 sites matched (16 non-numeric records and 3 unmatched sites excluded). Model midpoint available only for sites with a
single certified bracket (75); 59 sites (24 ASP, 21 GLU, 4 HIS, 3 C-term below; 9 above) have no pKa in the 1..13 grid at R=8.
Predictions of the model without experimentally tuned parameters:

| set | n | RMSE model | RMSE null | Delta | 95% bootstrap CI of Delta | Pearson r (mid, exp) | r (shift, shift) |
|---|---|---|---|---|---|---|---|
| primary (width<=1) | 6 | 0.72 | 0.58 | +0.14 | [-0.68, +0.71] | 0.98 | 0.56 |
| secondary (single bracket) | 75 | 1.59 | 0.74 | +0.84 | [+0.34, +1.34] | 0.83 | 0.02 (Spearman 0.16) |

**The model loses to the null (model pKa only)** on the secondary set, with a CI excluding 0; the primary set is too small to conclude. The shift
correlation is ~0: the model does not predict the sign or size of experimental shifts. The high Pearson r of midpoints is the
trivial between-type spread inherited from the model pKa values. Per type (secondary) RMSE model/null: ASP 1.27/0.61, GLU 1.09/0.72, HIS
1.97/0.84, LYS 0.73/0.21, N-term 1.01/0.53, TYR 2.99/1.20. Experimental values lie inside the interval for 73% of the single-bracket
sites, which is weak evidence given mean width 3.4. Likely causes (not tested): eps=4r gives large Coulomb shifts without
desolvation or screening, deposited structures with ligands/hemes omitted, 293-309 K data against a 298.15 K model.

## 5. Reproduction (from the worktree root; PYTHONPATH=src)
```
python3 experiments/coupled_pka/build_specs.py            # PDB download, manifest, specs
python3 experiments/coupled_pka/validate_exact.py data/coupled_pka/results/validate_exact.json 4PTI SUB:2LZT SUB:3RN3 SUB:1STN SUB:2RN2
python3 experiments/coupled_pka/run_system.py 4HHB OUT.json --radii 0,4,6,8 --frac-phs 3,5,7,9,11 --bracket-radii 0,4,6,8
python3 experiments/coupled_pka/compare_mid_exact.py run_4HHB.json OUT.json
# 1OEL: run_system.py 1OEL ... --frac-phs 7 --part p --nparts 8 (and --radii 8 --frac-phs 3,5,9,11 --bracket-radii "" --budget 256); mc_large.py 1OEL OUT; analyze_large.py
python3 experiments/coupled_pka/compare_experiment.py OUT.json results/run_*.json
PYTHONPATH=src pytest tests/test_coupled_pka_core.py tests/test_coupled_pka_certs.py   # certs tests need numpy
```
Raw PDB/PKAD files are not committed (hashes only); `build_specs.py` fetches PDBs, the PKAD page must be re-downloaded to `data/coupled_pka/raw/pkad2/pkad2_index.html`.

## 6. Deviations from the protocol
Addendum 1 (large protein 1OEL instead of 4HHB, because 4HHB has treewidth 17). Large-system pH 3,5,9,11 certificates use budget
256 and R=8 only (pH 7 uses 1024). Exact validation used 20-site subsets with budget list as stated. Nothing was changed after seeing experimental comparisons.
**Post-hoc deviation (disclosed after adversarial review P2-M1):** the MC-agreement analysis in `analyze_large.py` added an
absolute slack of +1e-3 to the registered "3 standard errors" criterion. §4 now reports both numbers, with the registered one
as primary. The experimental (PKAD) comparison is unaffected and matches the registered analysis.

## 7. Limitations
Hard cutoff and eps=4r are modelling choices, not validated physics; no desolvation, screening, tautomers, ligands, or temperature;
interior capped at 9 sites; pKa grid 1..13 so extreme shifts yield no interval (59/134 experimental sites); bracket
resolution is limited by early stopping (undecided midpoints); the float proposer and the float VE/MC references are untrusted (only the Fraction checker is trusted, and
it checks certificates, not the structure parsing or spec generation); the PKAD values are taken as given with per-record references and were not
individually re-verified against the papers; model radii other than R=8 were not used for experimental midpoints; single-bracket certification gives existence, not uniqueness, of a crossing.
No claim is made about non-uniqueness of the model's Gibbs measure.
