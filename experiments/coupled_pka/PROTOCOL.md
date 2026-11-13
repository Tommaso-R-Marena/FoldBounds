# PROTOCOL: coupled multi-site protonation model with separator certificates

Frozen: 2026-11-13, BEFORE any comparison of model output with experimental pKa
values. Experimental values were only inspected for (i) which PDB entries exist,
(ii) record counts and which records are non-numeric (inclusion rules below).
No model output has been computed at the time of freezing.
Deviations after the freeze are listed in PKA_APPLICATION.md under "Deviations".

## 0. What is certified

The certificate concerns the MODEL (a pairwise Ising/MRF on binary protonation
states with the rational weights defined below), never experiment. Q_i =
P(sigma_i=1)/P(sigma_i=0) is enclosed by separator certificates; theta_i =
Q/(1+Q) is the protonation fraction.

## 1. Sites, states, charges

* Structure: PDB file from files.rcsb.org (URL, retrieval date, SHA-256 in
  `data/coupled_pka/MANIFEST.json`). Model 1 only. Altloc ' ' or 'A'. Protein
  ATOM records only; waters, ligands, ions, hemes are ignored (limitation).
* Chains: for every protein except 4HHB, the chains that carry at least one
  PKAD2 record for that entry; for 4HHB all four chains (A, B, C, D).
* sigma_i = 1 protonated, 0 deprotonated. Acid sites: ASP, GLU, TYR, CYS,
  C-terminus (charge q = sigma-1). Base sites: HIS, LYS, ARG, N-terminus
  (q = sigma). q = sigma - a_i with a_i = 1 (acid) or 0 (base).
* Charge centre (mean of listed atoms, Angstrom): ASP OD1,OD2; GLU OE1,OE2;
  HIS ND1,NE2; CYS SG; TYR OH; LYS NZ; ARG CZ; N-terminus N of the first
  residue of the chain; C-terminus mean of O and OXT of the last residue (O only
  if OXT is absent). Sites with a missing charge-centre atom are dropped and
  listed. A CYS whose SG lies within 2.5 A of another SG is a disulfide and is
  not a site.
* Termini of a chain = first and last residue with ATOM records in that chain.

## 2. Model pKa values (provenance)

From `propka/propka.cfg`, lines `model_pkas`, of PROPKA 3.5.1 (PyPI sdist
`propka-3.5.1.tar.gz`, SHA-256
d8083b58d8c4c553ff8704c706190fe37ff82253adcd6e53b31d8ad7d6d52257, re-downloaded
with `pip download propka==3.5.1 --no-binary :all:` and compared on 2026-11-13):
C- 3.20, ASP 3.80, GLU 4.50, HIS 6.50, CYS 9.00, TYR 10.00, LYS 10.50,
ARG 12.50, N+ 8.00. Only these nine values are used. They are the NULL MODEL.

## 3. Energy model

* T = 298.15 K. k_B = 0.001987204 kcal/mol/K, so kT ln10 = 1.364239... kcal/mol.
* Pair couplings: for sites i != j in different residues with centre distance
  r_ij <= R_c = 12.0 A (HARD cutoff), Coulomb energy with distance-dependent
  dielectric eps(r) = 4 r:
  E_ij = 332.0637 q_i q_j / (eps(r) r) = 83.01593 q_i q_j / r^2   kcal/mol.
  In pK units: w_ij = E_ij-coefficient / (kT ln10), quantised to integer
  thousandths: k_ij = round_half_even(1000 * 83.01593 / (r_ij^2 * kT ln10)).
  No screening, no desolvation or background terms, no self terms: the only
  protein effect is site-site Coulomb. Positive k_ij means like charges repel.
* State weight at pH p: W(sigma) = prod_i u_i(sigma_i) * prod_{i<j} phi_ij(sigma_i,sigma_j),
  u_i(0)=1, u_i(1)=10^(pKa_i - p), phi_ij = 10^(-(k_ij/1000) q_i q_j).
  In the sigma basis the coefficient of sigma_i sigma_j is +w_ij > 0 for every
  pair, i.e. the model is ANTIFERROMAGNETIC; the sign of q_i q_j only moves linear terms.
* Exact rationals. Every 10^(e/1000) with integer e is replaced by the rational
  R15(10^(e/1000)): evaluate with Python `decimal` at 60 digits, round to 15
  significant digits (ROUND_HALF_EVEN), convert to Fraction. THE MODEL IS DEFINED
  as the MRF with these rational weights (relative perturbation <= 5e-15 of
  the real-valued model; pKa are multiples of 1/100 and pH of 1/1000, so e is an
  integer). The certificate is exact for this rational MRF.
  The Monte Carlo baseline uses the float values of the same rationals.

## 4. Certificates

* Target: for site i at pH p, Q_i = Z_B/Z_A with A: sigma_i=0, B: sigma_i=1.
* Interior S_R(i) = {i} plus sites with centre distance <= R from site i,
  nearest first, capped at |S| <= 9. Hierarchy R in {0, 4, 6, 8} A (R=0: S={i};
  there the enclosure is closed-form and exact for the separator).
  Boundary dB = sites outside S with a nonzero coupling to S. Sites not in
  S or dB are the exterior (no factor touches both S and the exterior, because of the hard cutoff).
* Lower/upper enclosure on a cube C of boundary states (some boundary
  sites fixed, others free): regions xi = interior configuration of S\{i}; for
  region xi the A-mass d_xi(tau) = Psi(0,xi) prod_t h_t(0,xi,tau_t) is bracketed
  by [L,U] using, for each free t, min/max over tau_t, and the regional ratio
  r_xi(tau) = Psi(1,xi)/Psi(0,xi) prod_t h_t(1,xi,tau_t)/h_t(0,xi,tau_t) is
  bracketed by [r_min, r_max] likewise (this is separable in tau). Then
  Q = sum r d / sum d, and Theorem 1.1/1.2 of CEC_THEORY.md (T1) gives lambda_- (lower) and
  lambda_+ (upper) valid for all tau in C. The certified interval of Q is
  [min over leaves lambda_-, max over leaves lambda_+] over a partition of
  {0,1}^dB into cubes (decision tree, boundary variables split in a FIXED
  order: strongest coupling to S first). With every boundary variable fixed the
  bound is exact; hence exhaustive enumeration is the limit and the
  intermediate trees are sound relaxations. Leaf budgets: decisions 256 leaves per side,
  reported fractions 1024 leaves, exact-validation budgets {1, 8, 64, 2^14}.
* The float proposer builds trees; claimed bounds are rationals backed off by a relative 1e-9.
  The INDEPENDENT CHECKER (exact Fractions, own derivation of the weights from the
  integer spec) verifies: tree partition validity, and for every leaf the exact T1
  test g_-(lambda)>=0 (resp. g_+(lambda)<=0) at the claimed lambda.
* theta bounds: theta in [lam_-/(1+lam_-), lam_+/(1+lam_+)].
* Model pKa_{1/2}: the pH where Q=1. Per-site monotonicity of theta in pH is NOT
  assumed (antiferromagnetic couplings can give anomalous, non-monotone titration;
  only the total protonation is monotone). Procedure: integer pH grid 1..13;
  at each grid pH decide certifiably A (lambda_- > 1: protonated-majority), B
  (lambda_+ < 1) or undecided. Crossing brackets = adjacent decided grid points with
  different labels, skipping undecided ones; each is bisected 3 times (resolution
  0.125) with certificate decisions; an undecided midpoint stops bisection.
  If there is exactly one bracket the certified pKa interval is that
  bracket; with several, the hull and a flag `multi`; with none, no
  interval ('none-below', 'none-above' or 'undecided').
  Width = hi - lo. "Converged" = single bracket with width <= 1.0 pH unit.
  Under monotonicity of theta (checked exactly only on the small exact systems; see
  below) the bracket contains the unique crossing; in general it contains at least one crossing.
* Exact validation (≤ 22 sites): brute force by exact integer depth-first
  enumeration of all 2^n states, a different algorithm from the certificates.
  Systems: BPTI 4PTI (all sites) and four 18-20 site subsets of other proteins
  (nearest sites to a deterministic seed site, listed). Criterion: every exact
  Q_i and theta_i lies in every certified interval (all radii, all budgets,
  several pH); required violations = 0.
* Corrupted-certificate tests: tampered claimed lambda, dropped cube, duplicated
  cube, wrong variable order, tampered spec weights must all be rejected.

## 5. Large protein

4HHB (human deoxyhemoglobin, four chains), >= 100 sites expected. Min-fill
treewidth estimate on the cutoff interaction graph justifies that exact
inference is out of reach. Report: certified widths of theta at pH 7.0 for all
sites and radii, runtime, convergence with radius, sites that fail to converge
(width of theta > 0.1 at largest radius, or pKa interval not converged), pKa
brackets. Monte Carlo baseline: Gibbs sampling, 3 seeds (20261113, 20261114,
20261115), 64 parallel chains per seed, 1000 burn-in and 4000 sampling sweeps,
at pH 3, 5, 7, 9, 11; report the fraction of (site, pH) with MC mean inside the
certified theta interval (final radius, 1024 leaves) and inside the interval widened by 3 standard errors.

## 6. Pre-registered comparison with experiment

* Data: PKAD-2 (https://compbio.clemson.edu/pkad-2/ , retrieved 2026-11-13, SHA-256
  a71dfa05dc7d8122c72de730c66d412432f5a3db6d30b9ef0862e06e045f8123 of the page; the
  `wildData` array embedded in the page; record DOI/URL stored with each record).
* Proteins/structures: 4PTI (BPTI), 2LZT (HEWL), 3RN3 (RNase A), 1UBQ
  (ubiquitin), 1PGA (GB1), 1STN (SNase), 1A2P (barnase), 2RN2 (RNase H), 4HHB (haemoglobin).
* Inclusion: wild-type PKAD2 records of these PDB IDs whose residue type is one of ASP,
  GLU, HIS, LYS, TYR, CYS, N-term, C-term (case-normalised), whose pKa is a plain
  decimal number (no '<', '>', '~', lists), and whose site exists in the model.
  Several records for one site are averaged. No filtering on temperature or
  method (reported). Model is evaluated at T=298.15 K regardless of the
  measurement temperature (limitation).
* Predictors: (a) NULL: the PROPKA model pKa of the residue type; (b) MODEL:
  midpoint of the certified pKa interval at R=8 A. Primary set: sites with a converged
  interval (width <= 1.0). Secondary: all sites with any single-bracket interval.
  Sites without an interval are excluded from the secondary set and counted.
* Metrics (primary): RMSE(model) and RMSE(null) with the same sites; Delta = RMSE(model) - RMSE(null);
  Pearson r and Spearman rho between (midpoint - model pKa) and (exp - model pKa);
  Pearson r of midpoint vs exp; 95% bootstrap CI (10000 resamples, seed 7)
  by resampling proteins then sites, for Delta; per-type RMSEs. No parameter is tuned
  with experimental data. Results are reported whatever the sign. Expectation recorded in advance:
  with eps=4r Coulomb-only and no desolvation the model may well LOSE to the null.
