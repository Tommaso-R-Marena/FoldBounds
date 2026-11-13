# Coupled separator certificates for Markov random fields: scaling experiments

Status: computational study of Theorems 2.1-2.8 of `CEC_THEORY.md` §2 on Ising grids and small
random models. No novelty is claimed: the enclosures are an elementary convex-combination
argument (Thm 2.1), and the point tested here is only *cost scaling* of the coupled ratio
enclosure versus separately computed bounds. All numbers below were produced by the scripts in
`experiments/coupled_mrf/` with fixed parameters (no randomness except seeded generators).

## 1. Methods

All trusted arithmetic is exact (`fractions.Fraction` or Python integers); model data must be
rationals (floats are rejected on input). Code: `src/foldbounds/coupled_mrf/`.

* **Model** (`model.py`, `io.py`): two factor graphs A and B over the same variables, stored as
  shared / only-A / only-B factors. Domains may differ between A and B for variables touched
  only by system-specific factors (mutation).
* **Separator** (`separator.py`): partition into interior / boundary / exterior. Verified:
  (S1) no factor touches both interior and exterior; (S2) every system-specific factor lies in
  interior + boundary; (S3) a variable whose domain differs between A and B is interior. Shell
  separators are graph-distance shells around the perturbation support.
* **Exact local ratio** (`certificate.py`, `tables.py`): `Phi_X(tau)` for every boundary state by
  exact integer variable elimination of the interior (min-fill / min-size order chosen by
  simulation); `rho = Phi_B / Phi_A`. Enclosure `[min rho, max rho]` computed exhaustively
  (chunked by clamping leading boundary variables) for up to 2^18 boundary states. If some
  `tau` has `Phi_A = 0 < Phi_B` the upper end is `+inf`.
* **Hierarchy and diagnostic**: `shell_hierarchy`, `check_nesting` (Thm 2.4, verified exactly),
  `annulus_contraction` (exact identity `Phi' = Phi K` and the Dobrushin coefficient of Thm 2.6).
* **Certificate** (`certificate_dict`): JSON with model hash, separator, claimed
  `lower`/`upper`/`tight`, witness boundary states, optional elimination-order hint.
* **Independent checker** (`checker.py`): imports only `fractions` and `json` (enforced by a
  test), own parser, own pure-Python SHA-256 canonical hash, own separator verification, own
  elimination; recomputes every `Phi_A, Phi_B`, accepts only if `claim.lower <= min rho` and
  `claim.upper >= max rho` (equality when `tight`; `"inf"` required when an unbounded state
  exists), checks witnesses, rejects unknown keys, floats, bad hashes, and fails closed on any
  exception or when the boundary exceeds its state limit.
* **Oracles** (`oracles.py`, independent of the producer): brute force; dict-based variable
  elimination; row-by-row ("broken row") transfer DP for binary pairwise grid models. Cross-
  validated in `tests/test_coupled_mrf_core.py`. The grid DP is feasible to about L = 18 (L=16:
  ~5 s, L=18: ~55 s with rational weights; exponential in L).
* **Separate baseline** (`baselines.py`): spanning-tree bounds, off-tree factors at their
  minimum (lower) / maximum (upper), exact tree DP; ratio interval
  `[LB_B/UB_A, UB_B/LB_A]`. Soundness: each off-tree factor lies between its constant min and
  max and all factors are nonnegative, so `LB <= Z <= UB`; exact when the graph is a tree
  (tested). Not implemented: TRW two-tree upper bound (optional). **Caveat:** FPTAS /
  correlation-decay methods (Weitz 2006) give much better separate bounds in uniqueness regimes
  at poly(n) cost. This baseline is the cheap tree bound only; the claim tested is
  equal-cost size-independence of the coupled bound, not that separate methods are hopeless.
* **Soft separator** (`constructions.py`, Thm 2.7): `soft_enclosure_from_tables` computes
  `min/max R(F)` over vertex fields `F_v(a) in {1, c_v}` exactly (integers, progressive
  contraction, 3^n candidates for n binary boundary variables). `soft_sensitivities` derives
  rational `c_v >= e^{s_v}` from the exterior factors and **refuses** (raises) if an exterior
  factor touches two boundary variables. For a free parameter s, `soft_rational_c` rounds
  `e^s` up to a rational, which keeps the enclosure sound (the width is monotone in c).
* **Sharpness** (`sharp_exterior_model`, Thm 2.2): appends a hidden variable `h` with domain
  `|T|`, hard couplings `[h_j = sigma_j]` to every boundary variable and unary `psi(h)` supported
  on argmin/argmax states; the resulting full model's exact `Z_B/Z_A` is compared with the
  predicted `theta*min + (1-theta)*max`.

## 2. Soundness arguments (why each bound is valid)

* Coupled enclosure: Lemma 2.0 `Z_X = sum_tau Phi_X(tau) Psi(tau)`, so for `Z_A > 0` and no
  `Phi_A = 0 < Phi_B`, `Z_B/Z_A` is a convex combination of `rho(tau)` with weights
  `Phi_A Psi / Z_A` (Thm 2.1). Shared factors inside the boundary are put on the local side;
  this only restricts the support and keeps the argument valid.
* Nesting: shared kernel argument (Thm 2.4); checked exactly in the tests/experiments (no
  violation found).
* Soft separator: Thm 2.7 (convex mixture over exterior states; per-coordinate linear-fractional
  monotonicity). Requires every exterior factor to touch at most one boundary variable.
  For Ising, each exterior neighbour edge contributes a factor `q` (or `1/q`) to `c_v`.
* Tree baseline: see §1.
* FKG shortcut (experiment ii only, R > 4): in a ferromagnetic Ising model with nonnegative
  field the centre marginal is monotone in the boundary condition (FKG/Holley) and `rho` is
  monotone in that marginal, so the extremes over boundary states are at the all-0 and all-1
  boundaries. This is a classical theorem, not enumeration; for R <= 4 it is checked against
  full exhaustive enumeration in every run (assertion in the script; 0 failures).

## 3. Results

Reproduction commands are in §6. Raw data: `experiments/coupled_mrf/results/*.csv`.

### (i) Size independence (Ising grid, q = 3/2, no bulk field, field factor [1, 1/2] on the centre)

Coupled width = `max rho - min rho` (radius R, exhaustive over `2^(4R)` boundary states).
Separate width is reported as `log10(UB_B/LB_A) - log10(LB_B/UB_A)`; the tree bounds are vacuous
(a ratio interval roughly `[10^-x, 10^x]`).

| L | coupled width R=2 | R=3 | separate log10 width | coupled time R=2 (s) | tree-bound time (s) | oracle (s) |
|---|---|---|---|---|---|---|
| 5 | 0.2052 | 0.0936* | 5.6 | 0.0015 | 0.002 | 0.001 |
| 8 | 0.2052 | 0.1059 | 17.3 | 0.0017 | 0.006 | 0.010 |
| 12 | 0.2052 | 0.1059 | 42.6 | 0.0021 | 0.016 | 0.40 |
| 18 | 0.2052 | 0.1059 | 101.8 | 0.0035 | 0.032 | 54.6 |
| 32 | 0.2052 | 0.1059 | 338.4 | 0.0070 | 0.108 | n/a |
| 64 | 0.2052 | 0.1059 | 1397.8 | 0.0292 | 0.534 | n/a |
| 128 | 0.2052 | 0.1059 | 5680.4 | 0.1361 | 3.30 | n/a |

(*at R=3 the widths for L <= 6 are smaller (0.048, 0.094, 0.0998 for L=4,5,6) because the shell is truncated by the grid edge. For all
L >= 7 the coupled width at each R is identical to 16 digits: R=0: 0.5, R=1: 0.3351, R=2: 0.2052,
R=3: 0.1059.) The coupled time shown is wall time for the whole pipeline including
separator verification, which is O(n) in the number of factors; the interior elimination itself
is independent of L. Hence the measured growth with L (0.0015 s -> 0.136 s) is entirely
Python overhead of linear scans of the model, whereas the *enclosure and its width* do not depend
on L. The exact value lay inside every coupled and separate enclosure for all L <= 18 (4 radii
x 11 sizes = 44 checks) and the independent checker accepted every certificate for L <= 24 (0
soundness failures, 0 checker rejections). Separate width grows like L^2 (log10 width ~ 0.35 L^2):
the spanning-tree bounds are vacuous at the sizes where the coupled bound is still sharp.

### (ii) Phase-transition obstruction (zero-field 2D Ising, J=1)

`q = round((1+sqrt2)^(beta/beta_c))`, rational with denominator <= 1000 . The exact rationals used are
q = 679/437 (0.5), 1010/499 (0.8), 2378/985 (1.0), 1613/536 (1.25), 5741/985 (2.0)
(also in the `q` column of `phase_ising.csv`; the approximation defines the model). Width of the hard-separator enclosure of the centre-field perturbation
`[1,1/2]` vs diamond radius R (R <= 4 exhaustive; R = 5..10 by FKG extremes):

| beta/beta_c | R=1 | R=2 | R=3 | R=4 | R=6 | R=8 | R=10 |
|---|---|---|---|---|---|---|---|
| 0.5 | 0.354 | 0.233 | 0.132 | 0.068 | 0.016 | 0.0035 | 0.0007 |
| 0.8 | 0.444 | 0.398 | 0.351 | 0.305 | 0.221 | 0.151 | 0.099 |
| 1.0 | 0.471 | 0.453 | 0.438 | 0.425 | 0.407 | 0.394 | 0.384 |
| 1.25 | 0.488 | 0.483 | 0.480 | 0.479 | 0.4775 | 0.4773 | 0.4772 |
| 2.0 | 0.4991 | 0.4990 | 0.4990 | 0.4990 | 0.4990 | 0.4990 | 0.4990 |

Matches Thm 2.8 predictions: exponential decay above T_c (factor ~ 2 per shell at 0.5; the rate
slows approaching T_c), slow decay at T_c (width 0.38 at R=10: power law, not resolved to
a fitted exponent here), and a plateau `(1-u) m*(beta)` below T_c (0.4772 at 1.25, 0.4990 at 2).
The "small field" repetition (bulk field factor `[1, 51/50]` on every site) is almost
indistinguishable for R <= 10: the field only starts to cut the plateau at larger R (0.8: 0.0956
vs 0.0988 at R=10; 1.0: 0.381 vs 0.384). A weak field does not restore decay within the radii
reachable here. Largest-R timing: about 2 s per (beta, R=10) point.

### (iii) Sharpness (Thm 2.2)

For the Ising local data (R=1: 4 boundary variables; R=2: 8) and a random 3-regular binary model
(R=1), the exterior with theta in {1, 0, 1/2, 1/3} gives full-model `Z_B/Z_A` exactly equal to the
min, the max, and the intermediate convex combinations (exact rational equality, 12 of 12 cases;
`sharpness.csv`). The independent checker also accepts the certificates of these extended models.

### (iv) Soft separator (Thm 2.7)

Ising L=7, centre field [1, 1/2]; exact ratio 0.75. Width (uniform `c = e^s`) vs s:

| s | 0 | 0.1 | 0.5 | 1 | 2 | 4 | 8 | 16 | hard |
|---|---|---|---|---|---|---|---|---|---|
| R=1 width | 0 | 0.020 | 0.097 | 0.177 | 0.274 | 0.327 | 0.3349 | 0.3351 | 0.3351 |
| R=2 width | 0 | 0.013 | 0.061 | 0.111 | 0.169 | 0.200 | 0.2051 | 0.2052 | 0.2052 |

As s grows the soft width tends to the hard width (Remark of Thm 2.7), and s=0 gives the exact
value. With the model-derived sensitivities (`c_v = q^{k_v}`, k_v the number of exterior
neighbours; here 3.375 and 2.25) the soft widths are 0.2052 (R=1) and 0.1059 (R=2), which equal
the *hard* widths at R+1 to all printed digits. For the ferromagnetic Ising grid the soft
enclosure with the true sensitivities coincides with the next-shell hard enclosure; it gives no
extra power there (negative finding, by FKG the extreme vertex fields are the uniform ones).
Exact soundness sweep: 24 random Ising grids (L in 5..8, q in {1/2,2/3,3/2,2,4}, bulk field,
u in {1/3,1/2,2,3}) x R in {1,2}: all 48 enclosures contain the exact ratio and lie inside the hard
enclosure (0 violations).

**Hypothesis needed (counterexample).** Two binary boundary variables, `Phi_A = 1`, `Phi_B =
[1,0,0,1]` (B rewards agreement; here a system-specific pairwise factor) and one exterior factor
`g(b1,b2,o)` equal to `[[2,1],[1,2]]` on `(b1,b2)` (flipping either variable changes g by at most
2). Exact `Z_B/Z_A = 2/3`, but the vertex-field enclosure with `c=2` per variable is
`[4/9, 5/9]`. Using `g = [[1,2],[2,1]]` gives `1/3`, also outside. The hard enclosure `[0,1]`
contains both. `soft_sensitivities` raises `SoftHypothesisError` on this model.

## 4. Negative findings and limitations

* Local certificates do not converge in the ordered phase: width plateaus at `(1-u) m*` (Thm 2.8
  made concrete). Size independence is cost-per-width, not convergence to the exact value.
* Even above T_c at beta/beta_c = 0.8 the width is still 0.10 at R=10; reaching 1e-3 needs R in
  the 20s, i.e. boundaries far beyond exhaustive enumeration. Exhaustive enclosures are limited
  to 4R <= 18, i.e. R <= 4, for 2D shells. R > 4 relied on the FKG reduction, which is specific to
  ferromagnetic Ising and is not a general method; it cannot be checked by the minimal checker.
* The independent checker enumerates all boundary states with Fractions (default limit 2^18), so
  it is slow for large boundaries; it re-hashes the whole model (pure-Python SHA-256).
* Timing of the coupled pipeline grows with L here because separator verification scans all
  factors (O(n)); only the interior elimination and the enclosure are L-independent.
* The tree baseline is weak by design; a Weitz-type or TRW separate bound would be far tighter
  in the high-temperature regime and is not compared. The TRW two-tree bound was not implemented.
* Theorem 2.6: the Dobrushin coefficient was computed only for inner R = 1, 2 (kernel size
  limit); observed oscillation ratios: 0.66 (R 1->2, beta/beta_c = 0.5) with delta = 0.804; at
  beta/beta_c = 1 and 2, delta = 0.990 and 0.99999, i.e. nearly vacuous, matching the lack of
  decay there. The bound held in every computed case; the identity `Phi' = Phi K` held exactly.
* Contraction/diagnostic and soft-separator enumeration cost are exponential in the boundary
  size (soft: 3^|boundary| candidates for binary variables; run here for <= 8 boundary variables).
* Zero soundness failures were observed in any run; the grid DP oracle was cross-validated
  against brute force and VE on small models, and enclosures against the oracle for L <= 18.

## 5. Tests

`tests/test_coupled_mrf_core.py` (30 tests: oracles, separator verification, enclosure containment
on random/Ising/rotamer-with-domain-change models, nesting, finite termination, unbounded case,
contraction identity, tree bounds, sharpness, soft soundness and counterexample) and
`tests/test_coupled_mrf_checker.py` (26 tests, more than 10 corruption / fail-closed cases). No numpy needed.

## 6. Reproduction

```
cd /home/user/wt-mrf
PYTHONPATH=src python3 -m pytest tests/test_coupled_mrf_core.py tests/test_coupled_mrf_checker.py -q
cd experiments/coupled_mrf
PYTHONPATH=../../src python3 run_scaling.py        # ~80 s,  results/scaling_ising_field.csv
PYTHONPATH=../../src python3 run_phase.py 10       # ~50 s,  results/phase_ising.csv
PYTHONPATH=../../src python3 run_soft_sharp.py     # ~2 s,   sharpness.csv, soft_separator.csv
PYTHONPATH=../../src python3 run_contraction.py    # ~30 s,  contraction.csv
```
