"""Experiment (ii): phase-transition obstruction.  Width vs separator radius R.

Zero-field and small-field 2D ferromagnetic Ising, J = 1.  exp(2 beta_c) = 1 + sqrt 2
exactly, so q = exp(2 beta J) = (1 + sqrt 2)^(beta/beta_c); a rational approximation
(denominator <= 1000) defines the model.  Perturbation: B has an extra field factor
[1, 1/2] on the centre spin.  Window: L = 2R + 1 grid, diamond separator of radius R.

R <= 4: exhaustive over all 2^(4R) boundary states.  R > 4: the extremes of rho are
evaluated at the two uniform boundary states (all-0, all-1), justified by FKG
monotonicity (ferromagnet, nonnegative field: the centre marginal is monotone in the
boundary condition, and rho is monotone in that marginal).  The two methods are
cross-checked exactly at R <= 4.
"""
import sys, time
from fractions import Fraction as F
from common import write_csv
from foldbounds.coupled_mrf.generators import ising_grid, ising_q, center, BETA_C
from foldbounds.coupled_mrf.certificate import exhaustive_enclosure, local_ratio
from foldbounds.coupled_mrf.separator import radius_separator

RATIOS = [0.5, 0.8, 1.0, 1.25, 2.0]
FIELDS = [("zero", F(1)), ("small", F(51, 50))]
U = F(1, 2)
RMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
fields = ["beta_over_beta_c", "q", "q_float", "field", "w", "R", "boundary_vars", "method", "lower", "upper",
          "width", "seconds"]
rows = []
for ratio in RATIOS:
    q = ising_q(ratio)
    for fname, w in FIELDS:
        for R in range(0, RMAX + 1):
            L = 2 * R + 1 if R > 0 else 3
            m = ising_grid(L, q, w, ("field", center(L), U))
            sep = radius_separator(m, R)
            t0 = time.perf_counter()
            nb = len(sep.boundary)
            ends = []
            for st in ([0] * nb, [1] * nb):
                a, b = local_ratio(m, sep, st)
                ends.append(b / a)
            lo, hi = min(ends), max(ends)
            method = "fkg-uniform-boundary"
            if R <= 4:
                enc = exhaustive_enclosure(m, sep)
                assert (enc.lower, enc.upper) == (lo, hi), ("FKG reduction disagrees with exhaustive", ratio, fname, R)
                method = "exhaustive(+fkg-checked)"
            sec = time.perf_counter() - t0
            rows.append(dict(beta_over_beta_c=ratio, q=str(q), q_float=float(q), field=fname, w=str(w), R=R,
                             boundary_vars=nb, method=method, lower=float(lo), upper=float(hi),
                             width=float(hi - lo), seconds=sec))
            print(ratio, fname, R, float(hi - lo), round(sec, 2), flush=True)
write_csv("phase_ising.csv", rows, fields)
