"""Theorem 2.6 diagnostic: observed osc(rho')/osc(rho) vs Dobrushin coefficient delta(kappa)."""
from fractions import Fraction as F
from common import write_csv
from foldbounds.coupled_mrf.generators import ising_grid, ising_q, center
from foldbounds.coupled_mrf.certificate import annulus_contraction
from foldbounds.coupled_mrf.separator import radius_separator

rows = []
for ratio in (0.5, 1.0, 2.0):
    q = ising_q(ratio)
    L = 9
    m = ising_grid(L, q, F(1), ("field", center(L), F(1, 2)))
    for R in (1, 2):
        d = annulus_contraction(m, radius_separator(m, R), radius_separator(m, R + 1), exact_limit=3_000_000)
        rows.append(dict(beta_over_beta_c=ratio, inner_R=R, outer_R=R + 1, dobrushin=(None if d.dobrushin is None else float(d.dobrushin)),
                         exact=d.exact, osc_inner=float(d.osc_inner), osc_outer=float(d.osc_outer),
                         observed_ratio=d.observed_ratio, identity_ok=d.identity_ok, bound_ok=d.bound_ok))
        print(rows[-1], flush=True)
write_csv("contraction.csv", rows, list(rows[0]))
