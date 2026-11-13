"""Review follow-up (adversarial review M2): the headline scaling instance had no bulk field, so
Z_B/Z_A = 3/4 exactly for every L by spin-flip symmetry, and "exact value inside the enclosure"
tested little. This companion run breaks the symmetry with a bulk field (W = 6/5 on every spin),
so the exact ratio varies with L. It checks containment against the exact row-DP oracle for
L <= ORACLE_MAX and records the exact values. Fixed parameters, no randomness."""
import sys, time
from fractions import Fraction as F
from common import write_csv
from foldbounds.coupled_mrf.generators import ising_grid, center
from foldbounds.coupled_mrf.certificate import exhaustive_enclosure
from foldbounds.coupled_mrf.separator import radius_separator
from foldbounds.coupled_mrf.oracles import grid_transfer_z

Q = F(3, 2)
W = F(6, 5)            # bulk field breaks the spin-flip symmetry
U = F(1, 2)
LS = [4, 5, 6, 7, 8, 10, 12, 14, 16]
RADII = [0, 1, 2, 3]
rows, violations = [], 0
for L in LS:
    m = ising_grid(L, Q, W, ("field", center(L), U))
    exact = grid_transfer_z(m, "B") / grid_transfer_z(m, "A")
    for R in RADII:
        enc = exhaustive_enclosure(m, radius_separator(m, R))
        ok = enc.contains(exact)
        violations += (not ok)
        rows.append(dict(L=L, R=R, exact=float(exact), exact_exact=str(exact)[:60],
                         lower=float(enc.lower), upper=float(enc.upper),
                         width=float(enc.upper - enc.lower), contains_exact=ok))
        print(L, R, f"{float(exact):.15f}", f"[{float(enc.lower):.6f},{float(enc.upper):.6f}]", ok, flush=True)
write_csv("scaling_ising_bulkfield.csv", rows, ["L", "R", "exact", "exact_exact", "lower", "upper", "width", "contains_exact"])
print("violations:", violations)
sys.exit(1 if violations else 0)
