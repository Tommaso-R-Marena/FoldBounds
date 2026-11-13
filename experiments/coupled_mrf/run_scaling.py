"""Experiment (i): size independence.  Ising grid, field perturbation at the centre.

For each L and radius R: coupled enclosure (exhaustive over 2^(4R) boundary states),
spanning-tree separate bounds, exact oracle (row DP, L <= ORACLE_MAX), the independent
checker (L <= CHECK_MAX).  Fixed parameters, no randomness.
"""
import json, sys, time
from fractions import Fraction as F
from common import write_csv, log10_frac
from foldbounds.coupled_mrf.generators import ising_grid, center
from foldbounds.coupled_mrf.certificate import exhaustive_enclosure, certificate_dict
from foldbounds.coupled_mrf.separator import radius_separator
from foldbounds.coupled_mrf.baselines import separate_ratio_interval
from foldbounds.coupled_mrf.oracles import grid_transfer_z
from foldbounds.coupled_mrf.io import model_to_dict
from foldbounds.coupled_mrf import checker

Q = F(3, 2)            # exp(2 beta J), beta = ln(1.5)/2 ~ 0.20 < beta_c ~ 0.44 (high temperature)
W = F(1)               # no bulk field
U = F(1, 2)            # B: extra field factor [1, 1/2] on the centre spin
LS = [4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 24, 32, 48, 64, 100, 128]
ORACLE_MAX = 18
CHECK_MAX = 24
RADII = [0, 1, 2, 3]
fields = ["L", "n_vars", "R", "boundary_states", "coupled_lower", "coupled_upper", "coupled_width",
          "coupled_total_s", "coupled_core_s", "sep_log10_lower", "sep_log10_upper", "sep_width_log10_ratio",
          "sep_s", "exact", "contains_exact", "oracle_s", "checker_ok", "checker_s", "oracle_value_float"]
rows = []
violations = 0
for L in LS:
    m = ising_grid(L, Q, W, ("field", center(L), U))
    t0 = time.perf_counter(); lo, hi, _ = separate_ratio_interval(m); sep_s = time.perf_counter() - t0
    exact = None; oracle_s = ""
    if L <= ORACLE_MAX:
        t0 = time.perf_counter(); exact = grid_transfer_z(m, "B") / grid_transfer_z(m, "A"); oracle_s = time.perf_counter() - t0
    for R in RADII:
        t0 = time.perf_counter()
        sp = radius_separator(m, R)
        enc = exhaustive_enclosure(m, sp)
        total = time.perf_counter() - t0
        row = dict(L=L, n_vars=L * L, R=R, boundary_states=enc.n_states, coupled_lower=float(enc.lower),
                   coupled_upper=float(enc.upper), coupled_width=float(enc.upper - enc.lower),
                   coupled_total_s=total, coupled_core_s=enc.seconds, sep_log10_lower=log10_frac(lo),
                   sep_log10_upper=log10_frac(hi), sep_width_log10_ratio=log10_frac(hi) - log10_frac(lo),
                   sep_s=sep_s, oracle_s=oracle_s)
        if exact is not None:
            ok = enc.contains(exact) and lo <= exact <= hi
            row.update(exact=float(exact), contains_exact=ok, oracle_value_float=float(exact))
            if not ok:
                violations += 1
                print("SOUNDNESS VIOLATION", L, R, file=sys.stderr)
        if L <= CHECK_MAX:
            md = json.loads(json.dumps(model_to_dict(m)))
            cert = json.loads(json.dumps(certificate_dict(m, sp, enc)))
            t0 = time.perf_counter(); ok, why = checker.check_certificate(md, cert)
            row.update(checker_ok=ok, checker_s=time.perf_counter() - t0)
            if not ok:
                violations += 1
                print("CHECKER REJECTED", L, R, why, file=sys.stderr)
        rows.append(row)
        print(L, R, row["coupled_width"], round(row["sep_width_log10_ratio"], 1), round(total, 3), flush=True)
write_csv("scaling_ising_field.csv", rows, fields)
print("violations:", violations)
