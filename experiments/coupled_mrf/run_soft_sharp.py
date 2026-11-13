"""Experiments (iii) sharpness and (iv) soft separator (Thm 2.7), all exact."""
import json, math, random, sys, time
from fractions import Fraction as F
from common import write_csv
from foldbounds.coupled_mrf import constructions as cons, checker
from foldbounds.coupled_mrf.certificate import (certificate_dict, exhaustive_enclosure, all_local_ratios)
from foldbounds.coupled_mrf.generators import ising_grid, center, random_sparse_binary
from foldbounds.coupled_mrf.io import model_to_dict
from foldbounds.coupled_mrf.model import CoupledModel, Factor, make_model
from foldbounds.coupled_mrf.oracles import grid_transfer_z, ve_z, brute_force_z
from foldbounds.coupled_mrf.separator import Separator, radius_separator

failures = 0


def local_submodel(m, sep):
    keep = [v for v in m.variables if v in sep.interior or v in set(sep.boundary)]
    ks = set(keep)
    f = lambda fs: tuple(x for x in fs if set(x.scope) <= ks)
    return CoupledModel(tuple(keep), {v: m.domains_a[v] for v in keep}, {v: m.domains_b[v] for v in keep},
                        f(m.shared), f(m.only_a), f(m.only_b), {})


# ------------------------------------------------------------ (iii) sharpness
rows3 = []
cases = []
for R in (1, 2):
    L = 2 * R + 1
    m = ising_grid(L, F(3, 2), F(6, 5), ("field", center(L), F(1, 2)))
    sep = radius_separator(m, R)
    cases.append((f"ising_R{R}", local_submodel(m, sep), Separator(sep.interior, sep.boundary, frozenset())))
m = random_sparse_binary(10, 3, seed=7)
sep = radius_separator(m, 1)
cases.append(("random_regular_R1", local_submodel(m, sep), Separator(sep.interior, sep.boundary, frozenset())))
for name, base, sep in cases:
    enc = exhaustive_enclosure(base, sep)
    for theta in (F(1), F(0), F(1, 2), F(1, 3)):
        model, sep2, lam = cons.sharp_exterior_model(base, sep, theta)
        ratio = ve_z(model, "B") / ve_z(model, "A")
        ok = ratio == lam
        if theta == 1:
            ok &= ratio == enc.lower
        if theta == 0:
            ok &= ratio == enc.upper
        cert = certificate_dict(model, sep2, exhaustive_enclosure(model, sep2))
        cok, _ = checker.check_certificate(json.loads(json.dumps(model_to_dict(model))), json.loads(json.dumps(cert)))
        ok &= cok and enc.contains(ratio)
        failures += (not ok)
        rows3.append(dict(case=name, boundary_vars=len(sep.boundary), theta=str(theta), enclosure_lower=str(enc.lower),
                          enclosure_upper=str(enc.upper), realised_ratio=str(ratio), lambda_target=str(lam),
                          exact_match=ok))
        print("sharp", name, theta, ok, flush=True)
write_csv("sharpness.csv", rows3, ["case", "boundary_vars", "theta", "enclosure_lower", "enclosure_upper",
                                    "realised_ratio", "lambda_target", "exact_match"])

# ------------------------------------------------------------ (iv) soft separator
rows4 = []
L = 7
m = ising_grid(L, F(3, 2), F(1), ("field", center(L), F(1, 2)))
exact = grid_transfer_z(m, "B") / grid_transfer_z(m, "A")
for R in (1, 2):
    sep = radius_separator(m, R)
    hard = exhaustive_enclosure(m, sep)
    pa, pb = all_local_ratios(m, sep)
    dims = [2] * len(sep.boundary)
    act = cons.soft_sensitivities(m, sep)
    for s in [0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16]:
        c = cons.soft_rational_c(s) if s > 0 else F(1)
        lo, hi, _, _ = cons.soft_enclosure_from_tables(pa, pb, dims, [c] * len(dims))
        rows4.append(dict(kind="width_vs_s", R=R, s=s, c=float(c), soft_lower=float(lo), soft_upper=float(hi),
                          soft_width=float(hi - lo), hard_width=float(hard.upper - hard.lower),
                          exact=float(exact)))
    lo, hi = cons.soft_enclosure(m, sep)
    ok = lo <= exact <= hi and hard.lower <= lo and hi <= hard.upper
    failures += (not ok)
    rows4.append(dict(kind="model_sensitivities", R=R, s="", c=";".join(f"{float(x):.4g}" for x in sorted(set(act))),
                      soft_lower=float(lo), soft_upper=float(hi), soft_width=float(hi - lo),
                      hard_width=float(hard.upper - hard.lower), exact=float(exact)))
    print("soft", R, ok, float(hi - lo), float(hard.upper - hard.lower), flush=True)

# exact soundness sweep with model-derived sensitivities (Ising grids; exterior in E_s by construction)
rng = random.Random(2026)
n_checked = n_bad = 0
for trial in range(24):
    L = rng.choice([5, 6, 7, 8])
    q = rng.choice([F(1, 2), F(2, 3), F(3, 2), F(2), F(4)])
    w = rng.choice([F(1), F(3, 4), F(5, 4)])
    u = rng.choice([F(1, 3), F(1, 2), F(2), F(3)])
    mm = ising_grid(L, q, w, ("field", center(L), u))
    ex = grid_transfer_z(mm, "B") / grid_transfer_z(mm, "A")
    for R in (1, 2):
        sep = radius_separator(mm, R)
        lo, hi = cons.soft_enclosure(mm, sep)
        hard = exhaustive_enclosure(mm, sep)
        ok = lo <= ex <= hi and hard.lower <= lo and hi <= hard.upper
        n_checked += 1
        n_bad += (not ok)
        rows4.append(dict(kind="soundness_sweep", R=R, s=f"L={L} q={q} w={w} u={u}", soft_lower=float(lo),
                          soft_upper=float(hi), soft_width=float(hi - lo), hard_width=float(hard.upper - hard.lower),
                          exact=float(ex), c="ok" if ok else "VIOLATION"))
failures += n_bad
print("soft soundness sweep:", n_checked, "checked,", n_bad, "violations")

# counterexample: exterior factor touching two boundary variables
c2 = F(2)
b = make_model({"b1": 2, "b2": 2, "o": 2}, shared=[Factor(("b1", "b2", "o"), tuple([c2, c2, F(1), F(1), F(1), F(1), c2, c2]))],
               only_b=[Factor(("b1", "b2"), (F(1), F(0), F(0), F(1)))])
sep = Separator.of([], ["b1", "b2"], ["o"])
exact_ce = ve_z(b, "B") / ve_z(b, "A")
hard = exhaustive_enclosure(b, sep)
try:
    cons.soft_sensitivities(b, sep)
    detected = False
except cons.SoftHypothesisError:
    detected = True
ce = cons.two_boundary_counterexample(c2)
single_flip = max(c2, F(1))
print("counterexample: exact", exact_ce, "soft enclosure", ce["soft_lo"], ce["soft_hi"], "hard", hard.lower, hard.upper,
      "hypothesis flagged:", detected)
ce_ok = (exact_ce > ce["soft_hi"]) and hard.lower <= exact_ce <= hard.upper and detected
rows4.append(dict(kind="counterexample_two_boundary_factor", R="", s="c=2 per variable", c="2", soft_lower=float(ce["soft_lo"]),
                  soft_upper=float(ce["soft_hi"]), soft_width=float(ce["soft_hi"] - ce["soft_lo"]),
                  hard_width=float(hard.upper - hard.lower), exact=float(exact_ce)))
failures += (not ce_ok)
write_csv("soft_separator.csv", rows4, ["kind", "R", "s", "c", "soft_lower", "soft_upper", "soft_width", "hard_width", "exact"])
print("FAILURES:", failures)
sys.exit(1 if failures else 0)
