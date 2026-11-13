"""4HHB (treewidth 17): compare certificates with a float64 variable-elimination reference.
usage: compare_mid_exact.py RUN.json OUT.json"""
import json, sys, math
from common import *
from foldbounds.coupled_pka import exact as ex
run = json.load(open(sys.argv[1]))
sp = load_spec(run["system"])
n = len(sp["sites"])
phs = sorted({float(k.split("|")[0]) for r in run["records"] for k in r["frac"]})
ref = {ph: ex.ve_theta(sp, ph) for ph in phs}
viol, tot, worst = 0, 0, 0.0
byR = {}
for r in run["records"]:
    for k, (lo, hi, nS, B, nl, nu) in r["frac"].items():
        ph, R = k.split("|"); ph = float(ph); R = int(R)
        x = ref[ph][r["i"]]
        tot += 1
        over = max(lo - x, x - hi, 0.0)
        worst = max(worst, over)
        viol += over > 1e-9
        byR.setdefault(R, []).append(hi - lo)
# pKa brackets: endpoints signs and monotone/crossing count on a 0.5 grid
grid = [1 + 0.5 * k for k in range(25)]
refg = {g: ex.ve_theta(sp, g) for g in grid}
bad_br, nbr_checked, multi_exact = 0, 0, []
for r in run["records"]:
    i = r["i"]
    th = [refg[g][i] for g in grid]
    cross = sum(1 for a, b in zip(th, th[1:]) if (a - .5) * (b - .5) < 0)
    if cross > 1:
        multi_exact.append(r["id"])
    for R, p in r["pka"].items():
        if p["kind"] in ("single", "multi"):
            nbr_checked += 1
            lo, hi = ex.ve_theta(sp, p["lo"], [i])[i], ex.ve_theta(sp, p["hi"], [i])[i]
            # a certified bracket needs theta>=1/2 at one end and <=1/2 at the other (for single)
            if p["kind"] == "single" and not ((lo >= .5 >= hi) or (lo <= .5 <= hi)):
                bad_br += 1
        if p["kind"].startswith("none"):
            pass
out = {"system": run["system"], "frac_checks": tot, "frac_violations": viol, "worst_overshoot": worst,
       "mean_theta_width_by_R": {R: sum(v) / len(v) for R, v in byR.items()},
       "bracket_checks": nbr_checked, "bracket_endpoint_violations": bad_br,
       "sites_with_multiple_exact_crossings_0.5grid": multi_exact, "cert_check": run["check"]}
json.dump(out, open(sys.argv[2], "w"), indent=1)
print(json.dumps(out, indent=1))
