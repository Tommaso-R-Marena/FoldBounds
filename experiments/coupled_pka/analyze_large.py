"""Large-system summary. usage: analyze_large.py SYSTEM OUT.json MC.json RESULT_GLOB..."""
import glob, json, sys
from collections import Counter
import numpy as np
from common import *
from foldbounds.coupled_pka.treewidth import spec_width

system, out, mcf = sys.argv[1:4]
files = sys.argv[4:]
sp = load_spec(system)
n = len(sp["sites"])
recs, secs, chk = {}, 0.0, Counter()
for pat in files:
    for f in sorted(glob.glob(pat)):
        d = json.load(open(f))
        secs += d["seconds"]; chk.update(d["check"])
        for r in d["records"]:
            e = recs.setdefault(r["i"], {"id": r["id"], "rtype": r["rtype"], "frac": {}, "pka": {}})
            e["frac"].update(r["frac"]); e["pka"].update(r["pka"])
assert len(recs) == n, (len(recs), n)
nbr = [dict() for _ in range(n)]
for k, v in sp["pairs_milli"].items():
    a, b = map(int, k.split(",")); nbr[a][b] = v; nbr[b][a] = v
strength = np.array([sum(nbr[i].values()) / 1000 for i in range(n)])
deg = np.array([len(nbr[i]) for i in range(n)])
res = {"system": system, "n_sites": n, "n_pairs": len(sp["pairs_milli"]), "min_fill_width_upper_bound": spec_width(sp),
       "cpu_seconds_total": secs, "exact_checker": dict(chk)}
w = {}
for R in (0, 4, 6, 8):
    w[R] = np.array([recs[i]["frac"][f"7.0|{R}"][1] - recs[i]["frac"][f"7.0|{R}"][0] for i in range(n)])
    res[f"theta_width_pH7_R{R}"] = {"mean": float(w[R].mean()), "median": float(np.median(w[R])), "p90": float(np.percentile(w[R], 90)),
                                    "frac_le_0.01": float((w[R] <= .01).mean()), "frac_le_0.1": float((w[R] <= .1).mean()),
                                    "frac_gt_0.5": float((w[R] > .5).mean())}
fail = w[8] > 0.1
res["fail_threshold"] = 0.1
res["n_fail_R8"] = int(fail.sum())
res["fail_by_type"] = dict(Counter(recs[i]["rtype"] for i in range(n) if fail[i]))
res["all_by_type"] = dict(Counter(recs[i]["rtype"] for i in range(n)))
res["strength_pK_mean_fail_vs_ok"] = [float(strength[fail].mean()) if fail.any() else None, float(strength[~fail].mean())]
res["degree_mean_fail_vs_ok"] = [float(deg[fail].mean()) if fail.any() else None, float(deg[~fail].mean())]
# failure rate by coupling-strength quartile
q = np.percentile(strength, [25, 50, 75])
res["fail_rate_by_strength_quartile"] = [float(fail[(lo <= strength) & (strength <= hi)].mean()) for lo, hi in
                                         zip([-1, q[0], q[1], q[2]], [q[0], q[1], q[2], 1e9])]
res["strength_quartile_edges"] = q.tolist()
# connected components of failing sites in the interaction graph
seen, comps = set(), []
for s in range(n):
    if fail[s] and s not in seen:
        st, c = [s], 0; seen.add(s)
        while st:
            v = st.pop(); c += 1
            for u in nbr[v]:
                if fail[u] and u not in seen:
                    seen.add(u); st.append(u)
        comps.append(c)
res["fail_components_sizes_top10"] = sorted(comps, reverse=True)[:10]
res["n_fail_components"] = len(comps)
kinds = Counter(recs[i]["pka"]["8"]["kind"] for i in range(n))
res["pka_kinds_R8"] = dict(kinds)
widths = [recs[i]["pka"]["8"]["hi"] - recs[i]["pka"]["8"]["lo"] for i in range(n) if recs[i]["pka"]["8"]["kind"] == "single"]
res["pka_single_width"] = {"n": len(widths), "mean": float(np.mean(widths)), "median": float(np.median(widths)),
                           "frac_le_1": float(np.mean(np.array(widths) <= 1.0)), "hist": dict(Counter(widths))}
# MC comparison
mc = json.load(open(mcf))
cmp = {}
for ph in mc["phs"]:
    key = f"{float(ph)}|8"
    m = np.array(mc["mean"][str(ph)]); se = np.array(mc["se"][str(ph)])
    lo = np.array([recs[i]["frac"][key][0] for i in range(n)]); hi = np.array([recs[i]["frac"][key][1] for i in range(n)])
    ins = (m >= lo) & (m <= hi)
    ins3 = (m >= lo - 3 * se - 1e-3) & (m <= hi + 3 * se + 1e-3)
    narrow = (hi - lo) <= 0.1
    over = np.maximum(lo - m, m - hi); over[ins] = 0
    cmp[str(ph)] = {"inside_strict": float(ins.mean()), "inside_3se_plus_1e-3": float(ins3.mean()), "n_outside_3se": int((~ins3).sum()),
                    "n_narrow": int(narrow.sum()), "inside_strict_among_narrow": float(ins[narrow].mean()) if narrow.any() else None,
                    "max_overshoot": float(over.max()), "mean_width": float((hi - lo).mean()),
                    "outside_3se_ids": [recs[i]["id"] for i in np.where(~ins3)[0]][:10]}
res["mc_vs_certificate_R8"] = cmp
json.dump(res, open(out, "w"), indent=1)
print(json.dumps(res, indent=1))
