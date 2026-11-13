"""Pre-registered comparison (PROTOCOL section 6). usage: compare_experiment.py OUT.json RUN.json [RUN.json ...]"""
import json, re, sys
import numpy as np
from common import *
from foldbounds.coupled_pka import model as M

wild = pkad_wild()
runs = {}
for f in sys.argv[2:]:
    r = json.load(open(f)); runs.setdefault(r["system"], []).extend(r["records"])
TYPES = {"ASP", "GLU", "HIS", "LYS", "TYR", "CYS", "N-term", "C-term"}
rows, excluded = [], {"non-numeric": 0, "type": 0, "no-site": 0}
per_site = {}
for r in wild:
    pdb = r[1]
    if pdb not in PROTEINS:
        continue
    rt = r[2] if r[2] in ("N-term", "C-term") else r[2].upper()
    if rt not in TYPES:
        excluded["type"] += 1; continue
    if not re.fullmatch(r"\d+(\.\d+)?", r[5]):
        excluded["non-numeric"] += 1; continue
    per_site.setdefault((pdb, r[3], int(r[4]), rt), []).append((float(r[5]), r[19], r[18]))
for (pdb, ch, num, rt), v in sorted(per_site.items()):
    sp = load_spec(pdb)
    idx = [k for k, s in enumerate(sp["sites"]) if s["res"][0] == ch and s["res"][1] == num and s["rtype"] == rt]
    if not idx:
        excluded["no-site"] += 1; continue
    i = idx[0]
    rec = next(x for x in runs[pdb] if x["i"] == i)
    p = rec["pka"]["8"]
    rows.append({"pdb": pdb, "chain": ch, "res": num, "type": rt, "exp": float(np.mean([x[0] for x in v])), "n_rec": len(v),
                 "temp": v[0][2], "ref": v[0][1], "null": float(M.MODEL_PKA[rt]), "kind": p["kind"], "lo": p["lo"], "hi": p["hi"],
                 "mid": None if p["lo"] is None else (p["lo"] + p["hi"]) / 2,
                 "width": None if p["lo"] is None else p["hi"] - p["lo"]})

def rmse(x): return float(np.sqrt(np.mean(np.square(x))))
def pear(a, b):
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 and np.std(a) > 0 and np.std(b) > 0 else float("nan")
def rank(a):
    o = np.argsort(a, kind="stable"); r = np.empty(len(a)); r[o] = np.arange(len(a))
    for v in np.unique(a):
        m = a == v; r[m] = r[m].mean()
    return r
def stats(sel):
    if len(sel) < 3:
        return {"n": len(sel)}
    e = np.array([x["exp"] for x in sel]); m = np.array([x["mid"] for x in sel]); n0 = np.array([x["null"] for x in sel])
    se, sn = m - e, n0 - e
    return {"n": len(sel), "rmse_model": rmse(se), "rmse_null": rmse(sn), "delta": rmse(se) - rmse(sn),
            "pearson_mid_exp": pear(m, e), "pearson_shift": pear(m - n0, e - n0),
            "spearman_shift": pear(rank(m - n0), rank(e - n0)), "mean_width": float(np.mean([x["width"] for x in sel])),
            "coverage_exp_in_interval": float(np.mean([x["lo"] <= x["exp"] <= x["hi"] for x in sel]))}
def boot(sel, B=10000, seed=7):
    rng = np.random.default_rng(seed)
    prot = sorted({x["pdb"] for x in sel}); by = {p: [x for x in sel if x["pdb"] == p] for p in prot}
    ds = []
    for _ in range(B):
        pick = rng.integers(0, len(prot), len(prot))
        s = []
        for k in pick:
            g = by[prot[k]]; s += [g[j] for j in rng.integers(0, len(g), len(g))]
        e = np.array([x["exp"] for x in s]); m = np.array([x["mid"] for x in s]); n0 = np.array([x["null"] for x in s])
        ds.append(rmse(m - e) - rmse(n0 - e))
    return [float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5))]
prim = [x for x in rows if x["kind"] == "single" and x["width"] <= 1.0]
sec = [x for x in rows if x["kind"] == "single"]
out = {"n_rows": len(rows), "excluded": excluded,
       "kinds": {k: sum(1 for x in rows if x["kind"] == k) for k in {x["kind"] for x in rows}},
       "primary": stats(prim), "secondary": stats(sec), "rows": rows}
for nm, sel in (("primary", prim), ("secondary", sec)):
    if len(sel) >= 3:
        out[nm]["delta_ci95_bootstrap"] = boot(sel)
out["by_type_secondary"] = {t: stats([x for x in sec if x["type"] == t]) for t in sorted({x["type"] for x in sec})}
out["by_protein_secondary"] = {p: stats([x for x in sec if x["pdb"] == p]) for p in sorted({x["pdb"] for x in sec})}
# null-model statistics for all rows (independent of the model), for context
out["null_rmse_all_rows"] = rmse([x["null"] - x["exp"] for x in rows])
json.dump(out, open(sys.argv[1], "w"), indent=1)
print(json.dumps({k: out[k] for k in ("n_rows", "excluded", "kinds", "primary", "secondary", "null_rmse_all_rows")}, indent=1))
