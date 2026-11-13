"""Gibbs Monte Carlo baseline for a system. usage: mc_large.py SYSTEM OUT.json"""
import json, sys, time
import numpy as np
from common import *
from foldbounds.coupled_pka import mc
sp = load_spec(sys.argv[1])
res = {"system": sys.argv[1], "seeds": [20261113, 20261114, 20261115], "phs": [3, 5, 7, 9, 11], "mean": {}, "se": {}}
t0 = time.time()
for ph in res["phs"]:
    chains = []
    for seed in res["seeds"]:
        m, pc = mc.gibbs_theta(sp, float(ph), seed, chains=64, burn=1000, sweeps=4000)
        chains.append(pc)
    allc = np.concatenate(chains, axis=1)  # n x 192
    res["mean"][str(ph)] = allc.mean(1).tolist()
    res["se"][str(ph)] = (allc.std(1, ddof=1) / np.sqrt(allc.shape[1])).tolist()
    res.setdefault("seedmeans", {})[str(ph)] = [c.mean(1).tolist() for c in chains]
    print(ph, "%.0fs" % (time.time() - t0), flush=True)
json.dump(res, open(sys.argv[2], "w"))
