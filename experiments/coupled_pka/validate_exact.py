"""Exact validation on <=20-site systems: every exact Q_i must lie in every certified interval.
usage: validate_exact.py OUT.json SYSTEM [SYSTEM...]   (SYSTEM = 4PTI or SUB:<PDB>)"""
import json, math, sys, time
from fractions import Fraction
from common import *
from foldbounds.coupled_pka import model as M, exact as ex, prop as P, checker as C

PHS = [3, 5, 7, 9, 11]
RADII = [0, 4, 6, 8, 10]
BUDGETS = [1, 8, 64, 16384]


def subset(pdb, k=20):
    sp = load_spec(pdb)
    xyz = [s["xyz"] for s in sp["sites"]]
    seed = 0
    keep = sorted(sorted(range(len(xyz)), key=lambda j: (math.dist(xyz[seed], xyz[j]), j))[:k])
    return M.spec_sub(sp, keep, f"SUB:{pdb}")


def main():
    out, systems = sys.argv[1], sys.argv[2:]
    res = {"violations": 0, "checks": 0, "cert_rejected": 0, "systems": {}, "records": []}
    for name in systems:
        sp = load_spec(name) if not name.startswith("SUB:") else subset(name[4:])
        n = len(sp["sites"])
        nbr = P.adjacency(sp)
        t0 = time.time()
        info = {"n": n, "pairs": len(sp["pairs_milli"]), "hash": M.spec_hash(sp)[:16], "sites": [s["id"] for s in sp["sites"]]}
        # float monotonicity scan (not trusted)
        import numpy as np
        grid = [x / 4 for x in range(0, 57)]
        th = ex.float_marginals_all_ph(sp, grid)
        d = np.diff(th, axis=0)
        info["nonmonotone_sites"] = [sp["sites"][i]["id"] for i in range(n) if (d[:, i] > 1e-12).any()]
        info["max_increase"] = float(d.max())
        for ph in PHS:
            Q = ex.exact_marginals(sp, Fraction(ph))
            for i in range(n):
                for R in RADII:
                    S = P.interior(sp, i, R)
                    prep = P.Prepared(sp, nbr, i, S, Fraction(ph))
                    for bud in BUDGETS:
                        lo = P.run_bound(prep, "lower", bud)
                        up = P.run_bound(prep, "upper", bud)
                        cl, cu = lo.claim(), up.claim()
                        ok = True
                        for sd, s, c in (("lower", lo, cl), ("upper", up, cu)):
                            try:
                                C.check(sp, s.certificate(sp, c, i, Fraction(ph)))
                            except C.CertificateError:
                                res["cert_rejected"] += 1
                                ok = False
                        res["checks"] += 1
                        inside = ok and cl <= Q[i] <= cu
                        if not inside:
                            res["violations"] += 1
                        thlo = float(cl / (1 + cl)); thhi = float(cu / (1 + cu)); thx = float(Q[i] / (1 + Q[i]))
                        res["records"].append([name, ph, i, R, len(S), prep.B, bud, thlo, thhi, thx, lo.nleaf, up.nleaf, int(inside)])
        info["seconds"] = time.time() - t0
        res["systems"][name] = info
        print(name, info["n"], "nonmono", info["nonmonotone_sites"], "viol", res["violations"], "rej", res["cert_rejected"], "%.0fs" % info["seconds"], flush=True)
    json.dump(res, open(out, "w"))


main()
