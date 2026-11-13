"""Certificates for every site of one system.
usage: run_system.py SYSTEM OUT.json --radii 0,4,6,8 --frac-phs 7 --bracket-radii 8 [--part k --nparts n] [--budget 1024] [--check-radii 8]"""
import argparse, json, time
from fractions import Fraction
from common import *
from foldbounds.coupled_pka import prop as P, checker as C

ap = argparse.ArgumentParser()
ap.add_argument("system"); ap.add_argument("out")
ap.add_argument("--radii", default="0,4,6,8")
ap.add_argument("--frac-phs", default="7")
ap.add_argument("--bracket-radii", default="8")
ap.add_argument("--check-radii", default="8")
ap.add_argument("--budget", type=int, default=1024)
ap.add_argument("--dec-budget", type=int, default=256)
ap.add_argument("--part", type=int, default=0); ap.add_argument("--nparts", type=int, default=1)
ap.add_argument("--limit", type=int, default=0)
a = ap.parse_args()
radii = [int(x) for x in a.radii.split(",")]
phs = [Fraction(x) for x in a.frac_phs.split(",")]
brad = [int(x) for x in a.bracket_radii.split(",") if x]
crad = {int(x) for x in a.check_radii.split(",") if x}
sp = load_spec(a.system)
nbr = P.adjacency(sp)
n = len(sp["sites"])
sites = [i for i in range(n) if i % a.nparts == a.part]
if a.limit:
    sites = sites[: a.limit]
res = {"system": a.system, "n": n, "records": [], "check": {"ok": 0, "bad": 0, "n": 0}, "args": vars(a)}
t0 = time.time()
tcheck = 0.0
for cnt, i in enumerate(sites):
    if cnt % 50 == 0:
        print(a.system, a.part, cnt, len(sites), '%.0fs' % (time.time() - t0), flush=True)
    rec = {"i": i, "id": sp["sites"][i]["id"], "rtype": sp["sites"][i]["rtype"], "frac": {}, "pka": {}}
    for R in sorted(set(radii) | set(brad)):
        S = P.interior(sp, i, R)
        prep = P.Prepared(sp, nbr, i, S, Fraction(7))
        if R in radii:
            for ph in phs:
                prep.set_ph(ph)
                lo = P.run_bound(prep, "lower", a.budget); up = P.run_bound(prep, "upper", a.budget)
                cl, cu = lo.claim(), up.claim()
                rec["frac"][f"{float(ph)}|{R}"] = [float(cl / (1 + cl)), float(cu / (1 + cu)), len(S), prep.B, lo.nleaf, up.nleaf]
                if R in crad:
                    tc = time.time()
                    for s, c in ((lo, cl), (up, cu)):
                        res["check"]["n"] += 1
                        try:
                            C.check(sp, s.certificate(sp, c, i, ph)); res["check"]["ok"] += 1
                        except C.CertificateError as e:
                            res["check"]["bad"] += 1; rec.setdefault("bad", []).append(str(e))
                    tcheck += time.time() - tc
        if R in brad:
            r, certs = P.pka_interval(sp, prep, a.dec_budget)
            rec["pka"][str(R)] = r
            if R in crad:
                tc = time.time()
                for c in certs:
                    res["check"]["n"] += 1
                    try:
                        C.check(sp, c); res["check"]["ok"] += 1
                    except C.CertificateError as e:
                        res["check"]["bad"] += 1; rec.setdefault("bad", []).append(str(e))
                tcheck += time.time() - tc
    res["records"].append(rec)
res["seconds"] = time.time() - t0
res["check_seconds"] = tcheck
json.dump(res, open(a.out, "w"))
print(a.system, len(sites), "sites", "%.1fs (check %.1fs)" % (res["seconds"], tcheck), res["check"])
