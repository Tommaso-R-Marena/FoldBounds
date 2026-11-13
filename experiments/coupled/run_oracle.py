"""Oracle experiment: ddF widths (trivial / separate / coupled T1) at matched frontiers vs exact.
usage: run_oracle.py N OUT.csv [CAP] [NSEQ]   (stdlib only; seed fixed)"""
import csv, math, random, sys, time
from fractions import Fraction as F
sys.path.insert(0, "src")
from foldbounds.coupled.hp_mutation import HPMutation, cert_json
from foldbounds.coupled.hp_oracle import Oracle
from foldbounds.coupled.checker import check

QS = [F(3, 2), F(2), F(3)]
CONFIGS = [("linear", "span"), ("linear", "ratio"), ("centred", "ratio")]


def ln(x):
    return math.log(x.numerator) - math.log(x.denominator)


def trie_nodes(words):
    pre = {w[:t] for w in words for t in range(1, len(w))}
    return len(pre) + len(words)


def main():
    n, out = int(sys.argv[1]), sys.argv[2]
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 6000
    nseq = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    only = sys.argv[5].split(",") if len(sys.argv) > 5 else None   # e.g. linear-span (deep run: no bfs)
    rng = random.Random(20261113 + n)
    seqs = ["".join(rng.choice("HP") for _ in range(n)) for _ in range(nseq)]
    O = Oracle(n)
    assert O.check_count()
    total_nodes = sum(len(Oracle(m).walks) for m in range(3, min(n, 11) + 1)) if False else None
    A = [1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292, 324932, 881500, 2374444, 6416596]
    total_nodes = sum(A[m - 1] // 4 for m in range(2, n + 1))
    spec = {}
    def spectrum(s):
        if s not in spec:
            spec[s] = O.spectrum(s)
        return spec[s]
    cols = ["n", "seq", "i", "q", "mut", "config", "regions", "trie_frac", "exact_lo_ddF", "exact_hi_ddF",
            "triv_w", "sep_w", "coup_w", "exact_ddF", "sound_triv", "sound_sep", "sound_coup",
            "checker", "cert_bytes", "t_sec"]
    f = open(out, "w", newline=""); W = csv.writer(f); W.writerow(cols)
    rowid = 0
    t_all = time.time()
    for s in seqs:
        for i in range(n):
            P0 = HPMutation(s, i, F(2))
            sa, sb = spectrum(s), spectrum(P0.B)
            for q in QS:
                ZA = sum((q ** k * c for k, c in sa.items()), F(0))
                ZB = sum((q ** k * c for k, c in sb.items()), F(0))
                exact = ZB / ZA
                rowid += 1
                for order, kind in (CONFIGS + [("linear", "bfs")]):
                    if only and (order + "-" + kind) not in only:
                        continue
                    P = HPMutation(s, i, q, order)
                    tl, th = P.trivial()
                    t0 = time.time()
                    def emit(words, cfg, verify):
                        lo, hi = P.coupled(words); sl, sh = P.separate(words)
                        ok = ""; nb = ""
                        if verify:
                            c = P.certificate(words, (lo, hi)); js = cert_json(c)
                            ok = int(check(js)[0]); nb = len(js)
                        W.writerow([n, s, i, str(q), s[i] + ">" + P.B[i], cfg, len(words),
                                    f"{trie_nodes(words) / total_nodes:.6f}", f"{-ln(hi):.6f}", f"{-ln(lo):.6f}",
                                    f"{ln(th / tl):.6f}", f"{ln(sh / sl):.6f}" if sh else "inf", f"{ln(hi / lo):.6f}",
                                    f"{-ln(exact):.6f}", int(tl <= exact <= th),
                                    int(sl <= exact <= (sh if sh else exact)), int(lo <= exact <= hi),
                                    ok, nb, f"{time.time() - t0:.3f}"])
                    if kind == "bfs":
                        for m in range(3, n):
                            words = P.bfs_frontier(m)
                            if len(words) > cap:
                                break
                            emit(words, "bfs", m in (5, 7))
                        continue
                    nxt, first = 4, True
                    cfg = f"{order}-{kind}"
                    last = None; did150 = False
                    for words in P.refine(kind=kind):
                        last = words
                        if len(words) >= nxt or first:
                            if not first:
                                v150 = (not did150) and len(words) >= 150
                                did150 = did150 or v150
                                emit(words, cfg, v150)
                            first = False
                            nxt = max(nxt + 1, int(nxt * 1.3))
                        if len(words) >= cap:
                            break
                    emit(last, cfg, rowid % 4 == 0)
        print(n, s, "done", f"{time.time() - t_all:.0f}s", file=sys.stderr, flush=True)
        f.flush()
    f.close()


main()
