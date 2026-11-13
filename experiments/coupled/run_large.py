"""n=30/50 without an oracle. Sound-by-construction: every row is verified by the independent checker;
additionally coupled must lie inside the trivial interval.  usage: run_large.py N OUT.csv"""
import csv, math, random, sys, time
from fractions import Fraction as F
sys.path.insert(0, "src")
from foldbounds.coupled.hp_mutation import HPMutation, cert_json
from foldbounds.coupled.checker import check

QS = [F(3, 2), F(2), F(3)]
CHECK = (100, 500, 1500)


def ln(x):
    return math.log(x.numerator) - math.log(x.denominator)


def main():
    n, out = int(sys.argv[1]), sys.argv[2]
    rng = random.Random(20261113 + n)
    seqs = ["".join(rng.choice("HP") for _ in range(n)) for _ in range(8)]
    W = csv.writer(open(out, "w", newline=""))
    W.writerow(["n", "seq", "i", "q", "mut", "config", "regions", "triv_w", "sep_w", "coup_w",
                "coup_in_triv", "checker", "cert_bytes", "t_sec"])
    t00 = time.time()
    for s in seqs:
        for i in sorted(rng.sample(range(n), 5)):
            for q in QS:
                for order, kind in (("linear", "ratio"), ("centred", "ratio"), ("linear", "span")):
                    P = HPMutation(s, i, q, order)
                    tl, th = P.trivial()
                    t0 = time.time()
                    todo = list(CHECK)
                    for words in P.refine(kind=kind):
                        if todo and len(words) >= todo[0]:
                            todo.pop(0)
                            lo, hi = P.coupled(words); sl, sh = P.separate(words)
                            js = cert_json(P.certificate(words, (lo, hi)))
                            W.writerow([n, s, i, str(q), s[i] + ">" + P.B[i], f"{order}-{kind}", len(words),
                                        f"{ln(th / tl):.6f}", f"{ln(sh / sl):.6f}", f"{ln(hi / lo):.6f}",
                                        int(tl <= lo and hi <= th), int(check(js)[0]), len(js), f"{time.time() - t0:.2f}"])
                            if not todo:
                                break
        print(n, s, f"{time.time() - t00:.0f}s", file=sys.stderr, flush=True)


main()
