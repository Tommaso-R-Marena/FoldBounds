import copy, json, random, unittest
from fractions import Fraction as F
from foldbounds.coupled.hp_oracle import Oracle
from foldbounds.coupled.hp_mutation import HPMutation, centred_order
from foldbounds.coupled.checker import check
from foldbounds.lattice import exact_spectrum

_O = {}
def oracle(n):
    if n not in _O:
        _O[n] = Oracle(n)
    return _O[n]


class OracleTests(unittest.TestCase):
    def test_counts(self):
        for n in range(3, 12):
            self.assertTrue(oracle(n).check_count(), n)

    def test_against_foldbounds_lattice(self):      # cross-check with the existing module (read-only)
        rng = random.Random(3)
        for _ in range(6):
            seq = "".join(rng.choice("HP") for _ in range(9))
            self.assertEqual(dict(sorted(oracle(9).spectrum(seq).items())), exact_spectrum(seq))


class ProducerCheckerTests(unittest.TestCase):
    def cases(self):
        rng = random.Random(11)
        for n in (7, 9, 10):
            for _ in range(3):
                seq = "".join(rng.choice("HP") for _ in range(n))
                yield n, seq, rng.randrange(n), rng.choice([F(3, 2), F(2), F(3)])

    def test_soundness_and_acceptance_all_orders(self):
        for n, seq, i, q in self.cases():
            O = oracle(n)
            for order in ("linear", "centred"):
                P = HPMutation(seq, i, q, order)
                exact = O.ratio(P.A, P.B, q)
                gen = P.refine()
                for step, words in enumerate(gen):
                    if step % 7 and step < 60:
                        continue
                    lo, hi = P.coupled(words)
                    slo, shi = P.separate(words)
                    self.assertTrue(lo <= exact <= hi, (seq, i, q, order))
                    self.assertTrue(slo <= exact <= shi)
                    self.assertTrue(slo <= lo and hi <= shi)        # Thm 1.4-style dominance (same data)
                    ok, why = check(P.certificate(words))
                    self.assertTrue(ok, why)
                    if step > 120:
                        break
                # complete exploration is exact
                words = P.bfs_frontier(n)
                lo, hi = P.coupled(words)
                self.assertEqual((lo, hi), (exact, exact))
                self.assertTrue(check(P.certificate(words))[0])

    def test_checker_tightness_matches_producer(self):
        # the checker accepts the producer's sharp endpoints and rejects any strictly tighter claim
        for n, seq, i, q in self.cases():
            P = HPMutation(seq, i, q)
            words = P.bfs_frontier(5)
            lo, hi = P.coupled(words)
            eps = F(1, 10**9)
            c = P.certificate(words)
            self.assertTrue(check(c)[0])
            self.assertFalse(check(P.certificate(words, (lo + eps, hi)))[0])
            self.assertFalse(check(P.certificate(words, (lo, hi - eps)))[0])

    def test_bfs_cover_counts(self):
        for n in (8, 9):
            P = HPMutation("H" * n, 2, F(2))
            for m in range(2, n + 1):
                W = P.bfs_frontier(m)
                self.assertTrue(check(P.certificate(W))[0])


class CorruptionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P = HPMutation("HPHHPHHPPH", 3, F(2))
        cls.words = cls.P.bfs_frontier(6)
        cls.good = cls.P.certificate(cls.words)
        assert check(cls.good)[0]

    def bad(self, mutate, msg=""):
        c = copy.deepcopy(self.good)
        mutate(c)
        ok, why = check(c)
        self.assertFalse(ok, msg or why)
        return why

    def test_corruptions(self):
        fr = lambda c: c["frontier"]
        self.bad(lambda c: fr(c).pop(), "dropped node (coverage)")
        self.bad(lambda c: fr(c).append(dict(fr(c)[0])), "duplicate")
        self.bad(lambda c: fr(c).append({"w": fr(c)[0]["w"][:-1]}), "overlap: ancestor added")
        self.bad(lambda c: fr(c).__setitem__(0, {"w": fr(c)[0]["w"] + "R"}), "child replaces node")
        self.bad(lambda c: fr(c).__setitem__(0, {"w": "U" + fr(c)[0]["w"][1:]}), "first step not R")
        self.bad(lambda c: fr(c).__setitem__(0, {"w": "RUL" + "D" + fr(c)[0]["w"][4:]}), "self-intersecting")
        self.bad(lambda c: fr(c).__setitem__(0, {"w": "RX" + fr(c)[0]["w"][2:]}), "bad letter")
        self.bad(lambda c: c["claim"].__setitem__("lo", str(F(c["claim"]["lo"]) * 2)), "lo too high")
        self.bad(lambda c: c["claim"].__setitem__("hi", str(F(c["claim"]["hi"]) / 2)), "hi too low")
        self.bad(lambda c: c["claim"].__setitem__("lo", 0.1), "float claim")
        self.bad(lambda c: c.__setitem__("q", "1/2"), "q<1")
        self.bad(lambda c: c.__setitem__("q", 2.0), "float q")
        self.bad(lambda c: c.__setitem__("seqB", c["seqA"]), "no mutation")
        self.bad(lambda c: c.__setitem__("seqB", c["seqB"][:-1] + ("P" if c["seqB"][-1] == "H" else "H")), "two mutations")
        self.bad(lambda c: c.__setitem__("i", (c["i"] + 1) % 10), "wrong mutated index")
        self.bad(lambda c: c.__setitem__("order", [1, 0] + list(range(2, 10))[::-1]), "invalid order")
        self.bad(lambda c: c.pop("claim"), "missing claim")
        self.bad(lambda c: c.__setitem__("frontier", []), "empty frontier")
        self.bad(lambda c: c.__setitem__("frontier", "RR"), "frontier wrong type")
        self.assertFalse(check("{not json")[0])
        self.assertFalse(check(None)[0])

    def test_witness_corruption(self):
        c = copy.deepcopy(self.good)
        k = next(j for j, nd in enumerate(c["frontier"]) if "wit" in nd)
        nd = c["frontier"][k]
        nd2 = dict(nd); nd2["wit"] = nd["wit"][:-1]            # too short
        c["frontier"][k] = nd2
        self.assertFalse(check(c)[0])
        nd3 = dict(nd); nd3["wit"] = "R" + "R" * (len(nd["wit"]) - 1)   # does not extend / not SAW-valid
        c["frontier"][k] = nd3
        self.assertFalse(check(c)[0])

    def test_claim_one_epsilon_beyond_sharp(self):
        c = copy.deepcopy(self.good)
        lo = F(c["claim"]["lo"])
        c["claim"]["lo"] = str(lo + F(1, 10**12))
        self.assertFalse(check(c)[0])
        c = copy.deepcopy(self.good)
        c["claim"]["hi"] = str(F(c["claim"]["hi"]) - F(1, 10**12))
        self.assertFalse(check(c)[0])


if __name__ == "__main__":
    unittest.main()
