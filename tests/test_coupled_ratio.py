import itertools, random, unittest
from fractions import Fraction as F
from foldbounds.coupled.ratio import *


def rnd(rng, kmax=6):
    K = rng.randint(1, kmax)
    regs = []
    for _ in range(K):
        mode = rng.random()
        L = F(rng.randint(0, 6), rng.randint(1, 4))
        U = L + (0 if mode < .2 else F(rng.randint(0, 8), rng.randint(1, 4)))
        r = F(rng.randint(0, 6), rng.randint(1, 3))
        R = r + (0 if rng.random() < .25 else F(rng.randint(0, 6), rng.randint(1, 3)))
        regs.append(Region(L, U, r, R))
    return regs


def brute(regs):
    lo = hi = None
    for ch in itertools.product((0, 1), repeat=len(regs)):
        a = [g.U if c else g.L for g, c in zip(regs, ch)]
        if sum(a) == 0:
            continue
        for key, which in (("r", 0), ("R", 1)):
            q = sum(getattr(g, key) * x for g, x in zip(regs, a)) / sum(a)
            if which == 0:
                lo = q if lo is None else min(lo, q)
            else:
                hi = q if hi is None else max(hi, q)
    return lo, hi


class RatioTests(unittest.TestCase):
    def test_bruteforce_random(self):
        rng = random.Random(12345)
        n = 0
        while n < 2500:
            regs = rnd(rng)
            if not admissible_nonempty(regs):
                with self.assertRaises(EmptyAdmissibleSet):
                    lower_enclosure(regs)
                continue
            n += 1
            lo, hi = brute(regs)
            el, eh = lower_enclosure(regs), upper_enclosure(regs)
            self.assertEqual(el.value, lo)
            self.assertEqual(eh.value, hi)
            self.assertEqual(g_minus(regs, lo), 0)
            self.assertEqual(g_plus(regs, hi), 0)
            self.assertGreaterEqual(g_minus(regs, lo), 0)
            self.assertLess(g_minus(regs, lo + F(1, 1000)), 0)   # Thm 1.1: lo+eps is not a valid lower bound
            self.assertGreater(g_plus(regs, hi - F(1, 1000)), 0)
            # argmin is admissible and attains the value
            for e, key in ((el, "r"), (eh, "R")):
                for g, a, b in zip(regs, e.alpha, e.beta):
                    self.assertTrue(g.L <= a <= g.U and g.r * a <= b <= g.R * a)
                self.assertEqual(objective(e.alpha, e.beta), e.value)
            # Theorem 1.4 dominance
            s = separate_lower(regs)
            self.assertGreaterEqual(lo, s)
            sp = separate_upper(regs)
            if sp is not None:
                self.assertLessEqual(hi, sp)
            # T1' dual certificates via the box polytope
            G, h = box_polytope(regs)
            yl = box_dual_certificate(regs, lo, "lower")
            self.assertTrue(verify_lp_lower_certificate(G, h, [g.r for g in regs], lo, yl))
            self.assertFalse(verify_lp_lower_certificate(G, h, [g.r for g in regs], lo + F(1, 100), box_dual_certificate(regs, lo + F(1, 100), "lower")))
            yu = box_dual_certificate(regs, hi, "upper")
            self.assertTrue(verify_lp_upper_certificate(G, h, [g.R for g in regs], hi, yu))
            self.assertFalse(verify_lp_upper_certificate(G, h, [g.R for g in regs], hi - F(1, 100), box_dual_certificate(regs, hi - F(1, 100), "upper")))

    def test_interior_points_inside(self):
        rng = random.Random(7)
        for _ in range(300):
            regs = rnd(rng)
            if not admissible_nonempty(regs):
                continue
            lo, hi = lower_enclosure(regs).value, upper_enclosure(regs).value
            for _ in range(5):
                a = [g.L + (g.U - g.L) * F(rng.randint(0, 4), 4) for g in regs]
                if sum(a) == 0:
                    continue
                b = [g.r * x + (g.R - g.r) * x * F(rng.randint(0, 4), 4) for g, x in zip(regs, a)]
                q = sum(b) / sum(a)
                self.assertTrue(lo <= q <= hi)

    def test_edge_cases(self):
        with self.assertRaises(EmptyAdmissibleSet):
            lower_enclosure([])
        with self.assertRaises(EmptyAdmissibleSet):
            upper_enclosure([Region(0, 0, 1, 2)])
        one = [Region(F(1), F(3), F(1, 2), F(2))]
        self.assertEqual(lower_enclosure(one).value, F(1, 2))      # Thm 1.5: gain unbounded
        self.assertEqual(separate_lower(one), F(1, 6))
        self.assertEqual(upper_enclosure(one).value, 2)
        self.assertEqual(separate_upper(one), 6)
        z = [Region(0, 5, 1, 2)]                                    # L all zero: S+ = inf
        self.assertIsNone(separate_upper(z))
        self.assertEqual(upper_enclosure(z).value, 2)
        u = [Region(1, 2, 3, 3), Region(0, 9, 3, 3)]                # uniform ratio (Thm 1.6)
        self.assertEqual((lower_enclosure(u).value, upper_enclosure(u).value), (3, 3))
        t = [Region(0, 1, 1, 1), Region(0, 1, 1, 1), Region(1, 1, 2, 2)]   # ties
        self.assertEqual(lower_enclosure(t).value, brute(t)[0])
        zeroU = [Region(0, 0, 0, 9), Region(2, 2, 1, 1)]            # inert region
        self.assertEqual(lower_enclosure(zeroU).value, 1)
        self.assertEqual(upper_enclosure(zeroU).value, 1)
        self.assertEqual(len(reciprocal_regions(zeroU)), 1)

    def test_exactness_rejections(self):
        for bad in (0.5, True, "1", None):
            with self.assertRaises(InvalidRegion):
                Region(bad, 1, 0, 1)
        with self.assertRaises(InvalidRegion):
            Region(2, 1, 0, 1)
        with self.assertRaises(InvalidRegion):
            Region(0, 1, 2, 1)
        with self.assertRaises(InvalidRegion):
            Region(-1, 1, 0, 1)

    def test_lp_verifier_rejects(self):
        G, h = [[F(1)]], [F(1)]
        self.assertTrue(verify_lp_lower_certificate(G, h, [F(1)], F(1), [F(0)]))
        self.assertFalse(verify_lp_lower_certificate(G, h, [F(1)], F(2), [F(0)]))   # reduced cost < 0
        self.assertFalse(verify_lp_lower_certificate(G, h, [F(1)], F(1), [F(-1)]))  # y < 0
        self.assertFalse(verify_lp_lower_certificate(G, h, [F(1)], F(0), [F(1)]))   # h^T y > 0
        with self.assertRaises(InvalidRegion):
            verify_lp_lower_certificate(G, h, [F(1)], 1.0, [F(0)])


if __name__ == "__main__":
    unittest.main()
