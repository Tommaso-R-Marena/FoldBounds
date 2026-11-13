"""Exact brute-force oracle for the 2D HP model (n <= ~16), independent of foldbounds.lattice.

Enumerates self-avoiding walks from (0,0) with first step R (one representative per rotation
orbit; reflections not quotiented). Count = A001411(n-1)/4. Stores, per walk, the
bitmasks of its nonconsecutive adjacent site pairs, so any sequence/q is a cheap pass.
"""
from collections import Counter
from fractions import Fraction

A001411 = [1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292, 324932,
           881500, 2374444, 6416596, 17245332]   # a(0..16)


class Oracle:
    def __init__(self, n: int):
        self.n = n
        self.walks = []          # (word, tuple of pair bitmasks)
        dirs = (("R", 1, 0), ("U", 0, 1), ("L", -1, 0), ("D", 0, -1))
        occ = {(0, 0): 0, (1, 0): 1}
        word = ["R"]
        path = [(0, 0), (1, 0)]
        pairs = []

        def rec():
            m = len(path)
            if m == n:
                self.walks.append(("".join(word), tuple(pairs)))
                return
            x, y = path[-1]
            for ch, dx, dy in dirs:
                p = (x + dx, y + dy)
                if p in occ:
                    continue
                added = 0
                for ex, ey in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    b = occ.get((p[0] + ex, p[1] + ey))
                    if b is not None and b != m - 1:
                        pairs.append((1 << b) | (1 << m)); added += 1
                occ[p] = m; path.append(p); word.append(ch)
                rec()
                word.pop(); path.pop(); del occ[p]
                for _ in range(added):
                    pairs.pop()

        if n >= 3:
            rec()

    def check_count(self) -> bool:
        return len(self.walks) * 4 == A001411[self.n - 1]

    def spectrum(self, seq: str) -> Counter:
        h = sum(1 << k for k, c in enumerate(seq) if c == "H")
        return Counter(sum((pm & h) == pm for pm in pms) for _, pms in self.walks)

    def Z(self, seq: str, q: Fraction) -> Fraction:
        return sum((q ** k * c for k, c in self.spectrum(seq).items()), Fraction(0))

    def ratio(self, seqA: str, seqB: str, q: Fraction) -> Fraction:
        return self.Z(seqB, q) / self.Z(seqA, q)
