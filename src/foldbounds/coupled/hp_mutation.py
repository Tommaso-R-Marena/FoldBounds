"""UNTRUSTED producer of T1 certificates for HP point mutations (see checker.py).

Nothing here is trusted: ``checker.check`` re-derives everything from the JSON.
Growth orders: ``linear`` (0,1,2,...) or ``centred`` (i, i+1, i-1, i+2, ... root at the
mutated residue). The ratio Z_B/Z_A is invariant to the choice because Z is a sum over
conformations modulo translation+rotation and the weights q^contacts are invariant;
fixing the first bond to R picks one representative per C4-orbit for any order, so
different orders enumerate different representatives of the same orbit set.
"""
from __future__ import annotations

import heapq
import json
import math
from fractions import Fraction

from .ratio import Region, lower_enclosure, upper_enclosure

STEP = {"R": (1, 0), "U": (0, 1), "L": (-1, 0), "D": (0, -1)}


def centred_order(n: int, i: int) -> list[int]:
    """i, i+1, i-1, i+2, i-2, ...: alternate the two arms, then finish the longer one."""
    order, lo, hi, right = [i], i, i, True
    while len(order) < n:
        if (right and hi + 1 < n) or lo == 0:
            hi += 1; order.append(hi)
        else:
            lo -= 1; order.append(lo)
        right = not right
    return order


def lnfrac(x: Fraction) -> float:
    return math.log(x.numerator) - math.log(x.denominator)


class Node:
    __slots__ = ("word", "L", "U", "r", "R", "LB", "UB", "wit", "dead", "complete")

    def region(self) -> Region:
        return Region(self.L, self.U, self.r, self.R)


class HPMutation:
    def __init__(self, seqA: str, i: int, q: Fraction, order: str | list[int] = "linear"):
        self.A = seqA
        self.n = n = len(seqA)
        self.i = i
        self.B = seqA[:i] + ("P" if seqA[i] == "H" else "H") + seqA[i + 1:]
        self.q = Fraction(q)
        self.sgn = -1 if seqA[i] == "H" else 1
        if order == "linear":
            self.order = list(range(n))
        elif order == "centred":
            self.order = centred_order(n, i)
        else:
            self.order = list(order)
        self.hpart = [k for k in range(n) if seqA[k] == "H" and k != i]
        self._cache: dict[str, Node] = {}

    # --- geometry
    def place(self, word: str):
        o = self.order
        pos = {o[0]: (0, 0)}
        occ = {(0, 0)}
        for t, ch in enumerate(word, start=1):
            s = o[t]
            c = s - 1 if (s - 1) in pos else s + 1
            d = STEP[ch]
            p = (pos[c][0] + d[0], pos[c][1] + d[1])
            if p in occ:
                return None
            occ.add(p)
            pos[s] = p
        return pos

    def child_letters(self, pos, m):
        s = self.order[m]
        c = s - 1 if (s - 1) in pos else s + 1
        occ = set(pos.values())
        x, y = pos[c]
        return [ch for ch, (dx, dy) in STEP.items() if (x + dx, y + dy) not in occ]

    @staticmethod
    def contacts(pos, seq):
        inv = {p: k for k, p in pos.items()}
        c = 0
        for a, p in pos.items():
            if seq[a] != "H":
                continue
            for dx, dy in STEP.values():
                b = inv.get((p[0] + dx, p[1] + dy))
                if b is not None and b > a + 1 and seq[b] == "H":
                    c += 1
        return c

    def ucap(self, seq, pos, m):
        placed = set(pos)
        n, tot = self.n, 0
        for t in range(m, n):
            j = self.order[t]
            if seq[j] == "H":
                before = placed | set(self.order[m:t])
                nj = sum(1 for a in before if seq[a] == "H" and abs(a - j) >= 3 and (a - j) % 2)
                tot += min(nj, 3 if j in (0, n - 1) else 2)
        return tot

    def witness(self, word, pos, limit=4000):
        """Greedy DFS for one complete self-avoiding extension (higher A-contacts first)."""
        n, o = self.n, self.order
        budget = [limit]
        occ = set(pos.values())
        inv = {p: k for k, p in pos.items()}

        def rec(m, out):
            if m == n:
                return True
            if budget[0] <= 0:
                return False
            budget[0] -= 1
            s = o[m]
            c = s - 1 if (s - 1) in pos else s + 1
            x, y = pos[c]
            opts = []
            for ch, (dx, dy) in STEP.items():
                p = (x + dx, y + dy)
                if p in occ:
                    continue
                gain = 0
                if self.A[s] == "H":
                    for ex, ey in STEP.values():
                        b = inv.get((p[0] + ex, p[1] + ey))
                        if b is not None and abs(b - s) >= 2 and self.A[b] == "H":
                            gain += 1
                opts.append((-gain, ch, p))
            for _, ch, p in sorted(opts):
                occ.add(p); pos[s] = p; inv[p] = s; out.append(ch)
                if rec(m + 1, out):
                    return True
                out.pop(); occ.discard(p); del pos[s]; del inv[p]
            return False

        out = []
        ok = rec(len(pos), out)
        return word + "".join(out) if ok else None

    # --- node data
    def node(self, word: str) -> Node:
        nd = self._cache.get(word)
        if nd is not None:
            return nd
        n, q, i = self.n, self.q, self.i
        pos = self.place(word)
        assert pos is not None
        m = len(word) + 1
        nd = Node(); nd.word = word; nd.wit = None; nd.dead = False; nd.complete = False
        if m == n:
            kA, kB = self.contacts(pos, self.A), self.contacts(pos, self.B)
            a, b = q ** kA, q ** kB
            nd.L = nd.U = a; nd.LB = nd.UB = b; nd.r = nd.R = b / a; nd.complete = True
        elif not self.child_letters(pos, m):
            nd.L = nd.U = nd.LB = nd.UB = nd.r = nd.R = Fraction(0); nd.dead = True
        else:
            r3 = Fraction(3) ** (n - m)
            nd.U = r3 * q ** (self.contacts(pos, self.A) + self.ucap(self.A, pos, m))
            nd.UB = r3 * q ** (self.contacts(pos, self.B) + self.ucap(self.B, pos, m))
            nd.L = nd.LB = Fraction(0)
            w = self.witness(word, dict(pos))
            if w is not None:
                wp = self.place(w)
                nd.wit = w
                nd.L = q ** self.contacts(wp, self.A)
                nd.LB = q ** self.contacts(wp, self.B)
            cmin, cmax = self.c_range(pos, m)
            d1, d2 = self.sgn * cmin, self.sgn * cmax
            nd.r, nd.R = q ** min(d1, d2), q ** max(d1, d2)
        self._cache[word] = nd
        return nd

    def c_range(self, pos, m):
        n, i = self.n, self.i
        if i in pos:
            inv = {p: k for k, p in pos.items()}
            x, y = pos[i]
            cmin = sum(1 for dx, dy in STEP.values()
                       if (b := inv.get((x + dx, y + dy))) is not None and abs(b - i) >= 2 and self.A[b] == "H")
            free = sum(1 for dx, dy in STEP.values() if (x + dx, y + dy) not in inv)
            free -= sum(1 for c in (i - 1, i + 1) if 0 <= c < n and c not in pos)
            nf = sum(1 for k in self.hpart if k not in pos and abs(k - i) >= 3 and (k - i) % 2)
            return cmin, cmin + max(0, min(nf, free))
        nf = sum(1 for k in self.hpart if abs(k - i) >= 3 and (k - i) % 2)
        return 0, min(nf, 3 if i in (0, n - 1) else 2)

    def children(self, word: str) -> list[str]:
        pos = self.place(word)
        m = len(word) + 1
        if m == self.n:
            return []
        return [word + ch for ch in self.child_letters(pos, m)]

    # --- frontiers
    def root(self) -> list[str]:
        return ["R"]

    def bfs_frontier(self, m: int) -> list[str]:
        """All valid prefixes with m placed sites (shorter if dead; complete if m >= n)."""
        out, stack = [], ["R"]
        while stack:
            w = stack.pop()
            if len(w) + 1 >= m:
                out.append(w)
                continue
            kids = self.children(w)
            if not kids:
                out.append(w)
            stack.extend(kids)
        return out

    def score(self, nd: Node, kind: str = "span") -> float:
        """Width contribution proxies for splitting (negative = never split).
        span : U*R - L*r, the width of the region's B-mass range (neutral; suits separate bounds)
        ratio: U*(R - r), the mass-weighted ratio uncertainty (what T1 cannot remove)"""
        if nd.complete or nd.dead:
            return -1.0
        if kind == "span":
            return float(nd.U * nd.R - nd.L * nd.r)
        return float(nd.U * (nd.R - nd.r)) if nd.R > nd.r else 0.0

    def refine(self, start: list[str] | None = None, kind: str = "span"):
        """Generator: yields the current frontier (list of words) after each split of the
        highest-scoring region. First yield is the start frontier."""
        words = set(start or self.root())
        heap = [(-self.score(self.node(w), kind), w) for w in words]
        heapq.heapify(heap)
        yield words
        while heap:
            s, w = heapq.heappop(heap)
            if -s < 0:
                return
            if -s == 0:      # nothing left that can narrow the ratio: finish remaining with span
                return
            words.discard(w)
            for c in self.children(w):
                words.add(c)
                heapq.heappush(heap, (-self.score(self.node(c), kind), c))
            yield words

    # --- bounds on ln(Z_B/Z_A) (kT units: ddF = -ln ratio, so width is the same)
    def trivial(self):
        """Global ||V||_inf-type bound: one region, ratio in [q^dmin, q^dmax] from the root."""
        nd = self.node("R") if self.n > 2 else None
        return nd.r, nd.R

    def coupled(self, words):
        regs = [self.node(w).region() for w in words]
        return lower_enclosure(regs).value, upper_enclosure(regs).value

    def separate(self, words):
        """Independent bounds on Z_A and Z_B from the same frontier (own caps for B)."""
        nodes = [self.node(w) for w in words]
        LA, UA = sum(x.L for x in nodes), sum(x.U for x in nodes)
        LB, UB = sum(x.LB for x in nodes), sum(x.UB for x in nodes)
        return (LB / UA, UB / LA if LA > 0 else None)

    def certificate(self, words, claim=None) -> dict:
        lo, hi = claim or self.coupled(words)
        fr = []
        for w in sorted(words):
            nd = self.node(w)
            e = {"w": w}
            if nd.wit is not None:
                e["wit"] = nd.wit
            fr.append(e)
        return {"seqA": self.A, "seqB": self.B, "q": str(self.q), "i": self.i,
                "order": self.order, "frontier": fr,
                "claim": {"lo": str(lo), "hi": str(hi)}}


def cert_json(cert: dict) -> str:
    return json.dumps(cert, separators=(",", ":"))
