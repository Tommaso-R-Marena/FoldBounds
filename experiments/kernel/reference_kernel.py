#!/usr/bin/env python3
"""Independent pure-Python reference implementation of ``FOLDBOUNDS-KERNEL 1``.

Written from the text of ``docs/research/CERTIFIED_FRONTIER_KERNEL_SPEC.md``
only (sections 2-6).  It does not import ``foldbounds``, the oracle or the
Lean checker, and it was not derived from the Lean sources.  Standard library
only.

Public API
----------
``parse(data: bytes) -> (seq, queries, tree_start)``
    Strict parser of section 2.  ``queries`` is a list of ``(Fraction, int)``;
    ``tree_start`` is the byte offset of ``<T>`` inside ``data`` (the tree is
    ``data[tree_start : len(data) - len(b"\\nend\\n")]``).  Raises
    ``Reject(category)`` on any deviation.
``check_tree(seq, tree_bytes) -> Acc``
    Sections 3-4: exact histogram (``2n + 2`` slots), recomputed bound entries
    ``(a, K, cnt, fam)`` and the node counters.  Raises ``Reject``.
``evaluate(acc, q, t) -> dict``
    Section 5, exact ``Fraction`` endpoints ``zLo .. pHi``.
``verify(data: bytes) -> dict``
    Composition, with the field names of the runner JSON (section 6).
``bound_entry(seq, sites) -> (a, K, cnt, fam)``
    Section 4 transcribed literally for one walk prefix (used by untrusted
    producers to score leaves and by the tests as a second implementation of
    the incremental computation inside ``check_tree``).

CLI: ``python reference_kernel.py check FILE...`` prints one JSON object per
file, like ``foldkernel check``.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from functools import lru_cache
from math import gcd
from typing import Iterable, Optional, Sequence

FORMAT = "FOLDBOUNDS-KERNEL 1"
_MAGIC = b"FOLDBOUNDS-KERNEL 1"
_TRAILER = b"\nend\n"

# Directions 0 = +x, 1 = -x, 2 = +y, 3 = -y.
_DX = (1, -1, 0, 0)
_DY = (0, 0, 1, -1)

_BYTE_N = ord("N")
_BYTE_E = ord("E")
_BYTE_X = ord("X")
_BYTE_B = ord("B")
_BYTE_BLOCKED = ord("-")

ENDPOINTS = ("zLo", "zHi", "eLo", "eHi", "nLo", "nHi", "pLo", "pHi")


class Reject(Exception):
    """The input is not an acceptable certificate; ``category`` says why."""

    def __init__(self, category: str, detail: str = ""):
        super().__init__(category if not detail else f"{category}: {detail}")
        self.category = category
        self.detail = detail


# ---------------------------------------------------------------------------
# Big-number helpers (CPython limits int <-> str conversions to 4300 digits)
# ---------------------------------------------------------------------------

def nat_from_digits(digits: bytes) -> int:
    """Value of a nonempty ASCII digit string of any length."""
    if len(digits) <= 2000:
        return int(digits)
    half = len(digits) // 2
    return nat_from_digits(digits[:half]) * 10 ** (len(digits) - half) + nat_from_digits(digits[half:])


def dec(n: int) -> str:
    """Decimal representation of an integer of any size."""
    if n < 0:
        return "-" + dec(-n)
    if n.bit_length() <= 8000:
        return str(n)
    half = (n.bit_length() * 30103 // 100000) // 2
    hi, lo = divmod(n, 10 ** half)
    return dec(hi) + dec(lo).rjust(half, "0")


def ratio_str(x: Fraction) -> str:
    """``num/den`` in lowest terms (``den`` may be 1)."""
    return f"{dec(x.numerator)}/{dec(x.denominator)}"


# ---------------------------------------------------------------------------
# Section 2: file format
# ---------------------------------------------------------------------------

def _canonical_nat(tok: bytes, category: str) -> int:
    """Canonical natural: ``0`` or a nonempty digit string without leading ``0``."""
    if len(tok) == 0:
        raise Reject(category, "empty number")
    for b in tok:
        if b < 0x30 or b > 0x39:
            raise Reject(category, "non-digit byte")
    if len(tok) > 1 and tok[0] == 0x30:
        raise Reject(category, "leading zero")
    return nat_from_digits(tok)


def _canonical_rat(tok: bytes) -> tuple[int, int]:
    """Canonical nonnegative rational: canonical natural, or ``num/den`` with
    canonical naturals, ``den >= 2`` and ``gcd(num, den) = 1``."""
    if b"/" not in tok:
        return _canonical_nat(tok, "bad-q"), 1
    parts = tok.split(b"/")
    if len(parts) != 2:
        raise Reject("bad-q", "more than one '/'")
    num = _canonical_nat(parts[0], "bad-q")
    den = _canonical_nat(parts[1], "bad-q")
    if den < 2:
        raise Reject("bad-q", "denominator must be >= 2")
    if gcd(num, den) != 1:
        raise Reject("bad-q", "not in lowest terms")
    return num, den


def parse(data: bytes) -> tuple[str, list[tuple[Fraction, int]], int]:
    """Strict parser of section 2.  Returns ``(seq, queries, tree_start)``."""
    if not isinstance(data, (bytes, bytearray, memoryview)):
        raise TypeError("parse() expects bytes")
    data = bytes(data)
    if len(data) == 0:
        raise Reject("empty-file")
    if not data.isascii():
        raise Reject("non-ascii-byte")
    if b"\r" in data:
        raise Reject("carriage-return")
    if data[-1:] != b"\n":
        raise Reject("missing-final-newline")
    lines = data[:-1].split(b"\n")

    if lines[0] != _MAGIC:
        raise Reject("bad-magic")

    if len(lines) < 2 or lines[1][:9] != b"sequence ":
        raise Reject("bad-sequence-line")
    seq = lines[1][9:]
    if len(seq) < 2:
        raise Reject("bad-sequence", "length must be >= 2")
    for b in seq:
        if b != 0x48 and b != 0x50:
            raise Reject("bad-sequence", "alphabet is {H, P}")

    queries: list[tuple[Fraction, int]] = []
    i = 2
    while i < len(lines) and lines[i][:6] == b"query ":
        fields = lines[i].split(b" ")
        if len(fields) != 3:
            raise Reject("bad-query-line", "expected 'query <q> <t>'")
        num, den = _canonical_rat(fields[1])
        t = _canonical_nat(fields[2], "bad-threshold")
        if num < den:
            raise Reject("q-below-one")
        queries.append((Fraction(num, den), t))
        i += 1
    if not queries:
        raise Reject("missing-query")

    if i >= len(lines) or lines[i][:5] != b"tree ":
        raise Reject("bad-tree-line")
    tree = lines[i][5:]
    tree_start = sum(len(line) + 1 for line in lines[:i]) + 5
    if i + 1 >= len(lines) or lines[i + 1] != b"end":
        raise Reject("bad-end-line")
    if i + 2 != len(lines):
        raise Reject("trailing-data")
    for b in tree:
        if b not in (_BYTE_N, _BYTE_E, _BYTE_X, _BYTE_B, _BYTE_BLOCKED):
            raise Reject("tree-bad-byte")
    if data[tree_start:len(data) - len(_TRAILER)] != tree:
        raise AssertionError("internal error: tree offset")
    return seq.decode("ascii"), queries, tree_start


# ---------------------------------------------------------------------------
# Section 4: bound entries
# ---------------------------------------------------------------------------

def cand(seq: str, i: int, j: int) -> bool:
    return i + 3 <= j and (j - i) % 2 == 1 and seq[i] == "H" and seq[j] == "H"


def slots(n: int, j: int) -> int:
    return 2 if j + 1 < n else 3


@lru_cache(maxsize=4096)
def _future_cap_table_literal(seq: str) -> tuple[int, ...]:
    """``futureCap(m)`` for ``m = 0 .. n`` straight from the definition."""
    n = len(seq)
    per_j = [min(slots(n, j), sum(1 for i in range(n) if cand(seq, i, j))) for j in range(n)]
    return tuple(sum(per_j[m:]) for m in range(n + 1))


def future_cap(seq: str, m: int) -> int:
    return _future_cap_table_literal(seq)[m]


def _future_cap_table_fast(seq: str) -> list[int]:
    """Same table in O(n) using per-parity prefix counts of H residues."""
    n = len(seq)
    even = [0] * (n + 1)
    odd = [0] * (n + 1)
    for i, ch in enumerate(seq):
        h = 1 if ch == "H" else 0
        even[i + 1] = even[i] + (h if i % 2 == 0 else 0)
        odd[i + 1] = odd[i] + (h if i % 2 == 1 else 0)
    table = [0] * (n + 1)
    for j in range(n - 1, -1, -1):
        c = 0
        if seq[j] == "H" and j >= 3:
            # partners i <= j - 3 of the opposite parity, i.e. indices < j - 2
            c = (odd if j % 2 == 0 else even)[j - 2]
        table[j] = table[j + 1] + min(slots(n, j), c)
    return table


def min_perim(total: int, w: int, h: int) -> int:
    best = 2 * total + 2
    for c in range(w, total + 1):
        v = c + max(h, (total + c - 1) // c)
        if v < best:
            best = v
    return best


def geom_cap(n: int, w: int, h: int) -> int:
    return max(0, (n + 1) - min_perim(n, w, h))


_PD: list[tuple[int, int, int]] = [(1, 1, 1)]


def pd_count(r: int) -> int:
    """``(f,u,d)_0 = (1,1,1)``, ``(f,u,d)_{r+1} = (f+u+d, f+u, f+d)``; returns ``f_r``."""
    while len(_PD) <= r:
        f, u, d = _PD[-1]
        _PD.append((f + u + d, f + u, f + d))
    return _PD[r][0]


def bound_entry(seq: str, sites: Sequence[tuple[int, int]]) -> tuple[int, int, int, int]:
    """Bound entry ``(a, K, cnt, fam)`` of the walk prefix ``sites`` (section 4).

    Literal transcription; ``sites`` are the placed sites in chain order.
    """
    n = len(seq)
    m = len(sites)
    if not 2 <= m <= n:
        raise ValueError("need 2 <= len(sites) <= len(seq)")
    r = n - m

    a = 0
    for i in range(m):
        if seq[i] != "H":
            continue
        for j in range(i + 2, m):
            if seq[j] == "H" and abs(sites[i][0] - sites[j][0]) + abs(sites[i][1] - sites[j][1]) == 1:
                a += 1

    xs = [s[0] for s in sites]
    ys = [s[1] for s in sites]
    w = max(xs) - min(xs) + 1
    h = max(ys) - min(ys) + 1
    cap = min(a + future_cap(seq, m), geom_cap(n, w, h))

    cnt = 3 ** r

    ex, ey = sites[-1]
    others = sites[:-1]
    strict_axis = (
        all(x < ex for x, _ in others)
        or all(x > ex for x, _ in others)
        or all(y < ey for _, y in others)
        or all(y > ey for _, y in others)
    )
    if strict_axis:
        fam = pd_count(r)
    elif any(
        all(sx * ex + sy * ey >= sx * x + sy * y for x, y in sites)
        for sx in (1, -1) for sy in (1, -1)
    ):
        fam = 2 ** r
    else:
        fam = 0
    return a, cap, cnt, fam


# ---------------------------------------------------------------------------
# Sections 3-4: tree check
# ---------------------------------------------------------------------------

class Acc:
    """Result of an accepted tree: exact histogram, bound entries, counters."""

    __slots__ = ("n", "hist", "entries", "split_nodes", "bound_leaves",
                 "enum_leaves", "exact_leaves")

    def __init__(self, n: int):
        self.n = n
        self.hist = [0] * (2 * n + 2)
        self.entries: list[tuple[int, int, int, int]] = []
        self.split_nodes = 0
        self.bound_leaves = 0
        self.enum_leaves = 0
        self.exact_leaves = 0

    @property
    def exact_count(self) -> int:
        return sum(self.hist)


def _enumerate_completions(is_h: list[bool], n: int, occ: dict, x: int, y: int,
                           m: int, k0: int, hist: list[int]) -> None:
    """Add the contact count of every completion of the current walk to ``hist``.

    ``occ`` maps placed sites to residue indices (``m`` sites placed, the last
    one at ``(x, y)``, ``k0`` contacts among them).  ``occ`` is restored.
    """
    slots_ = len(hist)
    if m == n:
        if k0 >= slots_:
            raise Reject("histogram-overflow")
        hist[k0] += 1
        return
    last = n - 1
    px = [x]
    py = [y]
    kk = [k0]
    nd = [0]
    depth = 0
    while depth >= 0:
        d = nd[depth]
        if d == 4:
            if depth > 0:
                del occ[(px[depth], py[depth])]
            px.pop()
            py.pop()
            kk.pop()
            nd.pop()
            depth -= 1
            continue
        nd[depth] = d + 1
        nx = px[depth] + _DX[d]
        ny = py[depth] + _DY[d]
        if (nx, ny) in occ:
            continue
        j = m + depth
        k = kk[depth]
        if is_h[j]:
            for e in range(4):
                i = occ.get((nx + _DX[e], ny + _DY[e]))
                if i is not None and i + 2 <= j and is_h[i]:
                    k += 1
        if j == last:
            if k >= slots_:
                raise Reject("histogram-overflow")
            hist[k] += 1
        else:
            occ[(nx, ny)] = j
            px.append(nx)
            py.append(ny)
            kk.append(k)
            nd.append(0)
            depth += 1


def check_tree(seq: str, tree_bytes: bytes, paranoid: bool = False) -> Acc:
    """Check a certificate tree against ``seq`` (sections 3 and 4).

    Iterative (no recursion), O(1) work per tree byte apart from ``X`` leaves.
    With ``paranoid=True`` every bound entry is recomputed with the literal
    ``bound_entry`` and compared (slow; used by the tests).
    """
    if isinstance(tree_bytes, str):
        tree_bytes = tree_bytes.encode("latin-1")
    tree = bytes(tree_bytes)
    n = len(seq)
    if n < 2 or any(ch not in "HP" for ch in seq):
        raise Reject("bad-sequence")
    is_h = [ch == "H" for ch in seq]
    acc = Acc(n)
    hist = acc.hist
    entries = acc.entries
    fcap = _future_cap_table_fast(seq)
    geom_memo: dict[tuple[int, int], int] = {}

    # Walk state.  Index m-1 of each list describes the prefix with m sites.
    xs = [0, 1]
    ys = [0, 0]
    occ: dict[tuple[int, int], int] = {(0, 0): 0, (1, 0): 1}
    contacts = [0, 0]
    box = [(0, 0, 0, 0), (0, 1, 0, 0)]              # (min x, max x, min y, max y)
    diag = [(0, 0, 0, 0), (1, 1, 0, 0)]             # max of x+y, x-y, -x+y, -x-y

    def push(x: int, y: int) -> None:
        j = len(xs)
        c = contacts[-1]
        if is_h[j]:
            for e in range(4):
                i = occ.get((x + _DX[e], y + _DY[e]))
                if i is not None and i + 2 <= j and is_h[i]:
                    c += 1
        x0, x1, y0, y1 = box[-1]
        d0, d1, d2, d3 = diag[-1]
        occ[(x, y)] = j
        xs.append(x)
        ys.append(y)
        contacts.append(c)
        box.append((min(x0, x), max(x1, x), min(y0, y), max(y1, y)))
        diag.append((max(d0, x + y), max(d1, x - y), max(d2, -x + y), max(d3, -x - y)))

    def pop() -> None:
        del occ[(xs.pop(), ys.pop())]
        contacts.pop()
        box.pop()
        diag.pop()

    def bound_leaf() -> None:
        m = len(xs)
        r = n - m
        a = contacts[-1]
        x0, x1, y0, y1 = box[-1]
        wh = (x1 - x0 + 1, y1 - y0 + 1)
        g = geom_memo.get(wh)
        if g is None:
            g = geom_cap(n, wh[0], wh[1])
            geom_memo[wh] = g
        cap = min(a + fcap[m], g)
        ex = xs[-1]
        ey = ys[-1]
        px0, px1, py0, py1 = box[-2]                # all other placed sites
        if ex > px1 or ex < px0 or ey > py1 or ey < py0:
            fam = pd_count(r)
        else:
            d0, d1, d2, d3 = diag[-2]
            if ex + ey >= d0 or ex - ey >= d1 or -ex + ey >= d2 or -ex - ey >= d3:
                fam = 2 ** r
            else:
                fam = 0
        entry = (a, cap, 3 ** r, fam)
        if paranoid:
            literal = bound_entry(seq, list(zip(xs, ys)))
            if literal != entry:
                raise AssertionError(
                    f"incremental bound entry {entry} != literal {literal} "
                    f"for {seq} at {list(zip(xs, ys))}")
        entries.append(entry)

    size = len(tree)
    pos = 0
    open_nodes: list[int] = []          # next direction of every open N node
    while True:
        # Read one node encoding for the current walk.
        if pos >= size:
            raise Reject("tree-truncated")
        b = tree[pos]
        pos += 1
        r = n - len(xs)
        if b == _BYTE_N:
            if r < 1:
                raise Reject("tree-split-at-terminal")
            acc.split_nodes += 1
            open_nodes.append(0)
        else:
            if b == _BYTE_E:
                if r != 0:
                    raise Reject("tree-exact-leaf-not-terminal")
                k = contacts[-1]
                if k >= len(hist):
                    raise Reject("histogram-overflow")
                hist[k] += 1
                acc.exact_leaves += 1
            elif b == _BYTE_X:
                _enumerate_completions(is_h, n, occ, xs[-1], ys[-1], len(xs), contacts[-1], hist)
                acc.enum_leaves += 1
            elif b == _BYTE_B:
                bound_leaf()
                acc.bound_leaves += 1
            elif b == _BYTE_BLOCKED:
                # '-' where a node is required: at the root, or for a free site.
                raise Reject("tree-blocked-marker-on-free-site")
            else:
                raise Reject("tree-bad-byte")
            if not open_nodes:
                break                   # the root itself was a leaf
            pop()
        # Move on to the next child that needs a node encoding.
        descended = False
        while open_nodes:
            d = open_nodes[-1]
            if d == 4:
                open_nodes.pop()
                if not open_nodes:
                    break               # the root split is complete
                pop()
                continue
            open_nodes[-1] = d + 1
            nx = xs[-1] + _DX[d]
            ny = ys[-1] + _DY[d]
            if (nx, ny) in occ:
                if pos >= size:
                    raise Reject("tree-truncated")
                if tree[pos] != _BYTE_BLOCKED:
                    if tree[pos] not in (_BYTE_N, _BYTE_E, _BYTE_X, _BYTE_B):
                        raise Reject("tree-bad-byte")
                    raise Reject("tree-node-on-occupied-site")
                pos += 1
            else:
                push(nx, ny)
                descended = True
                break
        if not descended:
            break
    if pos != size:
        raise Reject("tree-trailing-bytes")
    return acc


# ---------------------------------------------------------------------------
# Section 5: evaluation
# ---------------------------------------------------------------------------

def evaluate(acc: Acc, q, t: int) -> dict:
    """Exact endpoints for the query ``(q, t)`` (section 5)."""
    q = Fraction(q)
    if not isinstance(t, int) or isinstance(t, bool) or t < 0:
        raise ValueError("t must be a natural number")
    # Entries with the same (a, K) contribute sums of cnt and fam.
    groups: dict[tuple[int, int], list[int]] = {}
    for a, cap, cnt, fam in acc.entries:
        g = groups.get((a, cap))
        if g is None:
            groups[(a, cap)] = [cnt, fam]
        else:
            g[0] += cnt
            g[1] += fam

    powers: dict[int, Fraction] = {}

    def qpow(k: int) -> Fraction:
        v = powers.get(k)
        if v is None:
            v = q ** k
            powers[k] = v
        return v

    ex_all = Fraction(0)
    ex_event = Fraction(0)
    ex_rest = Fraction(0)
    for k, count in enumerate(acc.hist):
        if count:
            term = count * qpow(k)
            ex_all += term
            if k >= t:
                ex_event += term
            else:
                ex_rest += term

    z_lo = ex_all
    z_hi = ex_all
    e_lo = ex_event
    e_hi = ex_event
    n_lo = ex_rest
    n_hi = ex_rest
    for (a, cap), (cnt, fam) in groups.items():
        low = fam * qpow(a)
        high = cnt * qpow(cap)
        z_lo += low
        z_hi += high
        if t <= a:
            e_lo += low
        if t <= cap:
            e_hi += high
        if cap < t:
            n_lo += low
        if a < t:
            n_hi += cnt * qpow(min(cap, t - 1))
    p_lo = Fraction(0) if e_lo == 0 else e_lo / (e_lo + n_hi)
    if e_hi == 0:
        p_hi = Fraction(0)
    elif n_lo == 0:
        p_hi = Fraction(1)
    else:
        p_hi = e_hi / (e_hi + n_lo)
    return {"zLo": z_lo, "zHi": z_hi, "eLo": e_lo, "eHi": e_hi,
            "nLo": n_lo, "nHi": n_hi, "pLo": p_lo, "pHi": p_hi}


def evaluate_per_entry(acc: Acc, q, t: int) -> dict:
    """Section 5 transcribed entry by entry (slow; cross-check of ``evaluate``)."""
    q = Fraction(q)
    ex = acc.hist
    z_lo = sum((ex[k] * q ** k for k in range(len(ex))), Fraction(0)) \
        + sum((fam * q ** a for a, _, _, fam in acc.entries), Fraction(0))
    z_hi = sum((ex[k] * q ** k for k in range(len(ex))), Fraction(0)) \
        + sum((cnt * q ** cap for _, cap, cnt, _ in acc.entries), Fraction(0))
    e_lo = sum((ex[k] * q ** k for k in range(len(ex)) if k >= t), Fraction(0)) \
        + sum((fam * q ** a for a, _, _, fam in acc.entries if t <= a), Fraction(0))
    e_hi = sum((ex[k] * q ** k for k in range(len(ex)) if k >= t), Fraction(0)) \
        + sum((cnt * q ** cap for _, cap, cnt, _ in acc.entries if t <= cap), Fraction(0))
    n_lo = sum((ex[k] * q ** k for k in range(len(ex)) if k < t), Fraction(0)) \
        + sum((fam * q ** a for a, cap, _, fam in acc.entries if cap < t), Fraction(0))
    n_hi = sum((ex[k] * q ** k for k in range(len(ex)) if k < t), Fraction(0)) \
        + sum((cnt * q ** min(cap, t - 1) for a, cap, cnt, _ in acc.entries if a < t), Fraction(0))
    p_lo = Fraction(0) if e_lo == 0 else e_lo / (e_lo + n_hi)
    p_hi = Fraction(0) if e_hi == 0 else (Fraction(1) if n_lo == 0 else e_hi / (e_hi + n_lo))
    return {"zLo": z_lo, "zHi": z_hi, "eLo": e_lo, "eHi": e_hi,
            "nLo": n_lo, "nHi": n_hi, "pLo": p_lo, "pHi": p_hi}


# ---------------------------------------------------------------------------
# Section 6: composition
# ---------------------------------------------------------------------------

def verify(data: bytes) -> dict:
    """Full check of a certificate file; mirrors the runner JSON of section 6.

    Accepted: ``accepted, format, sequence, n, cert_bytes, split_nodes,
    bound_leaves, enum_leaves, exact_leaves, exact_count, exact_histogram,
    results``; every result has ``q`` (``num/den``), ``threshold`` (int) and
    the eight endpoints as ``num/den`` strings in lowest terms.
    ``exact_histogram`` is printed without trailing zeros, like the runner.
    Rejected: ``{"accepted": False, "reason": <category>}``.
    """
    try:
        seq, queries, tree_start = parse(data)
        tree = bytes(data)[tree_start:len(data) - len(_TRAILER)]
        acc = check_tree(seq, tree)
        results = []
        for q, t in queries:
            values = evaluate(acc, q, t)
            row = {"q": ratio_str(q), "threshold": t}
            for name in ENDPOINTS:
                row[name] = ratio_str(values[name])
            results.append(row)
    except Reject as rej:
        return {"accepted": False, "reason": rej.category}
    histogram = list(acc.hist)
    while histogram and histogram[-1] == 0:
        histogram.pop()
    return {
        "accepted": True,
        "format": FORMAT,
        "sequence": seq,
        "n": len(seq),
        "cert_bytes": len(tree),
        "split_nodes": acc.split_nodes,
        "bound_leaves": acc.bound_leaves,
        "enum_leaves": acc.enum_leaves,
        "exact_leaves": acc.exact_leaves,
        "exact_count": acc.exact_count,
        "exact_histogram": histogram,
        "results": results,
    }


def to_json(obj) -> str:
    """Minimal JSON writer that can print integers of any size."""
    if obj is True:
        return "true"
    if obj is False:
        return "false"
    if obj is None:
        return "null"
    if isinstance(obj, int):
        return dec(obj)
    if isinstance(obj, str):
        out = ['"']
        for ch in obj:
            if ch == '"' or ch == "\\":
                out.append("\\" + ch)
            elif ord(ch) < 0x20 or ord(ch) > 0x7E:
                out.append(f"\\u{ord(ch):04x}" if ord(ch) <= 0xFFFF else ch)
            else:
                out.append(ch)
        out.append('"')
        return "".join(out)
    if isinstance(obj, (list, tuple)):
        return "[" + ",".join(to_json(v) for v in obj) + "]"
    if isinstance(obj, dict):
        return "{" + ",".join(f"{to_json(str(k))}:{to_json(v)}" for k, v in obj.items()) + "}"
    raise TypeError(f"cannot serialise {type(obj).__name__}")


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) < 2 or args[0] != "check":
        print("usage: reference_kernel.py check CERT [CERT ...]", file=sys.stderr)
        return 2
    all_ok = True
    for name in args[1:]:
        try:
            with open(name, "rb") as fh:
                data = fh.read()
        except OSError as exc:
            print(f"cannot read {name}: {exc}", file=sys.stderr)
            return 2
        report = verify(data)
        all_ok = all_ok and report["accepted"]
        print(to_json({"file": name, **report}))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
