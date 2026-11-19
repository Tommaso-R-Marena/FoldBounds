#!/usr/bin/env python3
"""Matched certificate generators for foldkernel, foldcheck (v1) and foldcheck2 (v2).

UNTRUSTED producer code of the benchmark harness (standard library only; no
checker and no ``foldbounds`` import).  One depth-first walk over the legal
prefixes of a sequence yields a *frontier* (a prefix-free cover of the
conformation tree).  The same frontier object is then rendered for the three
checkers, so that every comparison is made on identical covers:

* kernel   ``FOLDBOUNDS-KERNEL 1``   preorder byte tree over ``N E X B -``
* v1       ``FOLDBOUNDS-FRONTIER 1`` flat lines ``E/D/B/T <word>``
* v2       ``FOLDBOUNDS-FRONTIER 2`` flat lines ``E/D/X/B/T/M<k> <word>``

Families
--------
``uniform(d)``        every legal prefix with ``d`` placed residues is a leaf
                      (complete conformations are exact leaves, trapped
                      shorter prefixes are dead ends); all depth-``d`` leaves
                      get the bounded leaf kind.
``hybrid(d, kappa)``  the uniform(d) frontier; a depth-``d`` leaf is the
                      exhaustive kind ``X`` when its prefix already has
                      ``>= kappa`` H-H contacts, else the bounded kind.
                      ``kappa = None`` means "never" (= uniform(d));
                      ``kappa = 0`` makes every leaf exhaustive.
``full``              the single root leaf ``X``.

Directions: 0 = +x, 1 = -x, 2 = +y, 3 = -y (the fixed first bond (0,0)->(1,0)
is not part of a word).  Depth ``d`` counts placed residues, as in
``foldbounds.lattice.build_frontier``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Optional

KERNEL_MAGIC = "FOLDBOUNDS-KERNEL 1"
V1_MAGIC = "FOLDBOUNDS-FRONTIER 1"
V2_MAGIC = "FOLDBOUNDS-FRONTIER 2"

_N, _E, _X, _B, _DASH = b"N"[0], b"E"[0], b"X"[0], b"B"[0], b"-"[0]
_LEAF0 = 128            # depth-d leaf with k known contacts is stored as byte 128 + k


def check_sequence(seq: str) -> int:
    if not isinstance(seq, str) or len(seq) < 2 or any(ch not in "HP" for ch in seq):
        raise ValueError("sequence must be an H/P string of length >= 2")
    return len(seq)


def q_text(q: Fraction) -> str:
    """Canonical rational accepted by all three parsers (``2`` or ``3/2``)."""
    q = Fraction(q)
    if q < 1:
        raise ValueError("q must be >= 1")
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


@dataclass
class Frontier:
    """Uniform-depth prefix frontier of one sequence (see module docstring)."""
    sequence: str
    depth: int                       # placed residues per leaf (clamped to n)
    base: bytes                      # kernel preorder tree, depth-d leaves as 128 + contacts
    entries: Optional[list]          # [(word, code, contacts)], code in "L" "E" "D"; DFS order
    n_depth_leaves: int = 0          # leaves with exactly ``depth`` < n placed residues
    n_exact: int = 0                 # complete conformations listed as exact leaves
    n_dead: int = 0                  # trapped prefixes above the frontier depth
    n_split: int = 0                 # kernel split nodes (dead ends included)
    contact_hist: dict = field(default_factory=dict)   # known contacts -> number of depth-d leaves

    @property
    def n_leaves(self) -> int:
        return self.n_depth_leaves + self.n_exact

    def count_exhaustive(self, kappa: Optional[int]) -> int:
        if kappa is None:
            return 0
        return sum(c for k, c in self.contact_hist.items() if k >= kappa)


def build_frontier(seq: str, depth: int, want_words: bool = True) -> Frontier:
    """Depth-first walk producing the uniform(depth) frontier of ``seq``.

    ``want_words=False`` skips the flat word list (only the kernel tree is
    produced); use it for deep frontiers that v1/v2 are not asked to check.
    """
    n = check_sequence(seq)
    if not isinstance(depth, int) or depth < 2:
        raise ValueError("depth must be an integer >= 2")
    d = min(depth, n)
    is_h = [1 if ch == "H" else 0 for ch in seq]
    width = 2 * n + 3
    occ = bytearray(width * width)
    hh = bytearray(width * width)
    start = (n + 1) * width + (n + 1)
    occ[start] = 1
    hh[start] = is_h[0]
    occ[start + 1] = 1
    hh[start + 1] = is_h[1]
    steps = (1, -1, width, -width)                # directions 0, 1, 2, 3
    out = bytearray()
    entries: Optional[list] = [] if want_words else None
    word: list[str] = []
    stats = {"L": 0, "E": 0, "D": 0, "N": 0}
    chist: dict[int, int] = {}

    def leaf(m: int, k: int) -> None:
        if m == n:
            out.append(_E)
            stats["E"] += 1
            if entries is not None:
                entries.append(("".join(word) or ".", "E", k))
        else:
            out.append(_LEAF0 + k)
            stats["L"] += 1
            chist[k] = chist.get(k, 0) + 1
            if entries is not None:
                entries.append(("".join(word) or ".", "L", k))

    def split(pos: int, m: int, k: int) -> None:
        """Node with m < d placed residues ending at ``pos``, k known contacts."""
        out.append(_N)
        stats["N"] += 1
        h = is_h[m]
        child_is_leaf = (m + 1 == d)
        alive = False
        for di in range(4):
            p = pos + steps[di]
            if occ[p]:
                out.append(_DASH)
                continue
            alive = True
            if h:
                kk = k + hh[p + 1] + hh[p - 1] + hh[p + width] + hh[p - width] - hh[pos]
            else:
                kk = k
            if entries is not None:
                word.append("0123"[di])
            if child_is_leaf:
                leaf(m + 1, kk)
            else:
                occ[p] = 1
                hh[p] = h
                split(p, m + 1, kk)
                occ[p] = 0
                hh[p] = 0
            if entries is not None:
                word.pop()
        if not alive:
            stats["D"] += 1
            if entries is not None:
                entries.append(("".join(word) or ".", "D", k))

    if d == 2:
        leaf(2, 0)
    else:
        split(start + 1, 2, 0)
    return Frontier(sequence=seq, depth=d, base=bytes(out), entries=entries,
                    n_depth_leaves=stats["L"], n_exact=stats["E"], n_dead=stats["D"],
                    n_split=stats["N"], contact_hist=dict(sorted(chist.items())))


# ----------------------------------------------------------------------------
# Rendering
# ----------------------------------------------------------------------------

def _leaf_table(kappa: Optional[int], bounded: int, exhaustive: int) -> bytes:
    table = bytearray(range(256))
    for b in range(_LEAF0, 256):
        k = b - _LEAF0
        table[b] = exhaustive if (kappa is not None and k >= kappa) else bounded
    return bytes(table)


def kernel_tree(fr: Frontier, kappa: Optional[int] = None, bounded: str = "B") -> bytes:
    """Kernel preorder tree of the frontier under the hybrid(depth, kappa) rule."""
    if bounded not in ("B", "X"):
        raise ValueError("kernel bounded leaf kind must be B (or X)")
    return fr.base.translate(_leaf_table(kappa, ord(bounded), _X))


def kernel_text(seq: str, queries, tree: bytes) -> bytes:
    """Complete kernel certificate; ``queries`` is an iterable of (q, t)."""
    check_sequence(seq)
    lines = [KERNEL_MAGIC, f"sequence {seq}"]
    nq = 0
    for q, t in queries:
        if not isinstance(t, int) or t < 0:
            raise ValueError("threshold must be a natural number")
        lines.append(f"query {q_text(q)} {t}")
        nq += 1
    if nq == 0:
        raise ValueError("at least one query is required")
    head = ("\n".join(lines) + "\ntree ").encode("ascii")
    return head + bytes(tree) + b"\nend\n"


def flat_text(fr: Frontier, q: Fraction, t: int, *, version: int, bounded: str,
              kappa: Optional[int] = None) -> bytes:
    """Flat v1/v2 certificate on the same frontier.

    ``bounded`` is the kind of non-exhaustive depth-d leaves: ``B`` or ``T``
    (both versions) or ``M<k>`` (version 2 only).  Exhaustive leaves (``X``)
    exist only in version 2, so version 1 requires ``kappa is None``.
    """
    if fr.entries is None:
        raise ValueError("frontier was built without words")
    if version == 1:
        if bounded not in ("B", "T"):
            raise ValueError("version 1 has the bounded kinds B and T only")
        if kappa is not None:
            raise ValueError("version 1 has no exhaustive leaf kind")
        magic = V1_MAGIC
    elif version == 2:
        ok = bounded in ("B", "T") or (bounded[:1] == "M" and bounded[1:].isdigit())
        if not ok:
            raise ValueError("version 2 bounded kinds are B, T and M<k>")
        magic = V2_MAGIC
    else:
        raise ValueError("version must be 1 or 2")
    lines = [magic, f"sequence {fr.sequence}", f"q {q_text(q)}", f"threshold {t}"]
    for w, code, k in fr.entries:
        if code == "E":
            lines.append(f"E {w}")
        elif code == "D":
            lines.append(f"D {w}")
        elif kappa is not None and k >= kappa:
            lines.append(f"X {w}")
        else:
            lines.append(f"{bounded} {w}")
    lines.append("end")
    return ("\n".join(lines) + "\n").encode("ascii")


# ----------------------------------------------------------------------------
# Convenience constructors for the three families
# ----------------------------------------------------------------------------

def kernel_full(seq: str, queries) -> bytes:
    return kernel_text(seq, queries, b"X")


def kernel_root_bound(seq: str, queries) -> bytes:
    return kernel_text(seq, queries, b"B")


def v2_full(seq: str, q: Fraction, t: int) -> bytes:
    check_sequence(seq)
    return (f"{V2_MAGIC}\nsequence {seq}\nq {q_text(q)}\nthreshold {t}\nX .\nend\n").encode("ascii")


def kernel_uniform(seq: str, d: int, queries) -> bytes:
    return kernel_text(seq, queries, kernel_tree(build_frontier(seq, d, want_words=False)))


def kernel_hybrid(seq: str, d: int, kappa: Optional[int], queries) -> bytes:
    return kernel_text(seq, queries, kernel_tree(build_frontier(seq, d, want_words=False), kappa))


def describe_family(family: str, d: Optional[int] = None, kappa="-") -> str:
    if family == "full":
        return "full"
    if family == "uniform":
        return f"uniform({d})"
    kap = "inf" if kappa is None else kappa
    return f"hybrid({d},{kap})"


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("usage: certgen.py SEQUENCE DEPTH [KAPPA|inf] [Q] [T]   (prints the kernel certificate)",
              file=sys.stderr)
        sys.exit(2)
    _seq, _d = sys.argv[1], int(sys.argv[2])
    _kap = None if len(sys.argv) < 4 or sys.argv[3] == "inf" else int(sys.argv[3])
    _q = Fraction(sys.argv[4]) if len(sys.argv) > 4 else Fraction(2)
    _t = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    sys.stdout.buffer.write(kernel_hybrid(_seq, _d, _kap, [(_q, _t)]))
