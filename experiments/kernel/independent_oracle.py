#!/usr/bin/env python3
"""Independent exact oracle for the anchored 2D HP lattice model.

Standalone: standard library only, imports nothing from this repository (in
particular neither ``foldbounds`` nor the reference kernel).  It enumerates
every conformation and counts contacts straight from the model definition
(section 1 of ``docs/research/CERTIFIED_FRONTIER_KERNEL_SPEC.md``):

* sequence ``S`` over ``{H, P}``, length ``n >= 2``;
* conformations are self-avoiding walks of ``n`` sites on ``Z^2`` whose first
  two sites are ``(0,0)`` and ``(1,0)``; reflections are NOT identified;
* ``k(c)`` = number of pairs ``i < j``, ``j >= i + 2``, ``S[i] = S[j] = H``,
  sites ``i`` and ``j`` at Manhattan distance 1.

The oracle knows nothing about certificates, caps or the checker; it is the
ground truth of the differential tests.

Two enumerators are provided on purpose:

* ``density_of_states`` - fast depth-first search on a flat occupancy grid
  with incremental contact counting;
* ``density_of_states_naive`` - coordinate/set based enumeration scored with
  the literal ``O(n^2)`` double loop ``naive_contacts``.

``self_check`` compares them for every sequence up to a given length.

CLI::

    python independent_oracle.py HPHPPHHPHH --q 3/2 --t 2
    python independent_oracle.py --self-check 8
"""

from __future__ import annotations

import argparse
import atexit
import json
import os
import sys
import tempfile
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Callable, Iterator, Optional

# Number of n-step self-avoiding walks on Z^2 (OEIS A001411), c_0 .. c_15.
SAW_COUNTS = (
    1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292,
    324932, 881500, 2374444, 6416596,
)

CACHE_MODEL_TAG = "anchored-2d-hp/contacts/v1"
_DEFAULT_CACHE = Path(__file__).resolve().parent / "results" / "oracle_cache.json"

Site = tuple[int, int]


def validate_sequence(seq: str) -> int:
    if not isinstance(seq, str):
        raise TypeError("sequence must be a str")
    if len(seq) < 2:
        raise ValueError("sequence must have length >= 2")
    if any(ch not in "HP" for ch in seq):
        raise ValueError("sequence must be over the alphabet {H, P}")
    return len(seq)


def expected_conformations(n: int) -> Optional[int]:
    """c_{n-1} / 4 for the tabulated range, else ``None``."""
    if 2 <= n <= len(SAW_COUNTS):
        return SAW_COUNTS[n - 1] // 4
    return None


# ---------------------------------------------------------------------------
# Fast enumerator (flat grid, incremental contacts)
# ---------------------------------------------------------------------------

def density_of_states(seq: str) -> list[int]:
    """``g[k]`` = number of conformations of ``seq`` with exactly ``k`` contacts.

    Trailing zeros are stripped, so ``len(g) - 1`` is the maximal contact count.
    """
    n = validate_sequence(seq)
    is_h = [ch == "H" for ch in seq]
    side = 2 * n + 1                      # coordinates -n .. n in both axes
    grid = bytearray(side * side)         # 0 empty, 1 P residue, 2 H residue
    steps = (1, -1, side, -side)          # +x, -x, +y, -y
    hist = [0] * (n + 2)
    origin = n * side + n
    grid[origin] = 2 if is_h[0] else 1
    grid[origin + 1] = 2 if is_h[1] else 1
    last = n - 1

    def extend(pos: int, j: int, k: int) -> None:
        # Place residue j on a free neighbour of ``pos`` (site of residue j-1).
        h = is_h[j]
        mark = 2 if h else 1
        # Residue j-1 is always adjacent; it is a chain neighbour, not a contact.
        chain = 1 if is_h[j - 1] else 0
        for st in steps:
            p = pos + st
            if grid[p]:
                continue
            if h:
                kk = (k + (grid[p + 1] == 2) + (grid[p - 1] == 2)
                      + (grid[p + side] == 2) + (grid[p - side] == 2) - chain)
            else:
                kk = k
            if j == last:
                hist[kk] += 1
            else:
                grid[p] = mark
                extend(p, j + 1, kk)
                grid[p] = 0

    if n == 2:
        hist[0] = 1
    else:
        old = sys.getrecursionlimit()
        if old < n + 200:
            sys.setrecursionlimit(n + 200)
        try:
            extend(origin + 1, 2, 0)
        finally:
            sys.setrecursionlimit(old)
    while len(hist) > 1 and hist[-1] == 0:
        hist.pop()
    return hist


def prefix_stats(seq: str,
                 visit: Callable[[list[Site], int, Optional[int], Optional[int]], None]
                 ) -> tuple[int, Optional[int], Optional[int]]:
    """Visit every valid prefix of ``seq`` with the statistics of its completions.

    A prefix is a self-avoiding walk with ``m`` sites, ``2 <= m <= n``, starting
    ``(0,0),(1,0)``.  ``visit(sites, count, kmin, kmax)`` is called once per
    prefix (post-order); ``sites`` is the live list of placed sites (copy it if
    it must be kept), ``count`` the number of complete conformations extending
    the prefix and ``kmin``/``kmax`` the extreme contact counts among them
    (``None`` when ``count == 0``).  Returns the statistics of the root.
    """
    n = validate_sequence(seq)
    is_h = [ch == "H" for ch in seq]
    side = 2 * n + 1
    grid = bytearray(side * side)
    moves = ((1, 1, 0), (-1, -1, 0), (side, 0, 1), (-side, 0, -1))
    origin = n * side + n
    grid[origin] = 2 if is_h[0] else 1
    grid[origin + 1] = 2 if is_h[1] else 1
    sites: list[Site] = [(0, 0), (1, 0)]

    def rec(pos: int, k: int) -> tuple[int, Optional[int], Optional[int]]:
        j = len(sites)
        if j == n:
            visit(sites, 1, k, k)
            return 1, k, k
        h = is_h[j]
        mark = 2 if h else 1
        chain = 1 if is_h[j - 1] else 0
        x, y = sites[-1]
        total = 0
        lo: Optional[int] = None
        hi: Optional[int] = None
        for st, dx, dy in moves:
            p = pos + st
            if grid[p]:
                continue
            if h:
                kk = (k + (grid[p + 1] == 2) + (grid[p - 1] == 2)
                      + (grid[p + side] == 2) + (grid[p - side] == 2) - chain)
            else:
                kk = k
            grid[p] = mark
            sites.append((x + dx, y + dy))
            c, l, u = rec(p, kk)
            sites.pop()
            grid[p] = 0
            if c:
                total += c
                lo = l if lo is None or l < lo else lo
                hi = u if hi is None or u > hi else hi
        visit(sites, total, lo, hi)
        return total, lo, hi

    old = sys.getrecursionlimit()
    if old < n + 200:
        sys.setrecursionlimit(n + 200)
    try:
        return rec(origin + 1, 0)
    finally:
        sys.setrecursionlimit(old)


# ---------------------------------------------------------------------------
# Naive enumerator / scorer (literal definitions)
# ---------------------------------------------------------------------------

def enumerate_conformations(n: int) -> Iterator[tuple[Site, ...]]:
    """All anchored self-avoiding walks with ``n`` sites, as tuples of sites."""
    if n < 2:
        raise ValueError("n must be >= 2")
    walk: list[Site] = [(0, 0), (1, 0)]
    used = {(0, 0), (1, 0)}

    def rec() -> Iterator[tuple[Site, ...]]:
        if len(walk) == n:
            yield tuple(walk)
            return
        x, y = walk[-1]
        for nxt in ((x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)):
            if nxt in used:
                continue
            used.add(nxt)
            walk.append(nxt)
            yield from rec()
            walk.pop()
            used.remove(nxt)

    yield from rec()


def naive_contacts(seq: str, sites) -> int:
    """Literal contact count over the first ``len(sites)`` residues: O(m^2)."""
    k = 0
    m = len(sites)
    for i in range(m):
        if seq[i] != "H":
            continue
        xi, yi = sites[i]
        for j in range(i + 2, m):
            if seq[j] != "H":
                continue
            xj, yj = sites[j]
            if abs(xi - xj) + abs(yi - yj) == 1:
                k += 1
    return k


def density_of_states_naive(seq: str) -> list[int]:
    n = validate_sequence(seq)
    hist = [0] * (n + 2)
    for conf in enumerate_conformations(n):
        hist[naive_contacts(seq, conf)] += 1
    while len(hist) > 1 and hist[-1] == 0:
        hist.pop()
    return hist


def prefix_stats_naive(seq: str) -> dict[tuple[Site, ...], tuple[int, Optional[int], Optional[int]]]:
    """Completion statistics of every prefix, by brute-force grouping."""
    n = validate_sequence(seq)
    out: dict[tuple[Site, ...], tuple[int, Optional[int], Optional[int]]] = {}
    # Prefixes without any completion still have to be listed (count 0).
    for m in range(2, n + 1):
        for pre in enumerate_conformations(m):
            out[pre] = (0, None, None)
    for conf in enumerate_conformations(n):
        k = naive_contacts(seq, conf)
        for m in range(2, n + 1):
            key = conf[:m]
            c, lo, hi = out[key]
            out[key] = (c + 1,
                        k if lo is None or k < lo else lo,
                        k if hi is None or k > hi else hi)
    return out


# ---------------------------------------------------------------------------
# On-disk cache (keyed by sequence)
# ---------------------------------------------------------------------------

_cache: Optional[dict[str, list[int]]] = None
_cache_dirty = 0
_cache_registered = False


def cache_path() -> Path:
    override = os.environ.get("FOLDBOUNDS_ORACLE_CACHE")
    return Path(override) if override else _DEFAULT_CACHE


def _entry_ok(seq, g) -> bool:
    if not isinstance(seq, str) or len(seq) < 2 or any(ch not in "HP" for ch in seq):
        return False
    if not isinstance(g, list) or not g:
        return False
    if any((not isinstance(v, int)) or isinstance(v, bool) or v < 0 for v in g):
        return False
    want = expected_conformations(len(seq))
    return want is None or sum(g) == want


def _load_cache() -> dict[str, list[int]]:
    global _cache
    if _cache is None:
        _cache = {}
        try:
            raw = json.loads(cache_path().read_text(encoding="ascii"))
            if isinstance(raw, dict) and raw.get("model") == CACHE_MODEL_TAG:
                for seq, g in raw.get("density_of_states", {}).items():
                    if _entry_ok(seq, g):
                        _cache[seq] = g
        except (OSError, ValueError):
            pass
    return _cache


def flush_cache() -> None:
    """Atomically write the cache (merging with what is on disk)."""
    global _cache_dirty
    if _cache is None or not _cache_dirty:
        return
    path = cache_path()
    merged: dict[str, list[int]] = {}
    try:
        raw = json.loads(path.read_text(encoding="ascii"))
        if isinstance(raw, dict) and raw.get("model") == CACHE_MODEL_TAG:
            for seq, g in raw.get("density_of_states", {}).items():
                if _entry_ok(seq, g):
                    merged[seq] = g
    except (OSError, ValueError):
        pass
    merged.update(_cache)
    payload = {
        "model": CACHE_MODEL_TAG,
        "note": "g[k] = number of anchored conformations with k H-H contacts; "
                "written by experiments/kernel/independent_oracle.py",
        "density_of_states": {s: merged[s] for s in sorted(merged, key=lambda s: (len(s), s))},
    }
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(prefix=".oracle_cache.", suffix=".tmp", dir=str(path.parent))
        with os.fdopen(fd, "w", encoding="ascii") as fh:
            json.dump(payload, fh, separators=(",", ":"))
            fh.write("\n")
        os.replace(tmp, path)
        _cache_dirty = 0
    except OSError:
        pass


def cached_density_of_states(seq: str) -> list[int]:
    """``density_of_states`` with the on-disk cache ``results/oracle_cache.json``."""
    global _cache_dirty, _cache_registered
    validate_sequence(seq)
    cache = _load_cache()
    hit = cache.get(seq)
    if hit is not None:
        return list(hit)
    g = density_of_states(seq)
    want = expected_conformations(len(seq))
    if want is not None and sum(g) != want:
        raise AssertionError(f"oracle count {sum(g)} != c_(n-1)/4 = {want} for {seq}")
    cache[seq] = list(g)
    _cache_dirty += 1
    if not _cache_registered:
        atexit.register(flush_cache)
        _cache_registered = True
    if _cache_dirty >= 64:
        flush_cache()
    return g


# ---------------------------------------------------------------------------
# Exact thermodynamic values
# ---------------------------------------------------------------------------

def exact_values(seq: str, q: Fraction, t: int, use_cache: bool = True) -> dict:
    """Exact ``Z``, ``ZA``, ``ZN`` and ``P`` for weight ``q`` and event ``k >= t``."""
    n = validate_sequence(seq)
    q = Fraction(q)
    if q < 0:
        raise ValueError("q must be nonnegative")
    if not isinstance(t, int) or isinstance(t, bool) or t < 0:
        raise ValueError("t must be a natural number")
    g = cached_density_of_states(seq) if use_cache else density_of_states(seq)
    z = Fraction(0)
    za = Fraction(0)
    zn = Fraction(0)
    power = Fraction(1)
    for k, count in enumerate(g):
        term = count * power
        z += term
        if k >= t:
            za += term
        else:
            zn += term
        power *= q
    return {
        "sequence": seq,
        "n": n,
        "q": q,
        "t": t,
        "conformations": sum(g),
        "density_of_states": g,
        "Z": z,
        "ZA": za,
        "ZN": zn,
        "P": (za / z) if z != 0 else None,
    }


# ---------------------------------------------------------------------------
# Self check
# ---------------------------------------------------------------------------

def self_check(max_n: int = 8, prefix_max_n: int = 6) -> dict:
    """Cross-check the fast enumerator against the naive one.

    For every sequence with ``2 <= n <= max_n``: totals equal ``c_{n-1}/4`` and
    the fast histogram equals the naive ``O(n^2)``-per-conformation histogram.
    For ``n <= prefix_max_n`` the per-prefix statistics are compared as well.
    Raises ``AssertionError`` on the first mismatch.
    """
    sequences = 0
    conformations = 0
    prefixes = 0
    for n in range(2, max_n + 1):
        confs = list(enumerate_conformations(n))
        want = expected_conformations(n)
        if want is not None and len(confs) != want:
            raise AssertionError(f"n={n}: {len(confs)} conformations, expected {want}")
        if len(set(confs)) != len(confs):
            raise AssertionError(f"n={n}: duplicate conformations")
        for letters in product("HP", repeat=n):
            seq = "".join(letters)
            naive = [0] * (n + 2)
            for conf in confs:
                naive[naive_contacts(seq, conf)] += 1
            while len(naive) > 1 and naive[-1] == 0:
                naive.pop()
            fast = density_of_states(seq)
            if fast != naive:
                raise AssertionError(f"{seq}: fast {fast} != naive {naive}")
            if sum(fast) != len(confs):
                raise AssertionError(f"{seq}: total {sum(fast)} != {len(confs)}")
            sequences += 1
            conformations += len(confs)
            if n <= prefix_max_n:
                seen: dict[tuple[Site, ...], tuple[int, Optional[int], Optional[int]]] = {}

                def visit(sites, count, lo, hi, _seen=seen):
                    key = tuple(sites)
                    if key in _seen:
                        raise AssertionError(f"{seq}: prefix visited twice: {key}")
                    _seen[key] = (count, lo, hi)

                root = prefix_stats(seq, visit)
                want_stats = prefix_stats_naive(seq)
                if seen != want_stats:
                    raise AssertionError(f"{seq}: prefix statistics differ from brute force")
                if root != want_stats[((0, 0), (1, 0))]:
                    raise AssertionError(f"{seq}: root statistics differ")
                prefixes += len(seen)
    return {"max_n": max_n, "sequences": sequences,
            "conformations_scored": conformations, "prefixes_compared": prefixes}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _dec(n: int) -> str:
    """Decimal string without tripping CPython's int->str digit limit."""
    if n < 0:
        return "-" + _dec(-n)
    if n.bit_length() <= 8000:
        return str(n)
    half = (n.bit_length() * 30103 // 100000) // 2
    hi, lo = divmod(n, 10 ** half)
    return _dec(hi) + _dec(lo).rjust(half, "0")


def _ratio(x: Optional[Fraction]) -> Optional[str]:
    if x is None:
        return None
    return f"{_dec(x.numerator)}/{_dec(x.denominator)}"


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("sequence", nargs="?", help="HP sequence, e.g. HPHPPHHPHH")
    ap.add_argument("--q", default="1", help="exact rational weight (default 1)")
    ap.add_argument("--t", type=int, default=0, help="event threshold k >= t (default 0)")
    ap.add_argument("--no-cache", action="store_true", help="do not read or write the cache")
    ap.add_argument("--self-check", type=int, metavar="N", nargs="?", const=8,
                    help="cross-check fast vs naive enumeration for all sequences with n <= N")
    args = ap.parse_args(argv)
    if args.self_check is not None:
        print(json.dumps({"self_check": "ok", **self_check(args.self_check)}))
        return 0
    if args.sequence is None:
        ap.error("a sequence is required (or use --self-check)")
    try:
        vals = exact_values(args.sequence, Fraction(args.q), args.t, use_cache=not args.no_cache)
    except (ValueError, ZeroDivisionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    out = {
        "sequence": vals["sequence"],
        "n": vals["n"],
        "conformations": vals["conformations"],
        "density_of_states": vals["density_of_states"],
        "q": _ratio(vals["q"]),
        "t": vals["t"],
        "Z": _ratio(vals["Z"]),
        "ZA": _ratio(vals["ZA"]),
        "ZN": _ratio(vals["ZN"]),
        "P": _ratio(vals["P"]),
    }
    print(json.dumps(out))
    flush_cache()
    return 0


if __name__ == "__main__":
    sys.exit(main())
