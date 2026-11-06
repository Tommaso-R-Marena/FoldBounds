"""Exact contact-score histogram of disjoint quadrant-monotone SAW families.

An additive polynomial-state dynamic program for genuine extensions of an
anchored square-lattice HP prefix. Quadrant inclusion-exclusion removes
overlapping cardinal rays. This is a model-relative Python research method,
not a Lean-verified program or prediction of physical protein structures.
"""
from __future__ import annotations
from collections import defaultdict
from .lattice import Path, _check_sequence, _contacts, _possible_hh_pairs, _validate_path


def monotone_spectral_histogram(sequence: str, prefix: Path) -> tuple[tuple[int, int], ...]:
    """Count disjoint full contact-score multiplicities, without listing paths.

    In each cone, a signed coordinate potential increases at each step.
    Consequently suffix sites cannot self-intersect or make nonconsecutive
    suffix-to-suffix contacts. Only prefix-suffix contacts need tracking.
    Every straight cardinal ray is in two cones; subtract its histogram once.
    """
    _check_sequence(sequence)
    n, m = len(sequence), len(prefix)
    _validate_path(prefix, n)
    pairs = _possible_hh_pairs(sequence)
    fixed_score = _contacts(prefix, pairs)
    if m == n:
        return ((fixed_score, 1),)
    occupied = frozenset(prefix)
    prefix_contacts = {
        j: tuple(prefix[i] for i, k in pairs if k == j and i < m)
        for j in range(m, n)
    }

    def histogram(directions):
        states = {(prefix[-1], 0): 1}
        for j in range(m, n):
            following = defaultdict(int)
            for ((x, y), cross), multiplicity in states.items():
                for dx, dy in directions:
                    site = (x+dx, y+dy)
                    if site in occupied:
                        continue
                    gain = sum(abs(site[0]-px)+abs(site[1]-py) == 1
                               for px, py in prefix_contacts[j])
                    following[(site, cross+gain)] += multiplicity
            states = following
            if not states:
                break
        counts = defaultdict(int)
        for (_, cross), multiplicity in states.items():
            counts[fixed_score+cross] += multiplicity
        return counts

    E, W, N, S = (1, 0), (-1, 0), (0, 1), (0, -1)
    union = defaultdict(int)
    for cone in ((E, N), (E, S), (W, N), (W, S)):
        for score, multiplicity in histogram(cone).items():
            union[score] += multiplicity
    for ray in ((E,), (W,), (N,), (S,)):
        for score, multiplicity in histogram(ray).items():
            union[score] -= multiplicity
    if any(multiplicity < 0 for multiplicity in union.values()):
        raise AssertionError("negative monotone-cone inclusion-exclusion count")
    return tuple(sorted((score, count) for score, count in union.items() if count))


def is_monotone_cone_suffix(prefix: Path, full_path: Path) -> bool:
    """Classify a valid full path's suffix for witness-overlap exclusion."""
    if len(full_path) <= len(prefix) or full_path[:len(prefix)] != prefix:
        return False
    horizontal, vertical = set(), set()
    for (x, y), (xx, yy) in zip(full_path[len(prefix)-1:],
                                 full_path[len(prefix):]):
        if xx != x:
            horizontal.add(1 if xx > x else -1)
        if yy != y:
            vertical.add(1 if yy > y else -1)
        if len(horizontal) > 1 or len(vertical) > 1:
            return False
    return True
