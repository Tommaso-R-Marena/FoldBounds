"""Conservative all-q sequence-dominance certificates for exact HP polynomials.

Coefficientwise dominance in t=q-1 is sufficient for global Z dominance for
all real q>=1, but not necessary. Absence of a certificate is not evidence
that the ordering ever reverses.
"""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations


def _validate(landscape, h_count):
    n = landscape.length
    if type(h_count) is not int or not 0 <= h_count <= n:
        raise ValueError("h_count must be an integer between 0 and n")
    if n > 16:
        raise ValueError("refusing exhaustive all-sequence scan for n > 16")


def sequence_from_mask(mask, n):
    if type(n) is not int or n < 2 or type(mask) is not int or not 0 <= mask < (1 << n):
        raise ValueError("invalid length or sequence mask")
    return "".join("H" if mask & (1 << i) else "P" for i in range(n))


def dominates(a, b):
    """Sufficient exact inequality for Z_a(q)>=Z_b(q), ALL q>=1.

    Returns True even for equal polynomials. Both share same zero-order
    partition coefficient in the same lattice length, but no such condition
    is required by this low-level inequality check.
    """
    if not isinstance(a, tuple) or not isinstance(b, tuple) or not a or not b:
        raise ValueError("expected nonempty polynomial coefficient tuples")
    if any(type(i) is not int or i < 0 for i in a + b):
        raise ValueError("coefficients must be nonnegative integers")
    return all(
        (a[i] if i < len(a) else 0) >= (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    )


@dataclass(frozen=True)
class DominanceAtlas:
    length: int
    h_count: int
    sequences_checked: int
    universally_dominant: tuple[str, ...]
    coefficientwise_nondominated: tuple[str, ...]
    dominated_count: int


def certify_fixed_composition_dominance(landscape, h_count):
    """Enumerate all H/P sequences at composition k and certify domination.

    A sequence listed in universally_dominant has polynomial coefficients
    >= every other k-H sequence: hence it maximizes Z for all real q>=1.
    The nondominated set can include false candidates, because the
    coefficientwise order is sufficient rather than necessary.
    """
    _validate(landscape, h_count)
    n = landscape.length
    masks = tuple(sum(1 << i for i in c)
                  for c in combinations(range(n), h_count))
    polys = {
        m: tuple(landscape.polynomial(sequence_from_mask(m, n)))
        for m in masks
    }
    maxima = []
    undecided = []
    for m in masks:
        a = polys[m]
        if all(dominates(a, polys[r]) for r in masks):
            maxima.append(sequence_from_mask(m, n))
        if not any(
            r != m and dominates(polys[r], a) and not dominates(a, polys[r])
            for r in masks
        ):
            undecided.append(sequence_from_mask(m, n))
    return DominanceAtlas(
        n, h_count, len(masks), tuple(maxima), tuple(undecided),
        len(masks) - len(undecided)
    )
