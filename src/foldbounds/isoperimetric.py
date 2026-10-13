"""Certified-integer isoperimetric contact caps for 2D HP lattice walks.

For a connected finite set S of n square-grid cells, write P for the
number of directed grid edges from S to its complement (edge perimeter).
Every occupied row and column has at least two outgoing boundary edges,
so P >= 2*(width(S)+height(S)). The number of internal grid edges is
2*n - P/2. A self-avoiding chain uses exactly n-1 of these edges as
backbone bonds; therefore the number of nonbonded grid contacts is
at most n+1-width(S)-height(S).

Future conformations containing a given prefix cannot shrink its bounding
box. The final dimensions w,h obey w>=prefix_w, h>=prefix_h and w*h>=n.
Minimizing w+h over these integer dimensions produces a sound upper
bound on *all* nonbonded contacts, hence on the H-H contact count.

These are elementary discrete-geometric bounds, not claimed original.
Python implementation not yet Lean verified.
"""
from __future__ import annotations

from math import isqrt
from .lattice import Path


def universal_contact_cap(n: int) -> int:
    """Upper bound on all noncovalent contacts for a connected n-cell walk."""
    if type(n) is not int or n < 2:
        raise ValueError('n must be an integer >= 2')
    root = isqrt(4*n)
    ceil_root = root if root*root == 4*n else root+1
    return max(0,n+1-ceil_root)


def minimum_final_box_sum(n: int, min_width: int, min_height: int) -> int:
    """Minimum w+h for w*h >= n, w,h at least specified positive integers.

    Search widths 1..n: an n-cell self-avoiding walk can span at most n
    integer columns, and there is always a feasible n-by-1 strip for the
    full length, subject to the specified prefix's dimensions.
    """
    if (type(n) is not int or n < 2 or
        type(min_width) is not int or type(min_height) is not int or
        not 1 <= min_width <= n or not 1 <= min_height <= n):
        raise ValueError('invalid n or minimum dimensions')
    return min(w+max(min_height,(n+w-1)//w)
               for w in range(min_width,n+1))


def prefix_contact_cap(path: Path, n: int) -> int:
    """Maximum possible number of *total* noncovalent contacts in any
    self-avoiding n-site extension of path, using the fixed prefix box.
    H-H contacts are a subset, so the same upper cap applies to them.
    """
    if type(n) is not int or n < len(path) or not path:
        raise ValueError('n must be at least the nonempty prefix length')
    xs=[p[0] for p in path]
    ys=[p[1] for p in path]
    width=max(xs)-min(xs)+1
    height=max(ys)-min(ys)+1
    return max(0,n+1-minimum_final_box_sum(n,width,height))
