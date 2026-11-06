"""Proof-carrying *joint* GNM rank intervals under contact uncertainty.

For connected contact graphs L_off <= L <= L_on, their zero-mean GNM
covariances satisfy C_on <= C <= C_off in Loewner order. Comparing
individual marginal bounds discards the cross-covariance of the *difference*
D = C_off-C_on. This module certifies the strictly sharper, rank-two
spectral envelope for C_ii-C_jj, using exact integer arithmetic only:

    C_on[ii]-C_on[jj] + (Dii-Djj ± sqrt((Dii+Djj)^2-4Dij^2))/2.

The producer of inverse matrices is UNTRUSTED. The checker first verifies
B*N=d*I entry by entry, B=n*L+11^T, n<=128 and d>0. No floating-point
inverses, eigensolvers, square roots or approximate residuals enter an
accepted certificate. This is standard inverse-monotonicity and rank-two
semidefinite optimization, not a claim of mathematical priority.

GNM is a simplified Gaussian network, not experimentally validated by this
module. Coordinate boxes are hypothetical independent Cartesian boxes.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Iterable

Edge = tuple[int, int]


def canonical_edges(n: int, edges: Iterable[Edge]) -> tuple[Edge, ...]:
    if type(n) is not int or not 2 <= n <= 128:
        raise ValueError("n must be an integer from 2 through 128")
    found = set()
    for edge in edges:
        if (not isinstance(edge, tuple) or len(edge) != 2
                or type(edge[0]) is not int or type(edge[1]) is not int):
            raise ValueError("edges must be integer (i,j) pairs")
        i, j = edge
        if not 0 <= i < j < n or edge in found:
            raise ValueError("edges must be sorted individually, distinct and in range")
        found.add(edge)
    for i in range(n - 1):
        if (i, i + 1) not in found:
            raise ValueError("all covalent backbone edges are mandatory")
    return tuple(sorted(found))


def scaled_grounded_laplacian(n: int, edges: Iterable[Edge]) -> tuple[tuple[int, ...], ...]:
    """B=nL+11^T, a symmetric positive-definite INTEGER matrix."""
    edges = canonical_edges(n, edges)
    matrix = [[1] * n for _ in range(n)]
    for i, j in edges:
        matrix[i][i] += n
        matrix[j][j] += n
        matrix[i][j] -= n
        matrix[j][i] -= n
    return tuple(tuple(row) for row in matrix)


@dataclass(frozen=True)
class IntegerInverse:
    """Untrusted proposed inverse N/d for a grounded integer matrix."""
    numerator: tuple[tuple[int, ...], ...]
    denominator: int


def verify_inverse(n: int, edges: Iterable[Edge], witness: IntegerInverse) -> bool:
    """Check all n^2 integer identities B*N=dI; never trust the producer."""
    b = scaled_grounded_laplacian(n, edges)
    if not isinstance(witness, IntegerInverse) or type(witness.denominator) is not int:
        return False
    if witness.denominator <= 0 or not isinstance(witness.numerator, tuple) or len(witness.numerator) != n:
        return False
    m = witness.numerator
    if any(not isinstance(row, tuple) or len(row) != n
           or any(type(x) is not int for x in row) for row in m):
        return False
    if any(m[i][j] != m[j][i] for i in range(n) for j in range(i + 1, n)):
        return False
    # Separate integer multiplication, no SymPy and no untrusted inverses.
    for i in range(n):
        for j in range(n):
            if sum(b[i][k] * m[k][j] for k in range(n)) != (witness.denominator if i == j else 0):
                return False
    return True


def untrusted_propose_inverse(n: int, edges: Iterable[Edge]) -> IntegerInverse:
    """Optional SymPy DomainMatrix producer; outputs ONLY untrusted integers.

    SymPy is intentionally not imported by any checker. Install it to
    generate real-data witnesses; the verification path needs only stdlib.
    """
    from sympy import ZZ, Matrix
    from sympy.polys.matrices import DomainMatrix

    mat = Matrix(scaled_grounded_laplacian(n, edges))
    proposed, denom = DomainMatrix.from_Matrix(mat).convert_to(ZZ).inv_den()
    p = proposed.to_Matrix()
    d = int(denom)
    if d < 0:
        p = -p
        d = -d
    vals = [[int(p[i, j]) for j in range(n)] for i in range(n)]
    common = d
    for row in vals:
        for v in row:
            common = gcd(common, abs(v))
    d //= common
    vals = [tuple(v // common for v in row) for row in vals]
    witness = IntegerInverse(tuple(vals), d)
    if not verify_inverse(n, edges, witness):
        raise ArithmeticError("untrusted proposal failed exact integer inverse check")
    return witness


@dataclass(frozen=True)
class RankCounts:
    total_pairs: int
    simple_certified: int
    joint_certified: int
    gained: int
    certified_orders: tuple[tuple[int, int, int], ...]
    # signed order = +1 if variance at i is strictly greater than at j.


def joint_rank_certificate(
    n: int, certain_edges: Iterable[Edge], possible_edges: Iterable[Edge],
    certain_inverse: IntegerInverse, possible_inverse: IntegerInverse,
    *, exact_verification: bool = True,
) -> RankCounts:
    """Exact rank-two Loewner interval certificate for every residue pair.

    The two inverse inputs must match B_certain and B_possible. A graph
    with any subset of optional edges is covered. If the two supplied
    matrices fail exact integer verification, refuse all output.
    """
    lo_edges = canonical_edges(n, certain_edges)
    hi_edges = canonical_edges(n, possible_edges)
    if not set(lo_edges) <= set(hi_edges):
        raise ValueError("certain contact edges must be a subset of possible edges")
    if exact_verification is not True:
        raise ValueError("independent exact verification cannot be disabled")
    if not verify_inverse(n, lo_edges, certain_inverse):
        raise ValueError("invalid certain-graph integer inverse witness")
    if not verify_inverse(n, hi_edges, possible_inverse):
        raise ValueError("invalid possible-graph integer inverse witness")

    den_lo = certain_inverse.denominator
    den_hi = possible_inverse.denominator
    num_lo = certain_inverse.numerator
    num_hi = possible_inverse.numerator
    simple = joint = gained = 0
    decisions = []
    for i in range(n):
        for j in range(i + 1, n):
            # Multiplicative n and subtractive 1/n translate-mode terms
            # cancel from these *sign* comparisons.
            simple_sign = (
                +1 if num_hi[i][i] * den_lo > num_lo[j][j] * den_hi else
                -1 if num_hi[j][j] * den_lo > num_lo[i][i] * den_hi else 0
            )
            # Work in one common positive denominator den_hi*den_lo.
            # D = C_certain - C_possible >= 0 in Loewner order.
            di = num_lo[i][i] * den_hi - num_hi[i][i] * den_lo
            dj = num_lo[j][j] * den_hi - num_hi[j][j] * den_lo
            dij = num_lo[i][j] * den_hi - num_hi[i][j] * den_lo
            if di < 0 or dj < 0 or di * dj < dij * dij:
                raise ArithmeticError("PSD covariance difference postcondition failed")
            f = (num_hi[i][i] - num_hi[j][j]) * den_lo
            s = 2 * f + di - dj
            radical_sq = (di + dj) ** 2 - 4 * dij ** 2
            if radical_sq < 0:
                raise ArithmeticError("negative rank-two discriminant")
            joint_sign = (+1 if s > 0 and s * s > radical_sq else
                          -1 if s < 0 and s * s > radical_sq else 0)
            if simple_sign:
                simple += 1
                if joint_sign != simple_sign:
                    raise ArithmeticError("joint bound failed to subsume diagonal bound")
            if joint_sign:
                joint += 1
                decisions.append((i, j, joint_sign))
                if not simple_sign:
                    gained += 1
    return RankCounts(n * (n - 1) // 2, simple, joint, gained, tuple(decisions))


def classify_contact_boxes(
    coordinates_milliangstrom: tuple[tuple[int, int, int], ...],
    radius_angstrom: Fraction, coordinate_error_angstrom: Fraction,
) -> tuple[tuple[Edge, ...], tuple[Edge, ...]]:
    """SOURCE ONLY: classify certain/possible contacts with exact integers.

    Each Cartesian coordinate of each residue may change independently by
    at most eps. The inter-residue distance shifts by <=2sqrt(3)*eps,
    which is strictly < 7eps/2. The rational 7/2 bound is conservative.
    Contacts use strict d<radius, and nonlocal pairs are |i-j|>=3.
    """
    coords = coordinates_milliangstrom
    n = len(coords)
    if type(radius_angstrom) is not Fraction or radius_angstrom <= 0:
        raise ValueError("contact radius must be an exact positive Fraction")
    if type(coordinate_error_angstrom) is not Fraction or coordinate_error_angstrom < 0:
        raise ValueError("coordinate error must be an exact nonnegative Fraction")
    if any(not isinstance(p, tuple) or len(p) != 3
           or any(type(v) is not int for v in p) for p in coords):
        raise ValueError("coordinates must be PDB-derived integer milliangstroms")
    # In milliangstroms; the scaled distance test stays purely rational.
    radius = radius_angstrom * 1000
    slack = Fraction(7, 2) * coordinate_error_angstrom * 1000
    lower_radius = max(Fraction(0), radius - slack)
    upper_radius = radius + slack
    certain = {(i, i+1) for i in range(n-1)}
    possible = set(certain)
    for i in range(n):
        for j in range(i+3, n):
            squared = sum((coords[i][k]-coords[j][k])**2 for k in range(3))
            if squared < lower_radius * lower_radius:
                certain.add((i,j))
            if squared < upper_radius * upper_radius:
                possible.add((i,j))
    return canonical_edges(n, certain), canonical_edges(n, possible)


def compare_certificate_to_labels(
    certificate: RankCounts, observed_values: tuple[float, ...]
) -> dict[str, float | int]:
    """Exploratory descriptive concordance; no calibration/performance claim."""
    if len(observed_values) * (len(observed_values) - 1) // 2 != certificate.total_pairs:
        raise ValueError("observed values must match residue population")
    agreeing = sum(1 for i, j, s in certificate.certified_orders
                   if (observed_values[i] - observed_values[j]) * s > 0)
    disagreeing = sum(1 for i, j, s in certificate.certified_orders
                      if (observed_values[i] - observed_values[j]) * s < 0)
    return {"certified_pairs": certificate.joint_certified,
            "agreement": agreeing, "disagreement": disagreeing,
            "agreement_fraction": agreeing / certificate.joint_certified if certificate.joint_certified else float('nan')}
