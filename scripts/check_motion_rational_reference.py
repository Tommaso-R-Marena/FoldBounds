#!/usr/bin/env python3
"""Dependency-free, independent exact-rational FoldBounds motion regression.

Tests two independently written finite matrix/energy formulas, exact Gaussian
solves, residual enclosures, and score-interval ranking. It is NOT a PDB
decoder audit or a replacement for the Lean kernel proofs.
"""
from fractions import Fraction as Q
import random


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def matvec(a, x):
    return [dot(row, x) for row in a]


def solve(a, z):
    n = len(z)
    t = [row[:] + [z[i]] for i, row in enumerate(a)]
    for k in range(n):
        pivot = next((j for j in range(k, n) if t[j][k]), None)
        assert pivot is not None, "ridge Hessian unexpectedly singular"
        t[k], t[pivot] = t[pivot], t[k]
        d = t[k][k]
        t[k] = [v / d for v in t[k]]
        for j in range(n):
            if j != k and t[j][k]:
                factor = t[j][k]
                t[j] = [x - factor*y for x, y in zip(t[j], t[k])]
    return [row[-1] for row in t]


def hessian_and_edges(coords, tau, cutoff_sq):
    n = len(coords)
    H = [[tau if i == j else Q(0) for j in range(3*n)]
         for i in range(3*n)]
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            dr = tuple(coords[i][k]-coords[j][k] for k in range(3))
            d2 = dot(dr, dr)
            assert d2 > 0, "colliding input must be rejected"
            if d2 > cutoff_sq:
                continue
            edges.append((i, j, dr, d2))
            stiffness = Q(49) / (d2*d2)
            for a in range(3):
                for b in range(3):
                    v = stiffness*dr[a]*dr[b]
                    H[3*i+a][3*i+b] += v
                    H[3*j+a][3*j+b] += v
                    H[3*i+a][3*j+b] -= v
                    H[3*j+a][3*i+b] -= v
    return H, edges


def energy_direct(u, edges, tau):
    total = tau * dot(u, u)
    for i, j, dr, d2 in edges:
        stretch = sum(dr[k]*(u[3*i+k]-u[3*j+k]) for k in range(3))
        total += Q(49)*stretch*stretch/(d2*d2)
    return total


def score_interval(H, z, x, tau, D):
    assert tau > 0 and D > 0
    r = [a-b for a, b in zip(z, matvec(H, x))]
    L = 2*dot(z, x)-dot(x, matvec(H, x))
    U = L + dot(r, r)/tau
    return L/D, U/D


def run():
    rng = random.Random(20261104)
    checked = 0
    ranking_checks = 0
    max_width = Q(0)
    for trial in range(60):
        n = 2 + trial % 4
        coords = []
        while len(coords) < n:
            c = tuple(Q(rng.randrange(-8, 9)) for _ in range(3))
            if c not in coords:
                coords.append(c)
        tau = Q(1, 5)
        H, edges = hessian_and_edges(coords, tau, Q(100))
        for _ in range(3):
            u = [Q(rng.randint(-3, 3)) for _ in range(3*n)]
            assert dot(u, matvec(H, u)) == energy_direct(u, edges, tau)
        scores = []
        for i, j in [(0, 1), (0, n-1)]:
            dr = [coords[i][k]-coords[j][k] for k in range(3)]
            D = dot(dr, dr)
            assert D > 0
            z = [Q(0)]*(3*n)
            for k in range(3):
                z[3*i+k] = dr[k]
                z[3*j+k] = -dr[k]
            y = solve(H, z)
            assert matvec(H, y) == z
            q = dot(z, y)/D
            approx = [Q(round(float(a)*10**5), 10**5) for a in y]
            lo, hi = score_interval(H, z, approx, tau, D)
            assert lo <= q <= hi
            max_width = max(max_width, hi-lo)
            scores.append((lo, hi, q))
            checked += 1
        (l1, u1, q1), (l2, u2, q2) = scores
        if u1 < l2:
            assert q1 < q2
            ranking_checks += 1
        elif u2 < l1:
            assert q2 < q1
            ranking_checks += 1
    print("PASS: 180 exact Hessian/energy identities")
    print(f"PASS: {checked} exact solves and response intervals")
    print(f"PASS: {ranking_checks} disjoint-interval rank comparisons")
    print(f"Max rational response-interval width: {float(max_width):.3g}")


if __name__ == "__main__":
    run()
