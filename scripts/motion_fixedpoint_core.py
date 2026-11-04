"""Exact integer/outward fixed-point interval kernel for FoldBounds motion.

Coordinates are exact integers / 2**24. A proposed linear-solve witness is
exact integers / 2**30. Every matrix entry and residual interval uses 2**96
fixed-point endpoints, rounded *outward* with integer floor/ceiling.

Untrusted SciPy never enters this module. The independent Lean implementation
of this specific integer rounding scheme remains a future trust boundary.
"""
D = 2 ** 24
B = 2 ** 96
Y_SCALE = 2 ** 30

def floorq(a, b):
    if b <= 0:
        raise ValueError("nonpositive denominator")
    return a // b

def ceilq(a, b):
    if b <= 0:
        raise ValueError("nonpositive denominator")
    return -((-a) // b)

def add(a, b):
    return a[0] + b[0], a[1] + b[1]

def negate(a):
    return -a[1], -a[0]

def mul_int_interval(a, y):
    """a is [lo,hi]/B, y is integer y/B; return rounded product /B."""
    if y >= 0:
        return floorq(a[0] * y, B), ceilq(a[1] * y, B)
    return floorq(a[1] * y, B), ceilq(a[0] * y, B)

def hessian_boxes(coords):
    """Return complete sparse H intervals and active contact count.

    Contacts include all i<j with nonzero squared distance <=100 Å².
    Entries use 49 * D² * delta_int_a * delta_int_b / norm_int^4,
    with the exact +1/5 ridge placed on every Cartesian diagonal.
    """
    n = len(coords)
    if n < 2 or any(len(v) != 3 for v in coords):
        raise ValueError("malformed coordinates")
    if any(not isinstance(q, int) for v in coords for q in v):
        raise ValueError("noninteger coordinates")
    rows = [{} for _ in range(3*n)]
    active = 0

    def inc(i, j, term):
        rows[i][j] = add(rows[i].get(j, (0, 0)), term)

    for a in range(n):
        for b in range(a+1, n):
            v = [coords[a][k] - coords[b][k] for k in range(3)]
            s = sum(t*t for t in v)
            if s <= 0:
                raise ValueError("colliding residues")
            if s > 100*D*D:
                continue
            active += 1
            for k in range(3):
                for l in range(3):
                    num = 49*D*D*v[k]*v[l]
                    den = s*s
                    term = floorq(num*B, den), ceilq(num*B, den)
                    ia, ja = 3*a+k, 3*a+l
                    ib, jb = 3*b+k, 3*b+l
                    inc(ia, ja, term)
                    inc(ib, jb, term)
                    inc(ia, jb, negate(term))
                    inc(ib, ja, negate(term))
    ridge = floorq(B, 5), ceilq(B, 5)
    for i in range(3*n):
        inc(i, i, ridge)
    assert all(lo <= hi for row in rows for lo, hi in row.values())
    return rows, active

def bound_score(rows, pair, coords, proposal):
    """Rigorous normalized ANM response [lower, upper]/2**96.

    The approximate linear solve is a vector of integers /2**30 and is
    NEVER trusted. Recompute the full residual against outward-rounded
    exact rational Hessian entry intervals.
    """
    n = len(coords)
    if len(rows) != 3*n or len(proposal) != 3*n:
        raise ValueError("wrong dimension")
    a, b = pair
    if not (0 <= a < b < n):
        raise ValueError("invalid 0-based residue-pair position")
    v = [coords[a][k] - coords[b][k] for k in range(3)]
    sq = sum(t*t for t in v)
    if sq <= 0:
        raise ValueError("zero-length pair")
    z = [0]*(3*n)
    for k, vk in enumerate(v):
        z[3*a+k] = vk*(B//D)
        z[3*b+k] = -vk*(B//D)
    y = [int(t)*(B//Y_SCALE) for t in proposal]
    Hy = []
    for row in rows:
        lo, hi = 0, 0
        for j, entry in row.items():
            term = mul_int_interval(entry, y[j])
            lo += term[0]
            hi += term[1]
        Hy.append((lo,hi))
    zTy = sum(a*b for a,b in zip(z,y))
    zl, zu = floorq(zTy,B), ceilq(zTy,B)
    ql, qu = 0, 0
    for yi, hy in zip(y,Hy):
        q = mul_int_interval(hy, yi)
        ql += q[0]
        qu += q[1]
    Llo, Lhi = 2*zl-qu, 2*zu-ql
    R = sum(max(abs(z[k]-Hy[k][0]), abs(z[k]-Hy[k][1]))**2
            for k in range(3*n))
    Rhi = ceilq(5*R, B)  # tau^-1=5, residual squared
    return floorq(Llo*D*D, sq), ceilq((Lhi+Rhi)*D*D, sq)
