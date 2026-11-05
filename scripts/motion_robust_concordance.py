"""Exact-rational interval ranking and directional-concordance checker.

Accepts intervals for nominal model scores that have ALREADY been separately
proved to enclose every allowed structural world. This does not validate the
interval producer, NMR target labels, coordinate bytes or the physical model.

The signed concordance is mean sign agreement over target-nontied pair
comparisons. Target ties are omitted; model-score ties count 0. This is not
Kendall tau-b, a physical confidence interval, or an independence test.
"""
from fractions import Fraction
from itertools import combinations


def Q(x):
    if isinstance(x, bool) or isinstance(x, float):
        raise ValueError("Only exact rational strings, integers, Fraction allowed")
    try:
        z = Fraction(x)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("Invalid rational input") from exc
    return z


def intervals_from_untrusted(data):
    if not isinstance(data, list) or len(data) < 2:
        raise ValueError("Need at least 2 intervals")
    result = []
    for raw in data:
        if not isinstance(raw, (tuple, list)) or len(raw) != 2:
            raise ValueError("Expected [lower,upper] rational interval")
        lo, hi = Q(raw[0]), Q(raw[1])
        if lo > hi:
            raise ValueError("Inverted interval")
        result.append((lo, hi))
    return result


def strict_sign(a, b):
    return (a > b) - (a < b)


def certified_sign(a, b):
    """1 if a is always higher, -1 if always lower, 0 if abstaining."""
    return 1 if a[0] > b[1] else (-1 if a[1] < b[0] else 0)


def check(data, target=None, k=None):
    iv = intervals_from_untrusted(data)
    n = len(iv)
    if k is not None and (type(k) is not int or not (1 <= k <= n)):
        raise ValueError("Invalid k")
    n_pairs = n * (n-1)//2
    certain = sum(certified_sign(iv[i], iv[j]) != 0
                  for i,j in combinations(range(n),2))
    # Descending ranks including every possible tie resolution.
    min_rank = [1 + sum(j != i and iv[j][0] > iv[i][1]
                        for j in range(n)) for i in range(n)]
    max_rank = [n - sum(j != i and iv[j][1] < iv[i][0]
                        for j in range(n)) for i in range(n)]
    result = dict(n=n, pair_count=n_pairs, certain_pair_orders=certain,
                  uncertain_pair_orders=n_pairs-certain,
                  minimum_rank=min_rank, maximum_rank=max_rank)
    if k is not None:
        result["k"] = k
        result["guaranteed_top_k"] = [i for i in range(n)
                                      if max_rank[i] <= k]
        result["possibly_top_k"] = [i for i in range(n)
                                    if min_rank[i] <= k]
    if target is not None:
        if not isinstance(target, list) or len(target) != n:
            raise ValueError("Missing or wrong-length target")
        t = [Q(v) for v in target]
        C = D = U = T = 0
        for i, j in combinations(range(n), 2):
            ts = strict_sign(t[i], t[j])
            if ts == 0:
                T += 1
                continue
            ps = certified_sign(iv[i], iv[j])
            if ps == 0: U += 1
            elif ps == ts: C += 1
            else: D += 1
        m = C+D+U
        if not m:
            raise ValueError("Target has no distinguishable comparisons")
        result.update(target_tied_pairs=T, certified_concordant=C,
                      certified_discordant=D, uncertain_target_pairs=U,
                      concordance_lower=str(Fraction(C-D-U,m)),
                      concordance_upper=str(Fraction(C-D+U,m)))
    return result


def verify(data, claim, target=None):
    """Reject false/incomplete/extra result claims; recompute everything."""
    if not isinstance(claim, dict):
        raise ValueError("Invalid claim")
    expected = check(data, target, claim.get("k"))
    if claim != expected:
        raise ValueError("Untrusted result differs from independent checker")
    return True
