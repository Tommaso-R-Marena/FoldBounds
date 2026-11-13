"""Exact integer factor tables and variable elimination (producer side).

Every rational factor ``f`` is stored as an integer table ``f_int`` together
with an exact scale ``s_f = 1 / lcm(denominators of f)`` so that
``f = s_f * f_int``.  Sums of products of scaled tables are therefore exact
integers times a known rational scale; no rounding ever happens.

Tables are flat Python lists in row-major order over their scope (first
variable most significant).  Variable elimination (sum-product) multiplies
all tables that mention a variable, laid out with that variable last, and
sums it out.  Elimination orders are chosen by simulating several greedy
heuristics (min-fill, min-size, and an "inside-out" order by a supplied
priority such as graph distance from the perturbation) and keeping the
cheapest; the choice of order affects only cost, never the exact result.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import lcm, prod
from operator import add, mul
from typing import Iterable, Mapping, Sequence

from .model import Factor


@dataclass
class Table:
    scope: tuple[str, ...]
    dims: tuple[int, ...]
    values: list[int]

    def size(self) -> int:
        return len(self.values)


def factor_to_table(f: Factor, domains: Mapping[str, int]) -> tuple[Table, Fraction]:
    """Integer table and exact scale with ``f == scale * table``."""
    den = 1
    for t in f.table:
        den = lcm(den, t.denominator)
    values = [t.numerator * (den // t.denominator) for t in f.table]
    dims = tuple(domains[v] for v in f.scope)
    return Table(tuple(f.scope), dims, values), Fraction(1, den)


def _strides(dims: Sequence[int]) -> list[int]:
    st = [1] * len(dims)
    for j in range(len(dims) - 2, -1, -1):
        st[j] = st[j + 1] * dims[j + 1]
    return st


def index_map(target_scope: Sequence[str], target_dims: Sequence[int],
              table: Table) -> list[int]:
    """Index into ``table`` for each row-major assignment of ``target_scope``.

    Variables of ``table`` must all occur in ``target_scope``.
    """
    strides = dict(zip(table.scope, _strides(table.dims)))
    idx = [0]
    for v, d in zip(target_scope, target_dims):
        st = strides.get(v, 0)
        if st == 0:
            if d != 1:
                idx = [i for i in idx for _ in range(d)]
        else:
            steps = [x * st for x in range(d)]
            idx = [i + s for i in idx for s in steps]
    return idx


def broadcast(table: Table, target_scope: Sequence[str], target_dims: Sequence[int]) -> list[int]:
    if tuple(target_scope) == table.scope:
        return list(table.values)
    vals = table.values
    return list(map(vals.__getitem__, index_map(target_scope, target_dims, table)))


def multiply_all(tables: Sequence[Table], target_scope: Sequence[str],
                 target_dims: Sequence[int]) -> list[int]:
    n = prod(target_dims)
    out: list[int] | None = None
    for t in tables:
        col = broadcast(t, target_scope, target_dims)
        out = col if out is None else list(map(mul, out, col))
    return out if out is not None else [1] * n


def sum_last(values: list[int], d: int) -> list[int]:
    if d == 1:
        return values
    if d == 2:
        return list(map(add, values[0::2], values[1::2]))
    return list(map(sum, zip(*(values[j::d] for j in range(d)))))


def clamp_table(table: Table, assignment: Mapping[str, int]) -> Table:
    """Restrict a table to fixed values of some of its variables."""
    fixed = [v for v in table.scope if v in assignment]
    if not fixed:
        return table
    strides = dict(zip(table.scope, _strides(table.dims)))
    offset = sum(assignment[v] * strides[v] for v in fixed)
    rest = [(v, d) for v, d in zip(table.scope, table.dims) if v not in assignment]
    idx = [offset]
    for v, d in rest:
        st = strides[v]
        steps = [x * st for x in range(d)]
        idx = [i + s for i in idx for s in steps]
    vals = table.values
    return Table(tuple(v for v, _ in rest), tuple(d for _, d in rest),
                 list(map(vals.__getitem__, idx)))


# ---------------------------------------------------------------- planning
@dataclass(frozen=True)
class Plan:
    order: tuple[str, ...]
    max_size: int
    total_size: int
    heuristic: str


def _simulate(adj: dict[str, set[str]], dims: Mapping[str, int], eliminate: set[str],
              heuristic: str, priority: Mapping[str, float] | None) -> Plan:
    adj = {v: set(n) for v, n in adj.items()}
    remaining = set(eliminate)
    order: list[str] = []
    max_size = 1
    total = 0
    while remaining:
        best = None
        best_key = None
        for v in remaining:
            nb = adj[v]
            size = dims[v]
            for w in nb:
                size *= dims[w]
            if heuristic == "min-size":
                key = (size,)
            else:
                nbl = list(nb)
                fill = 0
                for i in range(len(nbl)):
                    ai = adj[nbl[i]]
                    for j in range(i + 1, len(nbl)):
                        if nbl[j] not in ai:
                            fill += 1
                key = (fill, size) if heuristic == "min-fill" else (size, fill)
            if priority is not None:
                key = (priority.get(v, 0),) + key
            key = key + (v,)
            if best_key is None or key < best_key:
                best_key, best = key, v
        v = best
        nb = adj.pop(v)
        size = dims[v]
        for w in nb:
            size *= dims[w]
        max_size = max(max_size, size)
        total += size
        for w in nb:
            adj[w].discard(v)
            adj[w] |= nb - {w}
        remaining.discard(v)
        order.append(v)
    return Plan(tuple(order), max_size, total, heuristic)


def plan_elimination(scopes: Iterable[Sequence[str]], dims: Mapping[str, int],
                     eliminate: Iterable[str],
                     priority: Mapping[str, float] | None = None,
                     heuristics: Sequence[str] = ("min-fill", "min-size")) -> Plan:
    """Pick the cheapest of several greedy elimination orders (by simulation)."""
    eliminate = set(eliminate)
    adj: dict[str, set[str]] = {}
    for s in scopes:
        for u in s:
            adj.setdefault(u, set()).update(w for w in s if w != u)
    for v in eliminate:
        adj.setdefault(v, set())
    candidates = [_simulate(adj, dims, eliminate, h, None) for h in heuristics]
    if priority is not None:
        candidates += [_simulate(adj, dims, eliminate, h + "/priority", priority)
                       for h in ("min-size", "min-fill")]
    return min(candidates, key=lambda p: (p.total_size + p.max_size, p.max_size))


# ------------------------------------------------------------- elimination
def eliminate(tables: list[Table], order: Sequence[str], dims: Mapping[str, int]) -> list[Table]:
    """Sum out the variables of ``order`` (in that order); returns the rest."""
    pool = list(tables)
    for v in order:
        touching = [t for t in pool if v in t.scope]
        if not touching:
            # Variable without factors contributes a factor d (sum of ones).
            pool.append(Table((), (), [dims[v]]))
            continue
        pool = [t for t in pool if v not in t.scope]
        scope: list[str] = []
        for t in touching:
            for w in t.scope:
                if w != v and w not in scope:
                    scope.append(w)
        scope.append(v)
        tdims = [dims[w] for w in scope]
        values = multiply_all(touching, scope, tdims)
        pool.append(Table(tuple(scope[:-1]), tuple(tdims[:-1]), sum_last(values, dims[v])))
    return pool


def joint_table(tables: Sequence[Table], keep: Sequence[str], dims: Mapping[str, int]) -> list[int]:
    """Product of tables (all scopes inside ``keep``) as one row-major table."""
    kd = [dims[v] for v in keep]
    return multiply_all(tables, list(keep), kd)


def scaled_tables(factors: Iterable[Factor], domains: Mapping[str, int]) -> tuple[list[Table], Fraction]:
    tables, scale = [], Fraction(1)
    for f in factors:
        t, s = factor_to_table(f, domains)
        tables.append(t)
        scale *= s
    return tables, scale


def sum_product(factors: Iterable[Factor], domains: Mapping[str, int],
                keep: Sequence[str] = (), clamp: Mapping[str, int] | None = None,
                eliminate_vars: Iterable[str] | None = None,
                priority: Mapping[str, float] | None = None,
                plan: Plan | None = None) -> tuple[Fraction, list[int], Plan]:
    """Exact ``sum over eliminate_vars of prod factors`` as a table over ``keep``.

    ``clamp`` fixes values of some variables first.  Variables that are neither
    kept, clamped nor in ``eliminate_vars`` must not occur.  Returns
    ``(scale, integer_values, plan)`` with the exact table ``scale * values``.
    """
    tables, scale = scaled_tables(factors, domains)
    if clamp:
        tables = [clamp_table(t, clamp) for t in tables]
    keep = list(keep)
    occurring = {v for t in tables for v in t.scope}
    if eliminate_vars is None:
        elim = occurring - set(keep)
    else:
        elim = set(eliminate_vars)
        stray = occurring - elim - set(keep)
        if stray:
            raise ValueError(f"variables neither kept nor eliminated: {sorted(stray)[:5]}")
    if plan is None:
        plan = plan_elimination([t.scope for t in tables], domains, elim, priority)
    if set(plan.order) != elim:
        raise ValueError("elimination plan does not match the variables to eliminate")
    rest = eliminate(tables, plan.order, domains)
    return scale, joint_table(rest, keep, domains), plan
