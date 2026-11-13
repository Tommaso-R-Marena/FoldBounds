"""Coupled factor-graph (Markov random field) model holding two systems A and B.

Semantics
---------
Each variable ``v`` has a finite domain ``{0, ..., d_X(v) - 1}`` in system
X in {A, B}.  A *factor* is a scope (a tuple of distinct variable names) and a
table of nonnegative rationals stored in row-major order over the scope (the
first scope variable is the most significant digit).  The weight of a joint
state x of system X is the product of all factors of X evaluated at x, and

    Z_X = sum over joint states x of prod_{f in factors(X)} f(x_scope(f)).

The coupled model stores three factor lists:

* ``shared``  -- factors present, with identical tables, in BOTH systems;
* ``only_a``  -- factors present only in system A;
* ``only_b``  -- factors present only in system B;

so ``factors(A) = shared + only_a`` and ``factors(B) = shared + only_b``.
Both systems have the same variable names.  A variable may have different
domain sizes in A and B (for example a mutated residue with a different
rotamer set); such a variable may only appear in system-specific factors.

The *perturbation support* ``P`` is the union of the scopes of all
system-specific factors together with all variables whose domain differs.

Every number is an exact :class:`fractions.Fraction`.  No floating-point value
is ever accepted as model data.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import prod
from typing import Iterable, Mapping, Sequence


class ModelError(ValueError):
    """Raised for malformed coupled factor-graph models."""


def as_rational(value) -> Fraction:
    """Convert an int, Fraction or decimal/rational string to a Fraction.

    Floats (and bools) are rejected on purpose: model data must be exact.
    """
    if isinstance(value, bool):
        raise ModelError(f"boolean is not a rational number: {value!r}")
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, Fraction):
        return value
    if isinstance(value, str):
        text = value.strip()
        if not text or any(c in text for c in "eEjJ") or text.lower() in {"nan", "inf", "-inf", "+inf"}:
            raise ModelError(f"not an exact rational literal: {value!r}")
        try:
            return Fraction(text)
        except (ValueError, ZeroDivisionError) as exc:
            raise ModelError(f"not an exact rational literal: {value!r}") from exc
    raise ModelError(f"unsupported numeric type {type(value).__name__}: {value!r}")


@dataclass(frozen=True)
class Factor:
    """A nonnegative rational factor over an ordered scope of variable names."""

    scope: tuple[str, ...]
    table: tuple[Fraction, ...]
    name: str = ""

    def __post_init__(self):
        scope = tuple(self.scope)
        if len(set(scope)) != len(scope):
            raise ModelError(f"factor {self.name!r} has a repeated variable in scope {scope}")
        for v in scope:
            if not isinstance(v, str) or not v:
                raise ModelError(f"factor {self.name!r}: variable names must be nonempty strings")
        table = tuple(as_rational(t) for t in self.table)
        for t in table:
            if t < 0:
                raise ModelError(f"factor {self.name!r} has a negative entry {t}")
        object.__setattr__(self, "scope", scope)
        object.__setattr__(self, "table", table)

    def check_dims(self, domains: Mapping[str, int]) -> None:
        for v in self.scope:
            if v not in domains:
                raise ModelError(f"factor {self.name!r} uses unknown variable {v!r}")
        need = prod(domains[v] for v in self.scope)
        if len(self.table) != need:
            raise ModelError(
                f"factor {self.name!r} over {self.scope} needs {need} entries, got {len(self.table)}")

    def value(self, assignment: Mapping[str, int], domains: Mapping[str, int]) -> Fraction:
        index = 0
        for v in self.scope:
            index = index * domains[v] + assignment[v]
        return self.table[index]

    def key(self) -> tuple:
        return (self.scope, self.table)


@dataclass(frozen=True)
class SystemView:
    """One system (A or B) of a coupled model, as a plain factor graph."""

    domains: Mapping[str, int]
    factors: tuple[Factor, ...]


@dataclass(frozen=True)
class CoupledModel:
    """Two factor graphs A and B over the same variable names.

    ``variables`` fixes a canonical variable order (used for state encodings).
    """

    variables: tuple[str, ...]
    domains_a: Mapping[str, int]
    domains_b: Mapping[str, int]
    shared: tuple[Factor, ...] = ()
    only_a: tuple[Factor, ...] = ()
    only_b: tuple[Factor, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self):
        variables = tuple(self.variables)
        if len(set(variables)) != len(variables):
            raise ModelError("repeated variable name")
        da = dict(self.domains_a)
        db = dict(self.domains_b)
        if set(da) != set(variables) or set(db) != set(variables):
            raise ModelError("domains_a/domains_b must list exactly the model variables")
        for d in (da, db):
            for v, size in d.items():
                if isinstance(size, bool) or not isinstance(size, int) or size < 1:
                    raise ModelError(f"variable {v!r} has invalid domain size {size!r}")
        object.__setattr__(self, "variables", variables)
        object.__setattr__(self, "domains_a", da)
        object.__setattr__(self, "domains_b", db)
        object.__setattr__(self, "shared", tuple(self.shared))
        object.__setattr__(self, "only_a", tuple(self.only_a))
        object.__setattr__(self, "only_b", tuple(self.only_b))
        changed = self.domain_changed()
        for f in self.shared:
            f.check_dims(da)
            bad = [v for v in f.scope if v in changed]
            if bad:
                raise ModelError(
                    f"shared factor {f.name!r} touches variables {bad} whose domain differs between A and B")
        for f in self.only_a:
            f.check_dims(da)
        for f in self.only_b:
            f.check_dims(db)

    # ------------------------------------------------------------------ views
    def domain_changed(self) -> frozenset[str]:
        return frozenset(v for v in self.variables if self.domains_a[v] != self.domains_b[v])

    def system(self, which: str) -> SystemView:
        if which == "A":
            return SystemView(self.domains_a, self.shared + self.only_a)
        if which == "B":
            return SystemView(self.domains_b, self.shared + self.only_b)
        raise ValueError("system must be 'A' or 'B'")

    def perturbation_support(self) -> frozenset[str]:
        support = set(self.domain_changed())
        for f in self.only_a + self.only_b:
            support.update(f.scope)
        return frozenset(support)

    def adjacency(self) -> dict[str, set[str]]:
        """Primal (interaction) graph over the union of A and B factor scopes."""
        adj: dict[str, set[str]] = {v: set() for v in self.variables}
        for f in self.shared + self.only_a + self.only_b:
            for u in f.scope:
                adj[u].update(w for w in f.scope if w != u)
        return adj

    def num_states(self, which: str = "A") -> int:
        d = self.domains_a if which == "A" else self.domains_b
        return prod(d.values())


def make_model(domains_a: Mapping[str, int], shared: Iterable[Factor] = (),
               only_a: Iterable[Factor] = (), only_b: Iterable[Factor] = (),
               domains_b: Mapping[str, int] | None = None,
               variables: Sequence[str] | None = None,
               metadata: Mapping[str, object] | None = None) -> CoupledModel:
    """Convenience constructor; ``domains_b`` defaults to ``domains_a``."""
    if domains_b is None:
        domains_b = dict(domains_a)
    if variables is None:
        variables = tuple(domains_a)
    return CoupledModel(tuple(variables), dict(domains_a), dict(domains_b),
                        tuple(shared), tuple(only_a), tuple(only_b), dict(metadata or {}))


def from_two_systems(domains_a: Mapping[str, int], factors_a: Iterable[Factor],
                     domains_b: Mapping[str, int], factors_b: Iterable[Factor],
                     variables: Sequence[str] | None = None) -> CoupledModel:
    """Build a coupled model from two full factor lists by multiset difference.

    A factor (scope and table equal) present in both lists with multiplicity m
    in each is shared m times; the surplus copies become system-specific.
    Factors touching a variable whose domain differs are never shared.
    """
    if variables is None:
        variables = tuple(domains_a)
    changed = {v for v in variables if domains_a[v] != domains_b[v]}
    pool: dict[tuple, list[Factor]] = {}
    for f in factors_b:
        pool.setdefault(f.key(), []).append(f)
    shared, only_a = [], []
    for f in factors_a:
        bucket = pool.get(f.key())
        if bucket and not (set(f.scope) & changed):
            bucket.pop()
            shared.append(f)
        else:
            only_a.append(f)
    only_b = [f for bucket in pool.values() for f in bucket]
    return CoupledModel(tuple(variables), dict(domains_a), dict(domains_b),
                        tuple(shared), tuple(only_a), tuple(only_b))
