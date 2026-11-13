"""Separators (Markov blankets) for coupled factor graphs.

A separator partitions the variables into interior I, boundary dB and exterior
E such that

  (S1) no factor of A or B touches both I and E;
  (S2) every system-specific factor (only_A / only_B) has scope inside I u dB;
  (S3) every variable whose domain differs between A and B lies in I.

Under (S1)-(S3), the *local* factors (scope inside I u dB) and the *exterior*
factors (all others; they touch E, never I, and are shared by (S2)) give the
factorisation of CEC_THEORY.md section 2.1 with

    Phi_X(tau) = sum_{x_I} prod_{local f of X} f,   Psi(tau) = sum_{x_E} prod_{exterior f} f,

so that ``Z_X = sum_tau Phi_X(tau) Psi(tau)`` (Lemma 2.0).  Shared factors
whose scope lies inside dB are put on the local side: this never changes
rho(tau) where they are positive and makes the states where they vanish inert,
so the local enclosure can only get tighter (and stays sound).
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Iterable, Mapping

from .model import CoupledModel, Factor


class SeparatorError(ValueError):
    """Raised when a proposed separator violates (S1)-(S3)."""


@dataclass(frozen=True)
class Separator:
    interior: frozenset[str]
    boundary: tuple[str, ...]
    exterior: frozenset[str]

    @staticmethod
    def of(interior: Iterable[str], boundary: Iterable[str], exterior: Iterable[str]) -> "Separator":
        return Separator(frozenset(interior), tuple(boundary), frozenset(exterior))

    def to_dict(self, model: CoupledModel | None = None) -> dict:
        order = {v: i for i, v in enumerate(model.variables)} if model else None
        key = (lambda v: order[v]) if order else (lambda v: v)
        return {"interior": sorted(self.interior, key=key), "boundary": list(self.boundary),
                "exterior": sorted(self.exterior, key=key)}

    @staticmethod
    def from_dict(d: Mapping) -> "Separator":
        return Separator.of(d["interior"], d["boundary"], d["exterior"])


def separator_violations(model: CoupledModel, sep: Separator) -> list[str]:
    """Return a list of human-readable violations; empty means valid."""
    problems: list[str] = []
    boundary = set(sep.boundary)
    if len(boundary) != len(sep.boundary):
        problems.append("boundary lists a variable twice")
    parts = (sep.interior, boundary, sep.exterior)
    allv = set(model.variables)
    seen: set[str] = set()
    for part in parts:
        unknown = set(part) - allv
        if unknown:
            problems.append(f"unknown variables {sorted(unknown)}")
        overlap = seen & set(part)
        if overlap:
            problems.append(f"variables in two parts: {sorted(overlap)}")
        seen |= set(part)
    missing = allv - seen
    if missing:
        problems.append(f"variables in no part: {sorted(missing)[:10]}")
    for f in model.shared + model.only_a + model.only_b:
        s = set(f.scope)
        if s & sep.interior and s & sep.exterior:
            problems.append(f"factor {f.name or f.scope} touches interior and exterior")
    local = sep.interior | boundary
    for f in model.only_a + model.only_b:
        if not set(f.scope) <= local:
            problems.append(f"system-specific factor {f.name or f.scope} leaves interior+boundary")
    for v in model.domain_changed():
        if v not in sep.interior:
            problems.append(f"variable {v!r} has a different domain in A and B but is not interior")
    return problems


def require_separator(model: CoupledModel, sep: Separator) -> None:
    problems = separator_violations(model, sep)
    if problems:
        raise SeparatorError("; ".join(problems[:5]))


def graph_distances(adjacency: Mapping[str, set[str]], sources: Iterable[str]) -> dict[str, int]:
    """Breadth-first graph distance from a set of sources (unreachable: absent)."""
    dist: dict[str, int] = {}
    queue: deque[str] = deque()
    for s in sources:
        if s not in dist:
            dist[s] = 0
            queue.append(s)
    while queue:
        u = queue.popleft()
        for w in adjacency[u]:
            if w not in dist:
                dist[w] = dist[u] + 1
                queue.append(w)
    return dist


def radius_separator(model: CoupledModel, radius: int,
                     core: Iterable[str] | None = None,
                     adjacency: Mapping[str, set[str]] | None = None) -> Separator:
    """Shell separator around the perturbation support (or a given core).

    interior = {d < R}, boundary = {d = R}, exterior = {d > R} (and unreachable),
    where d is the interaction-graph distance from the core.  Every factor
    scope is a clique of the interaction graph, so it spans at most two
    consecutive shells; hence (S1) holds.  (S2) holds because system-specific
    factors lie in shell 0.  (S3) can fail at R = 0 when a domain changes;
    callers should then start at R = 1 (``require_separator`` reports it).
    """
    if radius < 0:
        raise ValueError("radius must be nonnegative")
    adj = adjacency if adjacency is not None else model.adjacency()
    core = list(core) if core is not None else sorted(model.perturbation_support(),
                                                       key=model.variables.index)
    dist = graph_distances(adj, core)
    order = {v: i for i, v in enumerate(model.variables)}
    interior = frozenset(v for v, d in dist.items() if d < radius)
    boundary = tuple(sorted((v for v, d in dist.items() if d == radius), key=order.__getitem__))
    exterior = frozenset(v for v in model.variables if v not in interior and v not in set(boundary))
    return Separator(interior, boundary, exterior)


def split_factors(factors: Iterable[Factor], sep: Separator) -> tuple[list[Factor], list[Factor]]:
    """Split a factor list into local (scope in I u dB) and exterior factors."""
    local_vars = sep.interior | set(sep.boundary)
    local, ext = [], []
    for f in factors:
        (local if set(f.scope) <= local_vars else ext).append(f)
    return local, ext
