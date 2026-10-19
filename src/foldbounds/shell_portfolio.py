"""Exact verified-by-recomputation budgeted antichains of extremal-quadrant shells.

Scientific trust boundary: Python reference logic with independent tests, not
formally verified by the Lean kernel. The planner is untrusted and separate
from the checker. Leaf-indexed shell families are disjoint by prefix ancestry.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Literal

from .lattice import (Path, CertificateError, _check_q, _check_sequence,
                      _neighbors, _validate_path, ROOT)
from .shell_family import extremal_shell_mass

Objective = Literal['event', 'complement']

@dataclass(frozen=True)
class PortfolioCertificate:
    sequence: str
    prefix: Path
    q: Fraction
    event_at_least: int
    leaves: tuple[Path, ...]

@dataclass(frozen=True)
class PortfolioMass:
    event: Fraction
    complement: Fraction
    distinct_paths: int
    leaf_count: int

    @property
    def partition(self) -> Fraction:
        return self.event + self.complement

@dataclass(frozen=True)
class PortfolioPlan:
    certificate: PortfolioCertificate
    mass: PortfolioMass
    objective: Objective
    max_depth: int
    leaf_budget: int
    candidate_nodes: int


def verify_shell_portfolio(cert: PortfolioCertificate) -> PortfolioMass:
    """Fail-closed rechecking: no proposer-supplied weights or claimed scores.

    Covers the full prefix tree exactly once by a prefix-free antichain; each
    selected leaf contributes an exactly counted restricted shell family.
    The family may be incomplete relative to ALL valid SAWs, so this returns
    rigorous lower masses, not the full partition mass.
    """
    _check_sequence(cert.sequence)
    _check_q(cert.q)
    n = len(cert.sequence)
    if type(cert.event_at_least) is not int:
        raise ValueError('event threshold must be an integer')
    _validate_path(cert.prefix, n)
    if not isinstance(cert.leaves,tuple) or not cert.leaves:
        raise CertificateError('nonempty tuple of portfolio leaves required')
    leaf_set = set()
    internal = set()
    for p in cert.leaves:
        _validate_path(p,n)
        if p[:len(cert.prefix)] != cert.prefix:
            raise CertificateError('leaf is not a descendant of the root prefix')
        if p in leaf_set:
            raise CertificateError('duplicate leaf')
        leaf_set.add(p)
        for j in range(len(cert.prefix),len(p)):
            internal.add(p[:j])
    if leaf_set & internal:
        raise CertificateError('shadowed leaf/ancestor overlap')

    e = c = Fraction(0)
    count = 0
    visited = set()
    def visit(p: Path) -> None:
        nonlocal e,c,count
        if p in leaf_set:
            visited.add(p)
            mass = extremal_shell_mass(cert.sequence,p,cert.q,
                                       cert.event_at_least,core_length=len(p))
            e += mass.event
            c += mass.complement
            count += mass.distinct_paths
            return
        if p not in internal or len(p)==n:
            raise CertificateError('missing or nonterminal uncovered prefix')
        children = _neighbors(p)
        if not children:
            raise CertificateError('missing dead-end leaf')
        for child in children:
            visit(child)
    visit(cert.prefix)
    if visited != leaf_set:
        raise CertificateError('unreachable leaf')
    return PortfolioMass(e,c,count,len(leaf_set))


def optimize_shell_portfolio(sequence: str, prefix: Path, q: Fraction,
                             event_at_least: int, *, leaf_budget: int,
                             max_core_length: int,
                             objective: Objective = 'event') -> PortfolioPlan:
    """Exact tree-knapsack optimum for the permitted portfolio certificate class.

    At a node, either STOP and claim the exact shell mass of that node,
    with cost one leaf, or EXPAND all SAW children with independently
    optimized disjoint frontier partitions. Convolution of child budget
    tables yields an exact optimum for each feasible number of leaves.
    This optimization is exponential in allowed core depth and polynomial
    in leaf_budget for the explicitly expanded prefix tree; it is not a
    polynomial-time unrestricted HP solver. The returned certificate is
    independently recomputed, not trusted on account of the optimization.
    """
    _check_sequence(sequence)
    _check_q(q)
    if type(event_at_least) is not int:
        raise ValueError('event threshold must be an integer')
    _validate_path(prefix,len(sequence))
    if type(leaf_budget) is not int or leaf_budget < 1:
        raise ValueError('leaf_budget must be a positive integer')
    if type(max_core_length) is not int or not len(prefix) <= max_core_length <= len(sequence):
        raise ValueError('max_core_length must be between prefix and sequence length')
    if objective not in ('event','complement'):
        raise ValueError('objective must be event or complement')

    examined=0
    @lru_cache(maxsize=None)
    def solve(p: Path):
        nonlocal examined
        examined+=1
        mass=extremal_shell_mass(sequence,p,q,event_at_least,core_length=len(p))
        own=mass.event if objective=='event' else mass.complement
        # budget: (score, ordered leaf tuple). All values are exact Fractions.
        best={1:(own,(p,))}
        if len(p)==max_core_length:
            return best
        children=_neighbors(p)
        if not children:
            return best
        products={0:(Fraction(0),())}
        for child in children:
            cands=solve(child)
            updated={}
            for cost0,(score0,paths0) in products.items():
                for cost1,(score1,paths1) in cands.items():
                    cost=cost0+cost1
                    if cost>leaf_budget:
                        continue
                    value=score0+score1
                    record=updated.get(cost)
                    if record is None or value>record[0]:
                        updated[cost]=(value,paths0+paths1)
            products=updated
            if not products:
                break
        for cost, (value, leaves) in products.items():
            record=best.get(cost)
            if record is None or value>record[0]:
                best[cost]=(value,leaves)
        return best

    candidates=solve(prefix)
    cost,(score,paths)=min(candidates.items(),key=lambda item:(-item[1][0],item[0]))
    cert=PortfolioCertificate(sequence,prefix,q,event_at_least,paths)
    verified=verify_shell_portfolio(cert)
    assert verified.leaf_count==cost
    assert (verified.event if objective=='event' else verified.complement)==score
    return PortfolioPlan(cert,verified,objective,max_core_length,leaf_budget,examined)