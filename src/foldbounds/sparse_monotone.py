"""Sparse prefix-free monotone-family lower certificates with an exact tree knapsack.

All masses are for the anchored 2D HP lattice, rational q>=1. The optimizer is
untrusted; the checker recomputes monotone-family masses. This Python program
is not Lean-verified, and no complete-conformation coverage is required.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Literal
from .lattice import (Path, CertificateError, _check_q, _check_sequence,
                      _neighbors, _validate_path)
from .monotone import monotone_family_mass

Objective = Literal['event', 'complement']

@dataclass(frozen=True)
class SparseMonotoneCertificate:
    sequence: str
    prefix: Path
    q: Fraction
    event_at_least: int
    leaves: tuple[Path, ...]

@dataclass(frozen=True)
class SparseMonotoneMass:
    event: Fraction
    complement: Fraction
    distinct_paths: int
    leaf_count: int
    @property
    def partition(self) -> Fraction:
        return self.event + self.complement

@dataclass(frozen=True)
class SparseMonotonePlan:
    certificate: SparseMonotoneCertificate
    mass: SparseMonotoneMass
    objective: Objective
    max_core_length: int
    leaf_budget: int
    candidate_prefixes: int

def verify_sparse_monotone(cert: SparseMonotoneCertificate) -> SparseMonotoneMass:
    """Check sparse prefix-free leaves and independently compute lower masses.

    Each monotone family contains only valid actual continuations; every
    complete path has a unique selected prefix because leaves are prefix-free.
    Skipped subtrees contribute zero and never invalidate a lower bound.
    """
    if not isinstance(cert, SparseMonotoneCertificate):
        raise CertificateError('sparse monotone certificate required')
    _check_sequence(cert.sequence)
    _check_q(cert.q)
    n=len(cert.sequence)
    if type(cert.event_at_least) is not int:
        raise ValueError('event threshold must be integer')
    _validate_path(cert.prefix,n)
    if type(cert.leaves) is not tuple:
        raise CertificateError('leaves must be a tuple')
    seen=set()
    event=complement=Fraction(0)
    count=0
    for leaf in cert.leaves:
        _validate_path(leaf,n)
        if leaf[:len(cert.prefix)]!=cert.prefix:
            raise CertificateError('leaf does not extend root prefix')
        if leaf in seen:
            raise CertificateError('duplicate leaf')
        if any(leaf[:i] in seen for i in range(len(cert.prefix),len(leaf))):
            raise CertificateError('selected ancestor shadows descendant')
        if any(old[:len(leaf)]==leaf for old in seen):
            raise CertificateError('selected descendant shadows ancestor')
        seen.add(leaf)
        mass=monotone_family_mass(cert.sequence,leaf,cert.q,cert.event_at_least)
        event+=mass.event
        complement+=mass.complement
        count+=mass.distinct_paths
    return SparseMonotoneMass(event,complement,count,len(seen))

def optimize_sparse_monotone(sequence: str,prefix: Path,q: Fraction,
                            event_at_least: int,*,leaf_budget: int,
                            max_core_length: int,
                            objective: Objective='event') -> SparseMonotonePlan:
    """Exact budget optimum within sparse monotone antichain family.

    Each prefix: SKIP (cost zero), STOP (cost one), or recursively EXPAND
    disjoint children. Convolve exact rational child budget tables. Work is
    exponential in maximum core depth; budget counts selected leaves, not CPU.
    """
    _check_sequence(sequence)
    _check_q(q)
    if type(event_at_least) is not int:
        raise ValueError('event threshold must be integer')
    _validate_path(prefix,len(sequence))
    if type(leaf_budget) is not int or leaf_budget<1:
        raise ValueError('leaf_budget must be positive integer')
    if type(max_core_length) is not int or not len(prefix)<=max_core_length<=len(sequence):
        raise ValueError('invalid maximum core length')
    if objective not in ('event','complement'):
        raise ValueError('objective must be event or complement')
    visited=0
    @lru_cache(maxsize=None)
    def solve(path: Path):
        nonlocal visited
        visited+=1
        family=monotone_family_mass(sequence,path,q,event_at_least)
        mass=family.event if objective=='event' else family.complement
        best={0:(Fraction(0),()),1:(mass,(path,))}
        if len(path)==max_core_length:
            return best
        children=_neighbors(path)
        if not children:
            return best
        products={0:(Fraction(0),())}
        for child in children:
            options=solve(child)
            updated={}
            for cost0,(value0,paths0) in products.items():
                for cost1,(value1,paths1) in options.items():
                    cost=cost0+cost1
                    if cost>leaf_budget:
                        continue
                    value=value0+value1
                    prior=updated.get(cost)
                    if prior is None or value>prior[0]:
                        updated[cost]=(value,paths0+paths1)
            products=updated
        for cost,(value,paths) in products.items():
            prior=best.get(cost)
            if prior is None or value>prior[0]:
                best[cost]=(value,paths)
        return best
    choices=solve(prefix)
    cost,(score,paths)=min(choices.items(),key=lambda item:(-item[1][0],item[0]))
    cert=SparseMonotoneCertificate(sequence,prefix,q,event_at_least,paths)
    checked=verify_sparse_monotone(cert)
    if checked.leaf_count!=cost or (checked.event if objective=='event' else checked.complement)!=score:
        raise AssertionError('proposed optimum disagrees with independent checker')
    return SparseMonotonePlan(cert,checked,objective,max_core_length,leaf_budget,visited)
