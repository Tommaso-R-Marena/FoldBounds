"""Independently checked conformational witnesses for rigorous lower bounds.

A heuristic proposer may submit ANY complete walks. The checker treats them as
untrusted, validates every vertex and prefix, rejects duplicates, and exactly
sums Boltzmann weights of distinct, real self-avoiding completions.

This is a finite 2D HP toy-model guarantee, NOT a Lean-executable proof.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .lattice import (
    CertificateError, DIRECTIONS, Path, _check_q, _check_sequence,
    _contacts, _possible_hh_pairs, _validate_path,
)


@dataclass(frozen=True)
class WitnessCertificate:
    sequence: str
    q: Fraction
    event_at_least: int
    prefix: Path
    completions: tuple[Path, ...]


@dataclass(frozen=True)
class CheckedWitnessMass:
    event: Fraction
    complement: Fraction
    paths: int

    @property
    def partition(self) -> Fraction:
        return self.event + self.complement


def verify_witness_certificate(cert: WitnessCertificate) -> CheckedWitnessMass:
    """Recompute all weights, without trusting any proposer-side scores.

    The complete paths must be unique, valid self-avoiding lattice walks
    extending the exact submitted prefix. Every path contributes once.
    """
    if not isinstance(cert, WitnessCertificate):
        raise CertificateError('not a witness certificate')
    _check_sequence(cert.sequence)
    q=_check_q(cert.q)
    if type(cert.event_at_least) is not int:
        raise CertificateError('event threshold must be an integer')
    n=len(cert.sequence)
    _validate_path(cert.prefix,n)
    if not isinstance(cert.completions,tuple):
        raise CertificateError('completions must be a tuple')
    seen=set()
    pairs=_possible_hh_pairs(cert.sequence)
    event=non_event=Fraction(0)
    for path in cert.completions:
        _validate_path(path,n)
        if len(path)!=n or path[:len(cert.prefix)]!=cert.prefix:
            raise CertificateError('witness does not complete the claimed prefix')
        if path in seen:
            raise CertificateError('repeated complete conformation')
        seen.add(path)
        weight=q**_contacts(path,pairs)
        if _contacts(path,pairs)>=cert.event_at_least:
            event+=weight
        else:
            non_event+=weight
    return CheckedWitnessMass(event,non_event,len(seen))


def propose_dfs_witnesses(prefix: Path, length: int,
                          *, node_budget: int) -> tuple[Path, ...]:
    """Deterministic, *untrusted* depth-first proposal within a node budget.

    One DFS pop costs one unit of budget (including internal nodes).
    Results need independent verification before affecting any claim.
    """
    if type(node_budget) is not int or node_budget<0:
        raise ValueError('node_budget must be a nonnegative integer')
    if type(length) is not int or length<2:
        raise ValueError('length must be an integer >=2')
    _validate_path(prefix,length)
    stack=[prefix]
    result=[]
    visited=0
    while stack and visited<node_budget:
        path=stack.pop()
        visited+=1
        if len(path)==length:
            result.append(path)
            continue
        x,y=path[-1]
        occupied=set(path)
        children=[]
        for dx,dy in DIRECTIONS:
            p=(x+dx,y+dy)
            if p not in occupied:
                children.append(path+(p,))
        stack.extend(reversed(children))
    return tuple(result)
