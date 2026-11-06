"""Untrusted rank-one integer Laplacian inverse proposals; exact rechecking.

Starting from verified (nL+11^T)N=dI, the producer proposes a witness
for adding unit contact (u,v) with B'=B+n*b*b^T, b=e_u-e_v:

    N'=(d+n*alpha)*N - n*(N*b)*(N*b)^T
    d'=d*(d+n*alpha), alpha=b^T N b.

GCD reduction shrinks certificates. Every accepted proposal must pass a
SEPARATE full integer multiplication check B'*N'=d'I. This is classical
Sherman-Morrison algebra; no new mathematical novelty is asserted.
"""
from __future__ import annotations
from functools import reduce
from math import gcd

from .gnm_joint_rank import (
    Edge, IntegerInverse, canonical_edges, verify_inverse,
    untrusted_propose_inverse,
)


def propose_add_contact(
    n: int, parent_edges: tuple[Edge, ...],
    parent: IntegerInverse, edge: Edge,
) -> IntegerInverse:
    edges=canonical_edges(n,parent_edges)
    if (not isinstance(edge,tuple) or len(edge)!=2 or
            type(edge[0]) is not int or type(edge[1]) is not int or
            not 0 <= edge[0] < edge[1] < n or edge in edges):
        raise ValueError('edge must be a new non-self contact')
    if not verify_inverse(n,edges,parent):
        raise ValueError('parent inverse witness fails exact check')
    u,v=edge
    N,d=parent.numerator,parent.denominator
    vector=[N[i][u]-N[i][v] for i in range(n)]
    alpha=vector[u]-vector[v]
    if alpha<0:
        raise ArithmeticError('negative effective resistance numerator')
    scale=d+n*alpha
    numerator=[[scale*N[i][j]-n*vector[i]*vector[j] for j in range(n)]
               for i in range(n)]
    den=d*scale
    common=den
    for row in numerator:
        for element in row:
            common=gcd(common,abs(element))
    proposal=IntegerInverse(tuple(tuple(element//common for element in row)
                                  for row in numerator),den//common)
    if not verify_inverse(n,tuple(sorted(set(edges)|{edge})),proposal):
        raise ArithmeticError('rank-one proposal failed complete exact integer identity')
    return proposal


class ExactEndpointCache:
    """Memoize all verified graph inverse witnesses while adding contacts.

    An untrusted proposal can be cached only after exact checking. Holds
    a single reference base graph and refuses targets with missing base edges.
    """
    def __init__(self,n:int,base_edges:tuple[Edge,...],base:IntegerInverse|None=None):
        self.n=n
        self.base_edges=canonical_edges(n,base_edges)
        witness=base if base is not None else untrusted_propose_inverse(n,self.base_edges)
        if not verify_inverse(n,self.base_edges,witness):
            raise ValueError('invalid base inverse')
        self._cache={self.base_edges:witness}

    def get(self,edges:tuple[Edge,...])->IntegerInverse:
        target=canonical_edges(self.n,edges)
        if not set(self.base_edges)<=set(target):
            raise ValueError('endpoint must include all fixed base edges')
        if target not in self._cache:
            extra=sorted(set(target)-set(self.base_edges))
            parent=tuple(sorted(set(target)-{extra[-1]}))
            self._cache[target]=propose_add_contact(
                self.n,parent,self.get(parent),extra[-1])
        return self._cache[target]

    @property
    def count(self)->int:
        return len(self._cache)
