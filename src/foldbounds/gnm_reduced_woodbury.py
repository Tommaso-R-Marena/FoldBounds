"""Exact, independently checked small-matrix Woodbury for GNM contact worlds.

All optional edges are rank-one updates of ONE verified 69x69 integer
inverse. Every proposed kxk inverse, k<=m optional edges, is checked by
stdlib integer multiplication. Accepted contrast decisions use integer
crossproducts, never floating point, square roots, or an unverified inverse.

C_S = (n/d)[N - n Q_S (dI + n E_S^T N E_S)^(-1) Q_S^T] - J/n,
where (nL_certain+J)N=dI and Q=N E. This is a standard Woodbury identity;
it is not a new mathematical theorem or a Lean-kernel byte parser.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import gcd
from typing import Sequence

from .gnm_joint_rank import Edge, IntegerInverse, canonical_edges, verify_inverse, untrusted_propose_inverse


def _small_inverse_proposal(matrix:tuple[tuple[int,...],...]) -> tuple[tuple[tuple[int,...],...],int]:
    k=len(matrix)
    if not k:
        return (),1
    from sympy import ZZ,Matrix
    from sympy.polys.matrices import DomainMatrix
    proposal,den=DomainMatrix.from_Matrix(Matrix(matrix)).convert_to(ZZ).inv_den()
    q=proposal.to_Matrix();delta=int(den)
    if delta<0:
        q=-q;delta=-delta
    numer=[[int(q[i,j]) for j in range(k)] for i in range(k)]
    common=delta
    for row in numer:
        for x in row:common=gcd(common,abs(x))
    return tuple(tuple(x//common for x in row) for row in numer),delta//common


def _verify_small_inverse(matrix:tuple[tuple[int,...],...],
                          numerator:tuple[tuple[int,...],...],denominator:int)->bool:
    k=len(matrix)
    if type(denominator) is not int or denominator<=0 or not isinstance(numerator,tuple) or len(numerator)!=k:
        return False
    if any(not isinstance(row,tuple) or len(row)!=k or any(type(v) is not int for v in row) for row in numerator):
        return False
    if any(numerator[i][j]!=numerator[j][i] for i in range(k) for j in range(i+1,k)):
        return False
    return all(sum(matrix[i][p]*numerator[p][j] for p in range(k))==(denominator if i==j else 0)
               for i in range(k) for j in range(k))


@dataclass
class ReducedCovariance:
    mask:int
    delta:int
    numerator:tuple[tuple[int,...],...]
    selected:tuple[int,...]
    base_num:tuple[tuple[int,...],...]
    edge_columns:tuple[tuple[int,...],...]
    n:int

    def _weighted(self,j:int):
        v=tuple(self.edge_columns[e][j] for e in self.selected)
        return tuple(sum(row[k]*v[k] for k in range(len(v))) for row in self.numerator)

    def entry_numerator(self,i:int,j:int)->int:
        # H_S[i,j] = entry_numerator / delta.
        if not self.selected:
            return self.base_num[i][j]
        weighted=self._weighted(j)
        q=sum(self.edge_columns[e][i]*weighted[k] for k,e in enumerate(self.selected))
        return self.delta*self.base_num[i][j]-self.n*q

    def diagonal_numerators(self)->tuple[int,...]:
        return tuple(self.entry_numerator(i,i) for i in range(self.n))


class ReducedWorldCache:
    """Check base n×n matrix ONCE; exact-check small updates per contact mask."""
    def __init__(self,n:int,certain:tuple[Edge,...],possible:tuple[Edge,...],
                 base_inverse:IntegerInverse|None=None,max_optional:int=20):
        self.n=n
        self.certain=canonical_edges(n,certain)
        self.possible=canonical_edges(n,possible)
        if not set(self.certain)<=set(self.possible):
            raise ValueError('invalid contact envelope')
        self.optional=tuple(sorted(set(self.possible)-set(self.certain)))
        if type(max_optional) is not int or len(self.optional)>max_optional:
            raise ValueError('too many optional edges')
        self.base=base_inverse if base_inverse is not None else untrusted_propose_inverse(n,self.certain)
        if not verify_inverse(n,self.certain,self.base):
            raise ValueError('invalid base integer inverse')
        N=self.base.numerator
        self._cols=tuple(tuple(N[i][u]-N[i][v] for i in range(n)) for u,v in self.optional)
        self._R=tuple(tuple(self._cols[e][u]-self._cols[e][v]
                            for u,v in self.optional) for e in range(len(self.optional)))
        if any(self._R[i][j]!=self._R[j][i] for i in range(len(self._R)) for j in range(len(self._R))):
            raise ArithmeticError('verified base inverse produced an asymmetric resistance matrix')
        self._cache={}

    def get(self,mask:int)->ReducedCovariance:
        m=len(self.optional)
        if type(mask) is not int or not 0<=mask<(1<<m):
            raise ValueError('contact mask out of range')
        if mask not in self._cache:
            sel=tuple(i for i in range(m) if mask&(1<<i))
            d=self.base.denominator
            M=tuple(tuple((d if i==j else 0)+self.n*self._R[sel[i]][sel[j]]
                          for j in range(len(sel))) for i in range(len(sel)))
            T,delta=_small_inverse_proposal(M)
            if not _verify_small_inverse(M,T,delta):
                raise ArithmeticError('UNTRUSTED small inverse failed exact integer check')
            self._cache[mask]=ReducedCovariance(mask,delta,T,sel,
                        self.base.numerator,self._cols,self.n)
        return self._cache[mask]

    @property
    def count(self):return len(self._cache)


def exact_joint_sign(lo:ReducedCovariance,hi:ReducedCovariance,i:int,j:int)->int:
    if lo.n!=hi.n or not 0<=i<j<lo.n:
        raise ValueError('invalid residue contrast')
    dl=lo.delta;dh=hi.delta
    lai=lo.entry_numerator(i,i);laj=lo.entry_numerator(j,j)
    hai=hi.entry_numerator(i,i);haj=hi.entry_numerator(j,j)
    a=lai*dh-hai*dl
    b=laj*dh-haj*dl
    c=lo.entry_numerator(i,j)*dh-hi.entry_numerator(i,j)*dl
    if a<0 or b<0 or a*b<c*c:
        raise ArithmeticError('PSD interval violated for verified endpoints')
    f=(hai-haj)*dl
    s=2*f+a-b
    disc=(a+b)**2-4*c*c
    if disc<0:raise ArithmeticError('invalid exact discriminant')
    return (1 if s>0 and s*s>disc else
           -1 if s<0 and s*s>disc else 0)
