#!/usr/bin/env python3
"""Source-only portability of exact coupled rank guarantees to CspB and crambin.

This does NOT read original NMR target evaluation data from the frozen ZIP.
Untrusted SymPy elimination proposes inverse witnesses, checked by all
integer entries of A*N=D*I before every exact rational cone comparison.
"""
from __future__ import annotations
from fractions import Fraction
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix
import hashlib,json,time
from pathlib import Path
from gnm_portability_audit import source, integer_graph, FROZEN, ARCHIVE_SHA


def contact_bounds(coords,epsilon=Fraction(1,20)):
    xyz=[[Fraction(*float(a).as_integer_ratio()) for a in row] for row in coords]
    assert epsilon>=0
    minimum=[];maximum=[];nominal=[]
    for i in range(len(coords)):
        for j in range(i+1,len(coords)):
            deltas=[abs(xyz[i][k]-xyz[j][k]) for k in range(3)]
            dlo=sum(max(d-2*epsilon,0)**2 for d in deltas)
            dhi=sum((d+2*epsilon)**2 for d in deltas)
            actual=sum(d*d for d in deltas)
            if dhi<=100:minimum.append((i,j))
            if dlo<=100:maximum.append((i,j))
            if actual<=100:nominal.append((i,j))
    return minimum,nominal,maximum


def exact_graph_inverse(edges,n,anchor):
    A,ix=integer_graph(edges,n,anchor)
    P,D=DomainMatrix.from_list(A,ZZ).inv_den()
    N=[[int(x) for x in row] for row in P.to_list()]
    D=int(D)
    if D<0:D=-D;N=[[-v for v in row] for row in N]
    if D<=0:raise ValueError('ABSTAIN: nonpositive inverse denominator')
    # independent all-entries integer oracle (not SymPy)
    for i in range(n-1):
        for j in range(n-1):
            if sum(A[i][k]*N[k][j] for k in range(n-1))!=(D if i==j else 0):
                raise ValueError('ABSTAIN: exact inverse does not satisfy integer identity')
    return N,D


def pairnumer(N,p,anchor):
    a,b=p
    x=None if a==anchor else (a if a<anchor else a-1)
    y=None if b==anchor else (b if b<anchor else b-1)
    return ((N[x][x] if x is not None else 0)+
            (N[y][y] if y is not None else 0)-
            (2*N[x][y] if x is not None and y is not None else 0))


def mixed(M,p,q,anchor):
    i,j=p;k,l=q
    terms=[]
    for x,y,s in ((i,k,1),(i,l,-1),(j,k,-1),(j,l,1)):
        if x==anchor or y==anchor:continue
        a=x if x<anchor else x-1;b=y if y<anchor else y-1
        terms.append(s*M[a][b])
    return sum(terms)


def audit(zipfile,protein):
    start=time.perf_counter()
    coords,pairs,frozen=source(zipfile,protein)
    n=len(coords);root=n-1
    certain,nominal,possible=contact_bounds(coords)
    if len(nominal)!=FROZEN[protein][3]:raise ValueError('ABSTAIN: nominal contact graph mismatch')
    Nc,Dc=exact_graph_inverse(certain,n,root)
    Nn,Dn=exact_graph_inverse(nominal,n,root)
    Np,Dp=exact_graph_inverse(possible,n,root)
    scores=[pairnumer(Nn,p,root) for p in pairs]
    order=sorted(range(len(pairs)),key=lambda x:scores[x])
    sorted_pairs=[pairs[k] for k in order]
    rlo=[pairnumer(Np,p,root)*Dc for p in sorted_pairs]
    rhi=[pairnumer(Nc,p,root)*Dp for p in sorted_pairs]
    if any(a>b for a,b in zip(rlo,rhi)):raise ValueError('ABSTAIN: endpoint monotonicity false')
    diff=[[Nc[i][j]*Dp-Np[i][j]*Dc for j in range(n-1)] for i in range(n-1)]
    equal=old=new=old_wrong=0
    examples=[]
    for a in range(len(order)-1):
        for b in range(a+1,len(order)):
            if scores[order[a]]==scores[order[b]]:
                equal+=1;continue
            if rlo[b]>rhi[a]:old+=1;new+=1;continue
            s=(rlo[b]+rhi[b])-(rlo[a]+rhi[a])
            if s<=0:continue
            da=rhi[a]-rlo[a];db=rhi[b]-rlo[b]
            cross=mixed(diff,sorted_pairs[a],sorted_pairs[b],root)
            if da+db+2*cross<=0:continue
            if s*s+4*cross*cross>(da+db)**2:
                new+=1
                if len(examples)<3:examples.append({'rank_indices':[a,b], 'pair_indices':[sorted_pairs[a],sorted_pairs[b]]})
    strict_total=len(pairs)*(len(pairs)-1)//2-equal
    result=dict(protein=protein,source_sha256=ARCHIVE_SHA,
        n_residues=n,n_pairs=len(pairs),n_nominal_contacts=len(nominal),
        n_certain_contacts=len(certain),n_possible_contacts=len(possible),
        n_uncertain_contacts=len(possible)-len(certain),
        n_nominal_exact_tied_comparisons=equal,n_strict_comparisons=strict_total,
        old_independent_certified_order_count=old,
        new_coupled_certified_order_count=new,
        additional_certified_strict_orders=new-old,
        guaranteed_strict_correct_order_fraction_old=str(Fraction(old,strict_total)),
        guaranteed_strict_correct_order_fraction_coupled=str(Fraction(new,strict_total)),
        illustrative_new_witnesses=examples,
        real_protein_integer_inverse_identities_checked=True,
        target_NMR_data_accessed=False,seconds=round(time.perf_counter()-start,3))
    return result

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('frozen_zip');a.add_argument('output_json');s=a.parse_args()
    data=[audit(s.frozen_zip,name) for name in FROZEN]
    Path(s.output_json).write_text(json.dumps({'format':'GNM_COUPLED_CROSSPROTEIN_SOURCE_ONLY_V1','results':data},sort_keys=True,indent=2)+'\n')
    print(json.dumps(data,indent=2))