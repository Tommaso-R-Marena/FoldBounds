#!/usr/bin/env python3
"""Exact integer Sherman-Morrison update for a grounded contact graph.

Inputs A*N=D*I over Z, D>0. Adding a contact incidence bb^T has
v=N*b, alpha=b^T*v, new D=D*(D+alpha),
new N=(D+alpha)*N-v*v^T. We divide the common integer gcd and verify.
NO floating point or target NMR inputs enter the updater. All computations
are exact. This is an algebraically known identity, not a new force law.
"""
from __future__ import annotations
import math,time
from pathlib import Path
import json


def add_edge_exact(N,D,u,v,anchor,n=124):
    """Add an unweighted edge; u,v are source-residue indices, 0<=u<v<n."""
    if not (type(D)==int and D>0 and 0<=u<v<n and 0<=anchor<n):
        raise ValueError('Invalid exact update input')
    m=n-1
    if len(N)!=m or any(len(row)!=m for row in N):
        raise ValueError('Invalid inverse witness shape')
    iu=(u if u<anchor else u-1) if u!=anchor else None
    iv=(v if v<anchor else v-1) if v!=anchor else None
    if iu is None and iv is None:raise ValueError('Zero contact incidence')
    # v=N*b for b=e_u-e_v, with the anchor coordinate omitted.
    vector=[(N[k][iu] if iu is not None else 0)-
            (N[k][iv] if iv is not None else 0) for k in range(m)]
    alpha=((vector[iu] if iu is not None else 0)
           -(vector[iv] if iv is not None else 0))
    if alpha<=0: raise ValueError('Expected positive effective resistance of a new edge')
    p=D+alpha
    nextD=D*p
    nextN=[[p*N[i][j]-vector[i]*vector[j] for j in range(m)]
           for i in range(m)]
    common=nextD
    for row in nextN:
        for x in row:
            common=math.gcd(common,x)
            if common==1:break
        if common==1:break
    if common>1:
        nextD//=common
        nextN=[[x//common for x in row] for row in nextN]
    return nextN,nextD,{'contact':[u,v],'alpha_numerator_digits':len(str(alpha)),
                       'next_D_digits':len(str(nextD))}


def update_certain_to_nominal(source,exact_certain_certificate,nominal_certificate=None):
    # These source/verification helpers are used only in the PROPOSER.
    from gnm_exact_inverse_producer import archive_source,exact_edges,grounded,check_scaled_inverse,pair_num
    coords,pairs,frozen,npzsha=archive_source(source)
    certain=exact_edges(coords,'certain')
    nominal=exact_edges(coords,'nominal')
    extras=sorted(set(nominal)-set(certain))
    if len(extras)!=48:raise ValueError('Frozen source has changed nominal-on optional contact count')
    c=json.loads(Path(exact_certain_certificate).read_text())
    if c['edges'] != [list(e) for e in certain]:raise ValueError('Wrong certain-graph source')
    N=[[int(x) for x in row] for row in c['N']];D=int(c['D'])
    A,idx=grounded(certain)
    check_scaled_inverse(A,N,D)
    t=time.perf_counter()
    details=[]
    for step,(i,j) in enumerate(extras):
        N,D,detail=add_edge_exact(N,D,i,j,45)
        details.append(detail)
        # Spot-check intermediate certificates; final check always required.
        if step in (0,7,23,47):
            Ai,_=grounded(sorted(set(certain).union(extras[:step+1])))
            check_scaled_inverse(Ai,N,D)
    Af,idx=grounded(nominal)
    check_scaled_inverse(Af,N,D)
    nums=[pair_num(i,j,N,idx) for i,j in pairs]
    if len(set(nums))!=2614:raise ValueError('Expected all distinct nominal pair scores')
    approximate=[z/D for z in nums]
    maxdiff=max(abs(z-f) for z,f in zip(approximate,frozen))
    if maxdiff>1e-10:raise ValueError('Exact update disagrees with frozen original GNM')
    elapsed=time.perf_counter()-t
    if nominal_certificate:
        frozen_nominal=json.loads(Path(nominal_certificate).read_text())
        D2=int(frozen_nominal['D'])
        N2=[[int(v) for v in row] for row in frozen_nominal['N']]
        if any(N[i][j]*D2!=N2[i][j]*D for i in range(123) for j in range(123)):
            raise ValueError('Updated exact inverse differs from independent exact recomputation')
    return {'n_added_contacts':len(extras),'selected_contacts':extras,'elapsed_seconds':elapsed,
            'inverse_denominator_digits':len(str(D)),'max_abs_float_error':maxdiff,
            'all_integer_identities':True,'all_nominal_scores_distinct':True,
            'independently_recomputed_exact_inverse_equal':bool(nominal_certificate)}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('stage1');ap.add_argument('certain');ap.add_argument('nominal',nargs='?')
    a=ap.parse_args();print(json.dumps(update_certain_to_nominal(a.stage1,a.certain,a.nominal),indent=2))


def update_contact_signed(N,D,u,v,anchor,sign,n=124):
    """Exact integer addition/removal of one graph contact (sign=+1/-1).

    No inversion, Cholesky or numerical solver; positive output D is checked.
    """
    if sign not in (-1,1) or type(sign) is not int:
        raise ValueError('ABSTAIN: edge update sign invalid')
    if not (type(D) is int and D>0 and 0<=u<v<n and 0<=anchor<n):
        raise ValueError('ABSTAIN: scaled inverse input invalid')
    m=n-1
    if len(N)!=m or any(len(row)!=m for row in N):
        raise ValueError('ABSTAIN: inverse matrix dimensions invalid')
    a=(u if u<anchor else u-1) if u!=anchor else None
    b=(v if v<anchor else v-1) if v!=anchor else None
    vec=[(N[k][a] if a is not None else 0)-
         (N[k][b] if b is not None else 0) for k in range(m)]
    alpha=(vec[a] if a is not None else 0)-(vec[b] if b is not None else 0)
    coefficient=D+sign*alpha
    if alpha<=0 or coefficient<=0:
        raise ValueError('ABSTAIN: contact update breaks positive definiteness')
    Dnew=D*coefficient
    Nnew=[[coefficient*N[i][j]-sign*vec[i]*vec[j] for j in range(m)]
          for i in range(m)]
    g=abs(Dnew)
    for row in Nnew:
        for x in row:
            g=math.gcd(g,x)
            if g==1:break
        if g==1:break
    if g>1:
        Dnew//=g;Nnew=[[v//g for v in row] for row in Nnew]
    return Nnew,Dnew