#!/usr/bin/env python3
"""Reference-only topology/mechanics models; no target labels are accepted.

Weighted generalized ANM (a **known** class of model, Zheng 2008):
  H_a = sum_contacts (7/r)^2 [(1-a)u*u.T + a*I3]
        tensored with [[1,-1],[-1,1]], plus tau I.
At a=0 this exactly recovers the original inverse-square weighted ANM.
At a=1, this is an isotropic distance-weighted *ridge* GNM variant,
NOT the original unweighted pseudoinverse GNM comparator.

A deliberately simple secondary rank-fusion predictor mixes frozen pure
GNM resistance and directional weighted ANM ranks 3:1. This is NOT a novel
physical theory. It is included as a transparency/ablation comparator.
"""
from __future__ import annotations
import hashlib
import json
import numpy as np

CUTOFF_SQ=100.0
RIDGE=0.2
POWER=2
SCALE=7.0
ALPHA=0.5  # fixed reference-only physical interpolation parameter
BLEND_GNM=0.75  # secondary rank fusion only


def contacts(x, cutoff_sq=CUTOFF_SQ):
    x=np.asarray(x,dtype=np.float64)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():
        raise ValueError('bad C-alpha coordinate matrix')
    edges=[]
    for i in range(len(x)):
        for j in range(i+1,len(x)):
            v=x[i]-x[j]; d2=float(v@v)
            if d2<=0:raise ValueError('coincident C-alpha residues')
            if d2<=cutoff_sq:edges.append((i,j,v,d2))
    if len(edges)<len(x)-1:raise ValueError('contact count too low')
    return edges


def rankdata(arr):
    x=np.asarray(arr,float)
    if not np.isfinite(x).all():raise ValueError('invalid score rank input')
    order=np.argsort(x,kind='stable')
    out=np.empty(len(x),dtype=float)
    k=0
    while k<len(order):
        j=k+1
        while j<len(order) and x[order[k]]==x[order[j]]:j+=1
        out[order[k:j]]=0.5*(k+j-1)+1
        k=j
    return out


def scalar_gnm(ref,ii,jj,edges=None):
    ref=np.asarray(ref,float)
    n=len(ref);edges=contacts(ref) if edges is None else edges
    L=np.zeros((n,n),dtype=np.float64)
    for i,j,v,d2 in edges:
        L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
    # Fully specified scalar model; choose anchor using deterministic BFS
    adjacency=[set() for _ in range(n)]
    for i,j,_,_ in edges:adjacency[i].add(j);adjacency[j].add(i)
    dist=np.full((n,n),np.inf)
    for start in range(n):
        dist[start,start]=0;queue=[start]
        for curr in queue:
            for nb in sorted(adjacency[curr]):
                if not np.isfinite(dist[start,nb]):
                    dist[start,nb]=dist[start,curr]+1
                    queue.append(nb)
    if not np.isfinite(dist).all():raise ValueError('GNM network not connected')
    anchor=int(np.argmin(dist.sum(axis=1)))
    ids=[i for i in range(n) if i!=anchor];inv=np.linalg.inv(L[np.ix_(ids,ids)])
    ind=np.full(n,-1,int)
    for m,i in enumerate(ids):ind[i]=m
    scores=[]
    for i,j in zip(ii,jj):
        a,b=int(ind[i]),int(ind[j]);
        q=(inv[a,a] if a!=-1 else 0)+(inv[b,b] if b!=-1 else 0)-2*(inv[a,b] if a!=-1 and b!=-1 else 0)
        scores.append(q)
    return np.asarray(scores,float), {'anchor':anchor,'edges':len(edges)}


def block_anm(ref,ii,jj,alpha,weighted=True,edges=None):
    ref=np.asarray(ref,float)
    if not (0 <= alpha <= 1):raise ValueError('anisotropy mixture must be in [0,1]')
    n=len(ref);edges=contacts(ref) if edges is None else edges
    H=np.zeros((n*3,n*3),dtype=np.float64)
    for i,j,v,d2 in edges:
        u=v/np.sqrt(d2)
        weight=(SCALE**2/d2 if weighted else 1)
        K=weight*((1-alpha)*np.outer(u,u)+alpha*np.eye(3))
        a=slice(3*i,3*i+3);b=slice(3*j,3*j+3)
        H[a,a]+=K;H[b,b]+=K;H[a,b]-=K;H[b,a]-=K
    H.flat[::3*n+1]+=RIDGE
    # Only the proposal uses floating-point; exact GNM certification separate.
    cov=np.linalg.inv(H)
    r=ref[ii]-ref[jj]
    d=np.linalg.norm(r,axis=1)
    if np.any(d<=0):raise ValueError('coincident scored pair')
    u=r/d[:,None]
    blocks=cov.reshape(n,3,n,3)
    C=blocks[ii,:,ii,:]+blocks[jj,:,jj,:]-blocks[ii,:,jj,:]-blocks[jj,:,ii,:]
    score=np.einsum('pi,pij,pj->p',u,C,u)
    if not np.isfinite(score).all():raise ValueError('invalid inverse response')
    return score


def frozen_methods(ref,ii,jj):
    ref=np.asarray(ref,float);ii=np.asarray(ii,int);jj=np.asarray(jj,int)
    edges=contacts(ref)
    gnm,meta=scalar_gnm(ref,ii,jj,edges)
    weighted=block_anm(ref,ii,jj,0,True,edges)
    plain=block_anm(ref,ii,jj,0,False,edges)
    hybrid=block_anm(ref,ii,jj,ALPHA,True,edges)
    deg=np.zeros(len(ref),int)
    for i,j,_,_ in edges:deg[i]+=1;deg[j]+=1
    d=np.linalg.norm(ref[ii]-ref[jj],axis=1)
    fusion=BLEND_GNM*rankdata(gnm)/len(gnm)+(1-BLEND_GNM)*rankdata(weighted)/len(gnm)
    return {
        'GNM_resistance':gnm,
        'ANM_weighted_inverse_square':weighted,
        'ANM_ridge':plain,
        'packing_inverse':1/(deg[ii]+1)+1/(deg[jj]+1),
        'reference_distance':d,
        'sequence_distance':(jj-ii).astype(float),
        'GANM_half_isotropic':hybrid,
        'GNM_ANM_rank_fusion_75_25':fusion,
    },meta