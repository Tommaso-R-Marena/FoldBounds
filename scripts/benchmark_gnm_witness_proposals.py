#!/usr/bin/env python3
"""Reference-only benchmarking; NOT part of certificate trust base."""
from fractions import Fraction
import time
import argparse
import numpy as np
import gnm_exact_certificates as g
from gnm_topology_boxes import classify_coord_box
from gnm_lowrank_witnesses import integer_proposal

def benchmark(stage1):
    arr,_=g.open_sealed_stage1(stage1)
    coords=arr['reference_coords'];n=len(coords)
    fixed,off,optional=classify_coord_box(coords,Fraction('0.05'))
    graph=[set() for _ in range(n)]
    for i,j in fixed:graph[i].add(j);graph[j].add(i)
    root,B,*_=g.choose_tree(graph)
    # Frozen reference-only selection order; benchmark is independent of target.
    selected=[(51,115),(1,9),(6,32),(43,121),(1,7),(52,76)]
    rem=[e for e in optional if e not in selected]
    base=[fixed+rem,fixed]
    seq=[]
    for mask in range(1<<len(selected)):
        forced=[selected[k] for k in range(len(selected)) if mask>>k&1]
        for j in range(2):
            es=base[j]+forced
            A,_=g.grounded_laplacian(es,n,root)
            seq.append((j,forced,A))
    Q=g.Q_DEFAULT
    start=time.perf_counter()
    direct=[]
    for typ,forced,A in seq:
        N=g.untrusted_inverse_proposal(A,Q)
        # Exact check is mandatory for both policies and is included here.
        g.verify_integer_witness(A,N,Q)
        direct.append(N)
    t_direct=time.perf_counter()-start
    start=time.perf_counter()
    inv_base=[]
    for es in base:
        A,_=g.grounded_laplacian(es,n,root)
        inv_base.append(np.linalg.inv(np.asarray(A,float)))
    updated=[]
    for typ,forced,A in seq:
        N=integer_proposal(inv_base[typ],root,n,forced,Q)
        g.verify_integer_witness(A,N,Q)
        updated.append(N)
    t_woodbury=time.perf_counter()-start
    # Floating proposals can differ by a few integer units; exact certificates
    # remain checked independently, not judged by agreement with direct inverse.
    changes=sum(sum(x!=y for a,b in zip(X,Y) for x,y in zip(a,b))
                for X,Y in zip(direct,updated))
    print({'endpoints':len(seq),'n':n,'direct_seconds':t_direct,
           'lowrank_seconds':t_woodbury,'speedup':t_direct/t_woodbury,
           'witness_entry_differences':changes,
           'both_exact_residual_checked':True})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage1');arg=p.parse_args()
    benchmark(arg.stage1)