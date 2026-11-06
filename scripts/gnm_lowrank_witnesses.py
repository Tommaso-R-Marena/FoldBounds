#!/usr/bin/env python3
"""Untrusted low-rank GNM inverse proposals from common base graphs.

Only used to propose rational witness matrices.  The exact integer residual
checker *always* rechecks each proposal for its particular contact graph.
Sherman-Morrison/Woodbury are classical, not newly discovered here.

For a grounded SPD scalar Laplacian A and k optional added contact incidences
U, the proposal is:
   (A+U U^T)^-1 = A^-1 - A^-1 U (I+U^T A^-1 U)^-1 U^T A^-1.
This form amortizes 2^k expensive base inversions into one base inverse and
small k-by-k updates, and never trusts the resulting floating values.
"""
from __future__ import annotations
import numpy as np

def grounded_incidence(root, nodes, edges):
    if not 0 <= root < nodes:
        raise ValueError('invalid root')
    idx={v:i for i,v in enumerate(w for w in range(nodes) if w!=root)}
    U=np.zeros((nodes-1,len(edges)),dtype=float)
    for col,(a,b) in enumerate(edges):
        if not (0<=a<b<nodes):raise ValueError('invalid directed contact')
        if a!=root:U[idx[a],col]=1
        if b!=root:U[idx[b],col]=-1
    return U

def woodbury_inverse(base_inverse, root, nodes, added_edges):
    X=np.asarray(base_inverse,dtype=float)
    if X.shape!=(nodes-1,nodes-1) or not np.isfinite(X).all():
        raise ValueError('invalid SPD base inverse')
    if len(added_edges)==0:return X.copy()
    U=grounded_incidence(root,nodes,added_edges)
    V=X@U
    T=np.eye(len(added_edges))+U.T@V
    W=np.linalg.solve(T,V.T)
    result=X - V@W
    if not np.isfinite(result).all():raise ValueError('nonfinite corrected inverse')
    return (result+result.T)/2

def integer_proposal(base_inverse,root,nodes,added_edges,Q):
    if type(Q) is not int or Q<=0:raise ValueError('invalid fixed scale')
    inv=woodbury_inverse(base_inverse,root,nodes,added_edges)
    return [[int(round(float(v)*Q)) for v in row] for row in inv]