#!/usr/bin/env python3
"""Source-only exact dyadic flow coercivity witness for a GNM contact graph.

Numerical inverse is an UNTRUSTED proposal. Rounded integer edge flows are
repaired by a source-derived spanning tree; final identity is checked with
unbounded integer arithmetic. For every grounded vector u, Cauchy-Schwarz
implies ||u||² <= (sum(P^2)/Q²) * edgeEnergy(u).

The same *single* certain-on flow works for every possible-on supergraph,
regardless of how uncertain contacts are correlated. Never reads NMR labels.
"""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction
import argparse, io, json, zipfile
import numpy as np

STAGE1_SHA='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
VERSION='FOLDBOUNDS_GNM_SOURCE_ONLY_FLOW_V1'
EPSILON=Fraction(1,20)

def frozen_reference(path):
    raw=Path(path).read_bytes()
    if sha256(raw).hexdigest()!=STAGE1_SHA:raise ValueError('ABSTAIN: frozen Stage1 checksum mismatch')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip():raise ValueError('ABSTAIN: broken reference ZIP')
        m=json.loads(z.read('SHA256SUMS.json'))
        if set(m)!=set(z.namelist())-{'SHA256SUMS.json'}:raise ValueError('ABSTAIN: missing/extra file')
        for name,h in m.items():
            if sha256(z.read(name)).hexdigest()!=h:raise ValueError('ABSTAIN: internal digest')
        data=z.read('FROZEN_REFERENCE_ONLY.npz')
        with np.load(io.BytesIO(data),allow_pickle=False) as f:
            x=f['reference_coords'].copy()
        return x,sha256(data).hexdigest()

def classify(x):
    if x.shape!=(124,3) or not np.all(np.isfinite(x)):
        raise ValueError('ABSTAIN: wrong input geometry')
    xyz=[[Fraction.from_float(float(v)) for v in row] for row in x]
    certain=[];optional=[];absent=[]
    for i in range(len(x)):
        for j in range(i+1,len(x)):
            dif=[abs(xyz[i][k]-xyz[j][k]) for k in range(3)]
            lo=sum(max(Fraction(0),d-2*EPSILON)**2 for d in dif)
            hi=sum((d+2*EPSILON)**2 for d in dif)
            if lo==0:raise ValueError('ABSTAIN: potential coordinate collision')
            if hi<=100:certain.append((i,j))
            elif lo>100:absent.append((i,j))
            else:optional.append((i,j))
    if (len(certain),len(optional),len(absent))!=(972,108,6546):
        raise ValueError('ABSTAIN: unexpected contact classification')
    return certain,optional

def tree(n,edges):
    from collections import deque
    nbr=[set() for _ in range(n)]
    for i,j in edges:nbr[i].add(j);nbr[j].add(i)
    def bfs(r):
        par=[-1]*n;dep=[-1]*n;dep[r]=0;q=deque([r])
        while q:
            v=q.popleft()
            for w in sorted(nbr[v]):
                if dep[w]==-1:
                    dep[w]=dep[v]+1;par[w]=v;q.append(w)
        if min(dep)<0:raise ValueError('ABSTAIN: disconnected certain-on graph')
        return par,dep
    root=min(range(n),key=lambda r:(sum(bfs(r)[1]),r))
    par,dep=bfs(root)
    return root,par,sum(dep)

def make_flow(source,output,q_bits=16):
    if type(q_bits)!=int or not 4<=q_bits<=24:raise ValueError('ABSTAIN: bad Q')
    coords,npz_sha=frozen_reference(source)
    n=len(coords); edges,optional=classify(coords)
    root,par,tree_B=tree(n,edges)
    nvar=n-1; E=len(edges); Q=1<<q_bits
    ix={v:k for k,v in enumerate(i for i in range(n) if i!=root)}
    edge_id={e:k for k,e in enumerate(edges)}
    D=np.zeros((E,nvar),dtype=np.int64)
    for k,(a,b) in enumerate(edges):
        if a!=root:D[k,ix[a]]=1
        if b!=root:D[k,ix[b]]=-1
    # Exact integer path inverse R_T with R_T D=I.
    T=np.zeros((nvar,E),dtype=np.int64)
    for v in ix:
        current=v
        while current!=root:
            parent=par[current]; a,b=sorted((current,parent))
            T[ix[v],edge_id[a,b]] += 1 if current==a else -1
            current=parent
    if not np.array_equal(T@D,np.eye(nvar,dtype=np.int64)):
        raise ValueError('ABSTAIN: tree path reconstruction failed')
    A=D.T@D
    # Untrusted electrical-flow proposal, whose error is REPAIRED exactly.
    P_float=np.linalg.solve(A.astype(np.float64),D.T.astype(np.float64))
    if not np.isfinite(P_float).all():raise ValueError('ABSTAIN: proposal invalid')
    W=np.rint(Q*P_float).astype(np.int64)
    residue=Q*np.eye(nvar,dtype=np.int64)-W@D
    W+=residue@T  # exact flow conservation, even for a bad numerical proposal
    if not np.array_equal(W@D,Q*np.eye(nvar,dtype=np.int64)):
        raise ValueError('ABSTAIN: reconstructed integer flows invalid')
    numerator=sum(int(v)*int(v) for v in W.flat)
    if numerator<=0 or numerator>=tree_B*Q*Q:
        raise ValueError('ABSTAIN: no strict source-only improvement')
    certificate={
      'format':VERSION,'stage1_sha256':STAGE1_SHA,'stage1_npz_sha256':npz_sha,
      'epsilon_A':'0.05','n_residues':n,'n_certain_on':len(edges),
      'n_uncertain':len(optional),'anchor':root,'tree_B':tree_B,
      'flow_denominator_Q':Q,'coercivity_numerator':numerator,
      'coercivity_denominator':Q*Q,
      'certain_contact_edges':[list(e) for e in edges],
      'tree_parent':par,'integer_flow_numerators':W.tolist(),
      'mathematical_claim':'A fixed set of source-only integer flows reconstructs each grounded vertex exactly; Cauchy-Schwarz yields ||u||^2 <= B_flow sum_fixed (u_i-u_j)^2 <= B_flow sum_current (u_i-u_j)^2 for ANY allowed uncertain contact network.',
      'limits':['GNM is a model, not ground truth about protein dynamics','Coordinate tolerance is hypothetical, not calibrated','Only source Stage1 7RSA was accessed; 2AAS and 1XOA are not used','Lean general theorem and Python finite witness checking are separate until end-to-end source parsing is formalized']
    }
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(certificate,sort_keys=True,separators=(',',':'))+'\n')
    summary={'status':'SOURCE_ONLY_FLOW_CERTIFICATE','sha256':sha256(output.read_bytes()).hexdigest(),
             'anchor':root,'Q':Q,'certain_edges':len(edges),'optional_edges':len(optional),
             'baseline_B':tree_B,'improved_B_exact':str(Fraction(numerator,Q*Q)),
             'improved_B':numerator/(Q*Q),'tightening_factor':tree_B*Q*Q/numerator,
             'flow_entries':W.size,'bytes':output.stat().st_size}
    print(json.dumps(summary,sort_keys=True,indent=2))
    return summary

if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('stage1');pa.add_argument('output');pa.add_argument('--q-bits',type=int,default=16)
    a=pa.parse_args();make_flow(a.stage1,a.output,a.q_bits)