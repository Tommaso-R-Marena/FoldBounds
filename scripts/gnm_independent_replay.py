#!/usr/bin/env python3
"""Independent exact-integer replay of all GNM residual certificates.

Uses integer common-denominator conversion for binary64 coordinates (NOT
Fraction-based source), direct spring-energy sums (NOT the producer's algebraic
shortcut), and an explicit tree-edge witness. Fail-closed on any mismatch.
No NMR target is read. No external web calls. No floating-point claim checks.
"""
import argparse
from collections import deque
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np

PIN = 'ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'

def sha(b):return hashlib.sha256(b).hexdigest()


def frozen_data(path):
    blob=Path(path).read_bytes()
    if sha(blob)!=PIN:raise ValueError('wrong Stage 1 input')
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        if z.testzip():raise ValueError('corrupt ZIP')
        m=json.loads(z.read('SHA256SUMS.json'))
        if set(m)!=set(z.namelist())-{'SHA256SUMS.json'}:raise ValueError('manifest incomplete')
        for name,value in m.items():
            if sha(z.read(name))!=value:raise ValueError('manifest mismatch '+name)
        raw=z.read('FROZEN_REFERENCE_ONLY.npz')
    with np.load(io.BytesIO(raw),allow_pickle=False) as d:
        return {k:d[k].copy() for k in d.files},sha(raw)


def integer_graph(x):
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():
        raise ValueError('invalid coordinates')
    ratios=[[float(z).as_integer_ratio() for z in row] for row in x]
    K=max(den.bit_length()-1 for row in ratios for _,den in row)
    S=1<<K
    xyz=[[int(num*(S//den)) for num,den in row] for row in ratios]
    edges=[];n=len(x)
    for a in range(n):
        for b in range(a+1,n):
            d2=sum((xyz[a][k]-xyz[b][k])**2 for k in range(3))
            if d2==0:raise ValueError('coincident residues')
            if d2<=100*S*S:edges.append((a,b))
    return edges


def verify_tree(n,edges,root,parent,depth,B):
    es=set(edges)
    if len(parent)!=n or len(depth)!=n or not(0<=root<n):raise ValueError('bad tree shape')
    if parent[root]!=-1 or depth[root]!=0:raise ValueError('bad root')
    for node in range(n):
        if node==root:continue
        p=parent[node]
        if not(type(p)==int and 0<=p<n and
               (min(node,p),max(node,p)) in es and
               depth[node]==depth[p]+1 and 1<=depth[node]<=n-1):
            raise ValueError('tree edge/depth invalid')
    if sum(depth)!=B:raise ValueError('bad Poincare path bound')


def reproduce(stage1_path, cert_path):
    dat,rawsha=frozen_data(stage1_path)
    c=json.loads(Path(cert_path).read_bytes())
    if c.get('format')!='FOLDBOUNDS_GNM_INTEGER_RESIDUAL_V1' or c.get('source_npz_sha256')!=rawsha or c.get('source_archive_sha256')!=PIN:
        raise ValueError('wrong experiment source')
    x=dat['reference_coords']; n=len(x)
    edges=integer_graph(x)
    root=c['anchor']; B=c['sum_tree_depths_B']
    verify_tree(n,edges,root,c['tree_parent'],c['tree_depth'],B)
    if len(edges)!=c['n_contacts'] or c['n_residues']!=n:raise ValueError('contact count mismatch')
    idx={v:k for k,v in enumerate(a for a in range(n) if a!=root)}
    Q=c['denominator_Q'];N=c['integer_inverse_proposal']
    if not(type(Q)==int and 0<Q<=1<<128):raise ValueError('invalid Q')
    if len(N)!=n-1 or any(len(row)!=n-1 for row in N):raise ValueError('bad witness dimensions')
    if any(type(v)!=int or v.bit_length()>256 for row in N for v in row):raise ValueError('bad witness entries')
    if sha(json.dumps(N,separators=(',',':')).encode('ascii'))!=c['integer_witness_sha256']:
        raise ValueError('witness hash mismatch')
    adj=[[] for _ in range(n-1)]
    for i,j in edges:
        if i!=root:adj[idx[i]].append((idx[j] if j!=root else -1))
        if j!=root:adj[idx[j]].append((idx[i] if i!=root else -1))
    pairs=list(zip(map(int,dat['pair_i']),map(int,dat['pair_j'])))
    if len(set(pairs))!=len(pairs):raise ValueError('duplicate score pair')
    bounds=[]
    for i,j in pairs:
        if not(0<=i<j<n):raise ValueError('invalid residue pair')
        ci,cj=idx.get(i),idx.get(j)
        u=[(N[k][ci] if ci is not None else 0) -
           (N[k][cj] if cj is not None else 0) for k in range(n-1)]
        u_with_root=[0 if k==root else u[idx[k]] for k in range(n)]
        energy=sum((u_with_root[a]-u_with_root[b])**2 for a,b in edges)
        bu=u_with_root[i]-u_with_root[j]
        residual=[]
        for row in range(n-1):
            lap_u=sum(u[row]-(u[k] if k!=-1 else 0) for k in adj[row])
            force=int(row==ci)-int(row==cj)
            residual.append(Q*force-lap_u)
        lo=2*Q*bu-energy
        hi=lo+B*sum(v*v for v in residual)
        bounds.append([lo,hi])
    if bounds!=c['score_numerators']:raise ValueError('replay score interval mismatch')
    order=sorted(range(len(bounds)),key=lambda k:(bounds[k][0],k))
    count=sum(bounds[order[k]][1]<bounds[order[k+1]][0] for k in range(len(bounds)-1))
    if count!=c['strict_adjacent_count'] or order!=c['ordered_indices']:
        raise ValueError('replay rank proof mismatch')
    if count==len(pairs)-1:print(f'Independent ALL-PASS: {len(edges)} contacts, {len(pairs)} scores, {count} strict adjacent ranks, B={B}')
    else:print(f'Independent PARTIAL: {count}/{len(pairs)-1} strict orders')
    return count


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage1');ap.add_argument('certificate');a=ap.parse_args()
    reproduce(a.stage1,a.certificate)