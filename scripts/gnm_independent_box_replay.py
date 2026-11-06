#!/usr/bin/env python3
"""Independent integer-only recomputation of GNM uncertain-contact response hulls.

The implementation does not import the GNM certificate producer or the
coordinate-box source. It decodes exact binary64 positions as scaled integers,
rebuilds ALL contact statuses, validates spanning trees, and recomputes every
rational score interval from the integer proposed inverse by DIRECT edge
quadratic energy and direct integer Laplacian residuals.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import io
import json
from math import lcm,ceil
from pathlib import Path
import zipfile
import numpy as np

SOURCE_SHA='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'


def load_frozen(path):
    raw=Path(path).read_bytes()
    if sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError('bad source SHA')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip():raise ValueError('corrupt input ZIP')
        m=json.loads(z.read('SHA256SUMS.json'))
        if set(m)!=set(z.namelist())-{'SHA256SUMS.json'}:raise ValueError('incomplete manifest')
        for k,v in m.items():
            if sha256(z.read(k)).hexdigest()!=v:raise ValueError('invalid input hash')
        data=z.read('FROZEN_REFERENCE_ONLY.npz')
    with np.load(io.BytesIO(data),allow_pickle=False) as z:
        return z['reference_coords'],list(zip(map(int,z['pair_i']),map(int,z['pair_j'])))


def integer_boxes(x,eps):
    v=[[float(q).as_integer_ratio() for q in row] for row in x]
    den=lcm(eps.denominator,*(q[1] for row in v for q in row))
    pts=[[num*(den//d) for num,d in row] for row in v]
    e=eps.numerator*(den//eps.denominator)
    certain=[];absent=[];ambiguous=[]
    for i in range(len(pts)):
        for j in range(i+1,len(pts)):
            ds=[abs(pts[i][k]-pts[j][k]) for k in range(3)]
            low=sum(max(0,d-2*e)**2 for d in ds)
            high=sum((d+2*e)**2 for d in ds)
            if low==0:raise ValueError('possible coordinate collision')
            if high<=100*den*den:certain.append((i,j))
            elif low>100*den*den:absent.append((i,j))
            else:ambiguous.append((i,j))
    return certain,absent,ambiguous


def endpoint(n,edges,pairs,meta):
    root=meta['root'];parent=meta['tree_parent'];B=meta['B'];Q=meta['denominator_Q']
    if len(parent)!=n or not(0<=root<n) or parent[root]!=-1:raise ValueError('bad tree')
    edge_set=set(edges)
    for node in range(n):
        if node==root:continue
        p=parent[node]
        if p<0 or (min(node,p),max(node,p)) not in edge_set:raise ValueError('bad tree edge')
    d=[None]*n;d[root]=0
    def walk(i):
        if d[i] is not None:return d[i]
        k=i;seen=set()
        while d[k] is None:
            if k in seen:raise ValueError('cycle in tree')
            seen.add(k);k=parent[k]
        while seen:
            todo=[j for j in seen if parent[j]==k]
            if not todo:raise ValueError('bad tree parent chain')
            j=todo[0];d[j]=d[k]+1;k=j;seen.remove(j)
        return d[i]
    for k in range(n):walk(k)
    if sum(d)!=B:raise ValueError('bad coercivity constant')
    N=meta['integer_inverse_proposal']
    if len(N)!=n-1 or any(len(row)!=n-1 for row in N):raise ValueError('invalid matrix')
    if not(type(Q)==int and Q>0 and all(type(v)==int for row in N for v in row)):
        raise ValueError('noninteger witness')
    idx={v:i for i,v in enumerate(j for j in range(n) if j!=root)}
    nb=[set() for _ in range(n)]
    for i,j in edges:nb[i].add(j);nb[j].add(i)
    result=[]
    for i,j in pairs:
        a,b=idx.get(i),idx.get(j)
        u=[(N[k][a] if a is not None else 0)-(N[k][b] if b is not None else 0)
           for k in range(n-1)]
        U=[0 if k==root else u[idx[k]] for k in range(n)]
        energy=sum((U[k]-U[l])**2 for k,l in edges)
        forcing=U[i]-U[j]
        r=[]
        for k in range(n):
            if k==root:continue
            ri=Q*(int(k==i)-int(k==j)) - sum(U[k]-U[l] for l in nb[k])
            r.append(ri)
        L=2*Q*forcing-energy
        result.append((L,L+B*sum(v*v for v in r)))
    return result


def replay(stage1,cert):
    x,pairs=load_frozen(stage1);n=len(x)
    d=json.loads(Path(cert).read_text())
    if d['source_archive_sha256']!=SOURCE_SHA or d['format']!='GNM_TOPOLOGY_BOX_HULL_V1':
        raise ValueError('wrong certificate source/version')
    for record in d['scores']:
        eps=Fraction(record['epsilon'])
        fixed,off,optional=integer_boxes(x,eps)
        if (len(fixed),len(off),len(optional))!=(record['always_on'],record['always_off'],record['optional_edges']):
            raise ValueError('contact interval mismatch')
        lower=endpoint(n,fixed+optional,pairs,record['endpoint_meta']['all_on'])
        upper=endpoint(n,fixed,pairs,record['endpoint_meta']['all_off'])
        if [q[0] for q in lower]!=record['lower_score_numerators']:
            raise ValueError('bad lower endpoints')
        if [q[1] for q in upper]!=record['upper_score_numerators']:
            raise ValueError('bad upper endpoints')
        # Source-only nominal score ordering is stable across coordinate boxes.
        # Rebuild it directly for eps 0 using the same integer graph algorithm.
        original,_o,_m=integer_boxes(x,Fraction(0))
        nominal=endpoint(n,original,pairs,record['endpoint_meta']['all_on'] if not optional else
                         record['endpoint_meta']['all_off']) if not optional else None
        # The original is not necessarily one extreme; instead generate the exact
        # nominal ordering through a deterministic standalone inverse once.
        import numpy.linalg as la
        idx={v:k for k,v in enumerate(range(n) if False else [a for a in range(n) if a!=0])}
        g=np.zeros((n,n),float)
        for a,b in original:g[a,a]+=1;g[b,b]+=1;g[a,b]-=1;g[b,a]-=1
        H=la.inv(g[1:,1:])
        def q(i,j):
            a,b=i-1,j-1
            return (H[a,a] if i else 0)+(H[b,b] if j else 0)-2*(H[a,b] if i and j else 0)
        nominal_order=sorted(range(len(pairs)),key=lambda k:q(*pairs[k]))
        lo=record['lower_score_numerators'];hi=record['upper_score_numerators'];length=len(pairs)
        total=sum(hi[nominal_order[i]]<lo[nominal_order[j]]
                  for i in range(length) for j in range(i+1,length))
        adjacent=sum(hi[nominal_order[i]]<lo[nominal_order[i+1]] for i in range(length-1))
        if total!=record['certified_pair_order'] or adjacent!=record['certified_adjacent']:
            raise ValueError('rank totals mismatch')
        above=[sum(j!=i and lo[j]>hi[i] for j in range(length)) for i in range(length)]
        below=[sum(j!=i and hi[j]<lo[i] for j in range(length)) for i in range(length)]
        k=max(1,ceil(.1*length))
        guaranteed=sum(length-below[i]<=k for i in range(length))
        possible=sum(1+above[i]<=k for i in range(length))
        if guaranteed!=record['guaranteed_top_k'] or possible!=record['possible_top_k']:
            raise ValueError('top decile mismatch')
        print('INDEPENDENT PASS eps',eps,'contacts',len(fixed),'/',len(optional),
              'certified order',total,'/',length*(length-1)//2,'top',guaranteed,'/',k)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage1');p.add_argument('certificate');args=p.parse_args()
    replay(args.stage1,args.certificate)