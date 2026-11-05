#!/usr/bin/env python3
"""Independent exact arithmetic replay of a source-only GNM flow witness.

No import of the producer, no float matrix inverse and no target/NMR file.
Checks all reference bytes and exact 119,556 integer flow entries against
the grounded incidence matrix. Rejected unless flow divergence is *exact*.
"""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction
from math import lcm
import io, json, zipfile
import numpy as np

EXPECTED_STAGE1='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
VERSION='FOLDBOUNDS_GNM_SOURCE_ONLY_FLOW_V1'

def independent_reference(path):
    b=Path(path).read_bytes()
    if sha256(b).hexdigest()!=EXPECTED_STAGE1:raise ValueError('ABSTAIN: Stage1 SHA changed')
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        if z.testzip() is not None:raise ValueError('ABSTAIN: ZIP damaged')
        m=json.loads(z.read('SHA256SUMS.json'))
        if set(m)!=(set(z.namelist())-{'SHA256SUMS.json'}):raise ValueError('ABSTAIN: incomplete manifest')
        for k,v in m.items():
            if sha256(z.read(k)).hexdigest()!=v:raise ValueError('ABSTAIN: content digest mismatch')
        raw=z.read('FROZEN_REFERENCE_ONLY.npz')
        with np.load(io.BytesIO(raw),allow_pickle=False) as f:coords=f['reference_coords'].copy()
        return coords,sha256(raw).hexdigest()

def classify_dyadics(coords):
    if coords.shape!=(124,3) or not np.isfinite(coords).all():
        raise ValueError('ABSTAIN: wrong residue/coordinate shape')
    ratios=[[float(v).as_integer_ratio() for v in r] for r in coords]
    den=lcm(20,*(d for row in ratios for _,d in row))
    xyz=[[p*(den//q) for p,q in row] for row in ratios]
    delta=den//20
    sure=[];possible=[];absent=[]
    for i in range(len(xyz)):
        for j in range(i+1,len(xyz)):
            dif=[abs(xyz[i][k]-xyz[j][k]) for k in range(3)]
            minimum=sum(max(0,d-2*delta)**2 for d in dif)
            maximum=sum((d+2*delta)**2 for d in dif)
            if minimum==0:raise ValueError('ABSTAIN: possible coordinate collision')
            if maximum<=100*den**2:sure.append((i,j))
            elif minimum>100*den**2:absent.append((i,j))
            else:possible.append((i,j))
    return sure,possible,absent

def integer_flow_check(edges,n,anchor,Q,flow):
    if type(anchor)!=int or not 0<=anchor<n or type(Q)!=int or Q<1 or Q>(1<<24):
        raise ValueError('ABSTAIN: invalid flow dimensions/Q')
    m=len(edges);verts=[v for v in range(n) if v!=anchor];idx={v:i for i,v in enumerate(verts)}
    if len(flow)!=len(verts):raise ValueError('ABSTAIN: incomplete flow rows')
    norm=0
    for i,row in enumerate(flow):
        if len(row)!=m:raise ValueError('ABSTAIN: incomplete edge columns')
        divergence=[0]*len(verts)
        for (a,b),x in zip(edges,row):
            if type(x)!=int or x.bit_length()>64:raise ValueError('ABSTAIN: invalid flow entry')
            if a!=anchor:divergence[idx[a]]+=x
            if b!=anchor:divergence[idx[b]]-=x
            norm+=x*x
        for j,value in enumerate(divergence):
            if value != (Q if i==j else 0):
                raise ValueError('ABSTAIN: edge flow conservation identity is false')
    if norm<=0:raise ValueError('ABSTAIN: zero/nonsensical energy')
    return norm

def verify_tree(edges,n,anchor,parent,B):
    if len(parent)!=n or parent[anchor]!=-1:raise ValueError('ABSTAIN: bad tree parent')
    es=set(edges);depths=[0]*n
    for v in range(n):
        if v==anchor:continue
        w=v;seen=set();d=0
        while w!=anchor:
            if w in seen:raise ValueError('ABSTAIN: cycle in provided spanning tree')
            seen.add(w);p=parent[w]
            if type(p)!=int or not 0<=p<n or (min(p,w),max(p,w)) not in es:
                raise ValueError('ABSTAIN: fake spanning tree link')
            w=p;d+=1
        depths[v]=d
    if sum(depths)!=B:raise ValueError('ABSTAIN: wrong tree B')

def verify(stage1,cert_file):
    x,npz_sha=independent_reference(stage1)
    certain,optional,absent=classify_dyadics(x)
    cert_raw=Path(cert_file).read_bytes()
    c=json.loads(cert_raw)
    if c.get('format')!=VERSION or c.get('stage1_sha256')!=EXPECTED_STAGE1 or c.get('stage1_npz_sha256')!=npz_sha:
        raise ValueError('ABSTAIN: unknown certificate/source')
    if (c.get('epsilon_A'),c.get('n_residues'),c.get('n_certain_on'),c.get('n_uncertain'))!=('0.05',124,972,108):
        raise ValueError('ABSTAIN: invalid graph metadata')
    if len(absent)!=6546 or c.get('certain_contact_edges')!=[list(e) for e in certain]:
        raise ValueError('ABSTAIN: graph not exact geometry')
    n=124;root=c['anchor'];Q=c['flow_denominator_Q'];parent=c['tree_parent']
    verify_tree(certain,n,root,parent,c['tree_B'])
    num=integer_flow_check(certain,n,root,Q,c['integer_flow_numerators'])
    if c.get('coercivity_numerator')!=num or c.get('coercivity_denominator')!=Q*Q:
        raise ValueError('ABSTAIN: fabricated B numerator/denominator')
    if not (0<num<c['tree_B']*Q*Q):raise ValueError('ABSTAIN: improvement invalid')
    res={'status':'VERIFIED_SOURCE_ONLY_FLOW','stage1_sha256':EXPECTED_STAGE1,
         'certificate_sha256':sha256(cert_raw).hexdigest(),
         'n_residues':n,'fixed_edges':len(certain),'optional_edges':len(optional),
         'flow_rows':n-1,'flow_entries':(n-1)*len(certain),
         'B_old':c['tree_B'],'B_exact':str(Fraction(num,Q*Q)),
         'B_decimal':num/(Q*Q),'coercivity_improvement':c['tree_B']*Q*Q/num}
    print(json.dumps(res,sort_keys=True,indent=2))
    return res

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('stage1');p.add_argument('certificate')
    a=p.parse_args();verify(a.stage1,a.certificate)