#!/usr/bin/env python3
"""Source-only exact GNM inverse certificates (zero floating inverse uncertainty).

Propose an integer matrix N and denominator D with A*N=D*I using SymPy
DomainMatrix.inv_den; then CHECK every identity with independent Python integer
multiplication before serialization. The proposal generator is NOT trusted.
The biological 2AAS target is NEVER accessed. Requires original 7RSA Stage1 ZIP.

python gnm_exact_inverse_producer.py STAGE1.zip out.json --mode nominal
python gnm_exact_inverse_producer.py STAGE1.zip out.json --mode certain
"""
from __future__ import annotations
import argparse, hashlib, io, json, time, zipfile
from pathlib import Path
from fractions import Fraction
import numpy as np
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

SOURCE_SHA = 'ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
VERSION = 'FOLDBOUNDS_GNM_EXACT_INTEGER_INVERSE_V1'
CUTOFF_SQ=Fraction(100)
ANCHOR=45

def archive_source(filename):
    b=Path(filename).read_bytes()
    if hashlib.sha256(b).hexdigest()!=SOURCE_SHA: raise ValueError('ABSTAIN: Stage1 bytes differ from frozen archive')
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        manifest=json.loads(z.read('SHA256SUMS.json'))
        if z.testzip() is not None or set(manifest)!=(set(z.namelist())-{'SHA256SUMS.json'}):
            raise ValueError('ABSTAIN: archive manifest invalid')
        for name,digest in manifest.items():
            if hashlib.sha256(z.read(name)).hexdigest()!=digest:
                raise ValueError('ABSTAIN: source component SHA256 mismatch')
        npz=z.read('FROZEN_REFERENCE_ONLY.npz')
        with np.load(io.BytesIO(npz),allow_pickle=False) as a:
            coords=np.asarray(a['reference_coords'],dtype=np.float64).copy()
            pair_i=a['pair_i'].copy(); pair_j=a['pair_j'].copy()
            original=a['prediction_GNM_resistance'].copy()
    if coords.shape!=(124,3) or not np.isfinite(coords).all() or pair_i.shape!=(2614,) or pair_j.shape!=(2614,):
        raise ValueError('ABSTAIN: invalid original reference geometry/pairs')
    pairs=list(zip([int(v) for v in pair_i],[int(v) for v in pair_j]))
    if len(set(pairs))!=2614 or any(i<0 or i>=j or j>=124 for i,j in pairs):
        raise ValueError('ABSTAIN: invalid reference pairs')
    return coords,pairs,original,hashlib.sha256(npz).hexdigest()

def exact_edges(coords,mode):
    if mode not in ('nominal','certain','possible'): raise ValueError('Invalid contact-graph mode')
    xyz=[[Fraction.from_float(float(v)) for v in row] for row in coords]
    epsilon=Fraction(1,20)
    edges=[]
    for i in range(124):
        for j in range(i+1,124):
            a=[abs(xyz[i][k]-xyz[j][k]) for k in range(3)]
            if mode=='nominal':
                sq=sum(t*t for t in a)
                if not sq: raise ValueError('ABSTAIN: colliding source CA coordinates')
                active=sq<=CUTOFF_SQ
            else:
                low=sum(max(Fraction(0),t-2*epsilon)**2 for t in a)
                high=sum((t+2*epsilon)**2 for t in a)
                if low==0:raise ValueError('ABSTAIN: structural box allows collision')
                active=high<=CUTOFF_SQ if mode=='certain' else low<=CUTOFF_SQ
            if active:edges.append((i,j))
    expected={'nominal':1020,'certain':972,'possible':1080}[mode]
    if len(edges)!=expected:raise ValueError(f'ABSTAIN: expected {expected} edges; found {len(edges)}')
    return edges

def grounded(edges,n=124,anchor=ANCHOR):
    idx={i:(i if i<anchor else i-1) for i in range(n) if i!=anchor}
    A=[[0]*(n-1) for _ in range(n-1)]
    for a,b in edges:
        if a!=anchor:A[idx[a]][idx[a]]+=1
        if b!=anchor:A[idx[b]][idx[b]]+=1
        if a!=anchor and b!=anchor:
            A[idx[a]][idx[b]]-=1;A[idx[b]][idx[a]]-=1
    return A,idx

def check_scaled_inverse(A,N,D):
    n=len(A)
    if type(D)!=int or D<=0 or len(N)!=n or any(len(row)!=n for row in N):
        raise ValueError('ABSTAIN: scaled inverse shape/sign invalid')
    if any(N[i][j]!=N[j][i] for i in range(n) for j in range(i+1,n)):
        raise ValueError('ABSTAIN: inverse not symmetric')
    sparse=[[(k,v) for k,v in enumerate(row) if v] for row in A]
    for i in range(n):
        for j in range(n):
            if sum(v*N[k][j] for k,v in sparse[i]) != (D if i==j else 0):
                raise ValueError('ABSTAIN: integer scaled-inverse identity failed')
    return True

def pair_num(i,j,N,idx):
    x=idx.get(i);y=idx.get(j)
    return (N[x][x] if x is not None else 0)+(N[y][y] if y is not None else 0)-2*(N[x][y] if x is not None and y is not None else 0)

def create(stage1,output,mode):
    t=time.perf_counter()
    coords,pairs,frozen,npzsha=archive_source(stage1)
    edges=exact_edges(coords,mode)
    A,idx=grounded(edges)
    dense=DomainMatrix.from_list(A,ZZ)
    Ni,Di=dense.inv_den()
    N=[[int(v) for v in row] for row in Ni.to_list()]
    D=int(Di)
    if D<0:
        D=-D;N=[[-v for v in row] for row in N]
    check_scaled_inverse(A,N,D)
    nums=[pair_num(i,j,N,idx) for i,j in pairs]
    if any(v<=0 for v in nums):raise ValueError('ABSTAIN: nonpositive effective resistance')
    order=sorted(range(len(nums)),key=lambda i:(nums[i],i))
    distinct=len(set(nums)); ties=len(nums)-distinct
    maxdiff=None
    rank_equal=None
    if mode=='nominal':
        predicted=np.array([v/D for v in nums])
        maxdiff=float(np.max(np.abs(predicted-frozen)))
        rank_equal=bool(np.array_equal(np.argsort(predicted,kind='stable'), np.argsort(frozen,kind='stable')))
        if maxdiff>1e-10 or not rank_equal:
            raise ValueError('ABSTAIN: nominal exact score disagrees with frozen GNM reference')
    record={
      'version':VERSION,'stage1_sha256':SOURCE_SHA,'stage1_npz_sha256':npzsha,
      'graph_mode':mode,'n_residues':124,'anchor':ANCHOR,'cutoff_sq_A2':'100',
      'epsilon_A':'0.05' if mode in ('certain','possible') else None,
      'edges':[list(pair) for pair in edges],
      'D':str(D),'N':[[str(v) for v in row] for row in N],
      'pair_numerators':[str(v) for v in nums],
      'canonical_rank_indices':order,'n_distinct_pair_scores':distinct,'n_exact_ties':ties,
      'nominal_float_max_absolute_difference':maxdiff,
      'nominal_frozen_rank_equal':rank_equal,
    }
    raw=(json.dumps(record,sort_keys=True,separators=(',',':'))+'\n').encode()
    Path(output).write_bytes(raw)
    report={'mode':mode,'edges':len(edges),'pairs':len(nums),'distinct':distinct,'ties':ties,
            'D_decimal_digits':len(str(D)),'certificate_bytes':len(raw),
            'certificate_sha256':hashlib.sha256(raw).hexdigest(),
            'seconds':round(time.perf_counter()-t,3),'nominal_max_abs_float_error':maxdiff,
            'all_integer_identities_checked':True}
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage1');ap.add_argument('output');ap.add_argument('--mode',choices=('nominal','certain','possible'),required=True)
    args=ap.parse_args()
    print(json.dumps(create(args.stage1,args.output,args.mode),indent=2))