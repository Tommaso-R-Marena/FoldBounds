#!/usr/bin/env python3
"""Independent integer replay of FoldBounds source-only EXACT GNM scores.

NO producer import, SymPy import, matrix inverse or floating solver. The complete
scaled-inverse equality A*N = D*I is checked with Python unbounded integers.
The source PDB data must be the byte-for-byte frozen 7RSA Stage1 archive.
2AAS and thioredoxin targets are NEVER read.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile,math
from pathlib import Path
import numpy as np

SOURCE_SHA='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
FORMAT='FOLDBOUNDS_GNM_EXACT_INTEGER_INVERSE_V1'

def loaded_source(path):
    b=Path(path).read_bytes()
    if hashlib.sha256(b).hexdigest()!=SOURCE_SHA:raise ValueError('ABSTAIN: source ZIP hash mismatch')
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        man=json.loads(z.read('SHA256SUMS.json'))
        if z.testzip() or set(man)!=set(z.namelist())-{'SHA256SUMS.json'}:
            raise ValueError('ABSTAIN: invalid frozen manifest')
        if any(hashlib.sha256(z.read(k)).hexdigest()!=v for k,v in man.items()):
            raise ValueError('ABSTAIN: corrupted source ZIP internals')
        raw=z.read('FROZEN_REFERENCE_ONLY.npz')
        with np.load(io.BytesIO(raw),allow_pickle=False) as f:
            coords=f['reference_coords'].copy()
            pi=f['pair_i'].copy();pj=f['pair_j'].copy()
            prediction=f['prediction_GNM_resistance'].copy()
    if coords.shape!=(124,3) or not np.isfinite(coords).all():
        raise ValueError('ABSTAIN: wrong coordinates')
    pairs=list(zip([int(v) for v in pi],[int(v) for v in pj]))
    if len(pairs)!=2614 or len(set(pairs))!=len(pairs) or any(i<0 or j<=i or j>=124 for i,j in pairs):
        raise ValueError('ABSTAIN: invalid frozen pair index')
    return coords,pairs,prediction,hashlib.sha256(raw).hexdigest()

def independent_edges(coords,mode):
    if mode not in ('nominal','certain','possible'):raise ValueError('ABSTAIN: bad mode')
    # Rebuild all dyadic binary64 positions using a single integer scale.
    frac=[[float(v).as_integer_ratio() for v in row] for row in coords]
    denominator=math.lcm(20,*(q for row in frac for _,q in row))
    points=[[p*(denominator//q) for p,q in row] for row in frac]
    two_eps=denominator//10
    limit=100*denominator**2
    edges=[]
    for i in range(124):
        for j in range(i+1,124):
            dif=[abs(points[i][k]-points[j][k]) for k in range(3)]
            if mode=='nominal':
                d2=sum(x*x for x in dif)
                if not d2:raise ValueError('ABSTAIN: coordinate collision')
                included=d2<=limit
            else:
                low=sum(max(0,t-two_eps)**2 for t in dif)
                high=sum((t+two_eps)**2 for t in dif)
                if not low:raise ValueError('ABSTAIN: box permits colliding coordinates')
                included=high<=limit if mode=='certain' else low<=limit
            if included:edges.append((i,j))
    if len(edges)!=({'nominal':1020,'certain':972,'possible':1080}[mode]):
        raise ValueError('ABSTAIN: original graph contact count changed')
    return edges

def strict_bigint(v):
    if not isinstance(v,str) or not v or len(v)>1024 or (v[0]=='-' and len(v)==1):
        raise ValueError('ABSTAIN: integer text invalid')
    s=v[1:] if v.startswith('-') else v
    if not s.isascii() or not s.isdecimal() or (len(s)>1 and s[0]=='0'):
        raise ValueError('ABSTAIN: noncanonical decimal integer')
    return int(v)

def validate_matrix(edges,N,D,anchor):
    n=123
    idx={v:(v if v<anchor else v-1) for v in range(124) if v!=anchor}
    A=[[0]*n for _ in range(n)]
    for i,j in edges:
        if i!=anchor:A[idx[i]][idx[i]]+=1
        if j!=anchor:A[idx[j]][idx[j]]+=1
        if i!=anchor and j!=anchor:
            A[idx[i]][idx[j]]-=1;A[idx[j]][idx[i]]-=1
    if len(N)!=n or any(len(row)!=n for row in N) or D<=0:
        raise ValueError('ABSTAIN: malformed scaled-inverse matrix')
    if any(N[i][j]!=N[j][i] for i in range(n) for j in range(i+1,n)):
        raise ValueError('ABSTAIN: inverse not symmetric')
    # This is the only trusted matrix arithmetic, entirely integer.
    sparse=[[(k,element) for k,element in enumerate(row) if element] for row in A]
    for i in range(n):
        for j in range(n):
            if sum(v*N[k][j] for k,v in sparse[i]) != (D if i==j else 0):
                raise ValueError(f'ABSTAIN: A*N=D*I violated at ({i},{j})')
    return idx

def check(stage1,certificate):
    coords,pairs,pred,npzsha=loaded_source(stage1)
    raw=Path(certificate).read_bytes(); c=json.loads(raw)
    mode=c.get('graph_mode')
    if (c.get('version')!=FORMAT or c.get('stage1_sha256')!=SOURCE_SHA or
        c.get('stage1_npz_sha256')!=npzsha or c.get('n_residues')!=124 or
        c.get('anchor')!=45 or c.get('cutoff_sq_A2')!='100' or
        c.get('epsilon_A')!=('0.05' if mode in ('certain','possible') else None)):
        raise ValueError('ABSTAIN: wrong certificate metadata')
    edges=independent_edges(coords,mode)
    if c.get('edges')!=[list(pair) for pair in edges]:
        raise ValueError('ABSTAIN: incomplete or changed contact network')
    D=strict_bigint(c['D'])
    rawN=c['N']
    if not isinstance(rawN,list) or len(rawN)!=123 or any(not isinstance(row,list) or len(row)!=123 for row in rawN):
        raise ValueError('ABSTAIN: inverse shape invalid')
    N=[[strict_bigint(v) for v in row] for row in rawN]
    idx=validate_matrix(edges,N,D,45)
    nums=[]
    for i,j in pairs:
        a,b=idx.get(i),idx.get(j)
        score=(N[a][a] if a is not None else 0)+(N[b][b] if b is not None else 0)-2*(N[a][b] if a is not None and b is not None else 0)
        if score<=0:raise ValueError('ABSTAIN: nonpositive exact pair resistance')
        nums.append(score)
    if c.get('pair_numerators')!=[str(x) for x in nums]:
        raise ValueError('ABSTAIN: scores do not match verified inverse')
    order=sorted(range(len(nums)),key=lambda k:(nums[k],k))
    if c.get('canonical_rank_indices')!=order:
        raise ValueError('ABSTAIN: incorrect exact rank order')
    distinct=len(set(nums))
    if c.get('n_distinct_pair_scores')!=distinct or c.get('n_exact_ties')!=len(nums)-distinct:
        raise ValueError('ABSTAIN: incorrect equality classification')
    if mode=='nominal':
        approx=np.array([v/D for v in nums])
        err=float(np.max(np.abs(approx-pred)))
        equal=bool(np.array_equal(np.argsort(approx,kind='stable'),np.argsort(pred,kind='stable')))
        if (err>1e-10 or not equal or c.get('nominal_frozen_rank_equal') is not True or
            abs(c['nominal_float_max_absolute_difference']-err)>1e-15):
            raise ValueError('ABSTAIN: exact GNM conflicts with frozen predictions')
    else:
        if c.get('nominal_frozen_rank_equal') is not None or c.get('nominal_float_max_absolute_difference') is not None:
            raise ValueError('ABSTAIN: improper nominal-mode report')
    return {'status':'PASS','mode':mode,'source_sha256':SOURCE_SHA,
            'certificate_sha256':hashlib.sha256(raw).hexdigest(),
            'n_contacts':len(edges),'n_pairs':len(nums),'strict_distinctions':distinct,
            'exact_ties':len(nums)-distinct,'scaled_inverse_identity':True,
            'matrix_size':123,'denominator_digits':len(str(D)),
            'n_approximate_residual_bounds':0}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage1');ap.add_argument('certificate')
    args=ap.parse_args();print(json.dumps(check(args.stage1,args.certificate),indent=2))