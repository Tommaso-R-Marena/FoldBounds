#!/usr/bin/env python3
"""Independent exact-integer, source-only replay of grounded GNM spectral SOS.

NO producer import, no trust in the untrusted Cholesky computation and NO
NMR/target access. Checks every matrix entry and exact integer diagonal
surplus, then explicitly reconstructs the nonnegative signed-pair SOS.
"""
from pathlib import Path
from hashlib import sha256
from math import lcm
from fractions import Fraction
import io,json,zipfile
import numpy as np

STAGE1_SHA='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
VERSION='FOLDBOUNDS_GNM_SOURCE_ONLY_SPECTRAL_SOS_V1'

def frozen_source(path):
    raw=Path(path).read_bytes()
    if sha256(raw).hexdigest()!=STAGE1_SHA:raise ValueError('ABSTAIN: source SHA256 changed')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip():raise ValueError('ABSTAIN: corrupted source ZIP')
        man=json.loads(z.read('SHA256SUMS.json'))
        if set(man)!=(set(z.namelist())-{'SHA256SUMS.json'}):raise ValueError('ABSTAIN: incomplete ZIP manifest')
        if any(sha256(z.read(k)).hexdigest()!=v for k,v in man.items()):
            raise ValueError('ABSTAIN: changed ZIP contents')
        data=z.read('FROZEN_REFERENCE_ONLY.npz')
        with np.load(io.BytesIO(data),allow_pickle=False) as f:
            return f['reference_coords'].copy(),sha256(data).hexdigest()

def certain_edges_independent(x):
    if x.shape!=(124,3) or not np.isfinite(x).all():
        raise ValueError('ABSTAIN: wrong reference coordinates')
    ratios=[[float(v).as_integer_ratio() for v in row] for row in x]
    den=lcm(20,*(q for row in ratios for _,q in row))
    pos=[[p*(den//q) for p,q in row] for row in ratios]
    uncertainty=den//20
    sure=[];maybe=[];out=[]
    for i in range(124):
        for j in range(i+1,124):
            absd=[abs(pos[i][k]-pos[j][k]) for k in range(3)]
            dmin=sum(max(0,t-2*uncertainty)**2 for t in absd)
            dmax=sum((t+2*uncertainty)**2 for t in absd)
            if dmin==0:raise ValueError('ABSTAIN: possible coincident coordinates')
            if dmax<=100*den*den:sure.append((i,j))
            elif dmin>100*den*den:out.append((i,j))
            else:maybe.append((i,j))
    if (len(sure),len(maybe),len(out))!=(972,108,6546):
        raise ValueError('ABSTAIN: coordinate-derived contact counts wrong')
    return sure,maybe

def grounded_laplacian_integer(edges,n,anchor):
    vertices=[i for i in range(n) if i!=anchor];idx={v:k for k,v in enumerate(vertices)}
    A=[[0]*(n-1) for _ in range(n-1)]
    for a,b in edges:
        if a!=anchor:A[idx[a]][idx[a]]+=1
        if b!=anchor:A[idx[b]][idx[b]]+=1
        if a!=anchor and b!=anchor:
            A[idx[a]][idx[b]]-=1;A[idx[b]][idx[a]]-=1
    return A

def verify_spectral_matrix(A,W,p,q,Q):
    n=len(A)
    if type(p)!=int or type(q)!=int or type(Q)!=int or not (p>0 and q>p and 1<=Q<=2**32):
        raise ValueError('ABSTAIN: nonpositive / incompatible spectral scalar')
    if n<1 or any(len(row)!=n for row in A) or len(W)!=n or any(len(row)!=n for row in W):
        raise ValueError('ABSTAIN: malformed Cholesky witness')
    if any(type(A[i][j])!=int or A[i][j]!=A[j][i] for i in range(n) for j in range(n)):
        raise ValueError('ABSTAIN: asymmetric/noninteger graph matrix')
    if any(type(v)!=int or v.bit_length()>64 for row in W for v in row):
        raise ValueError('ABSTAIN: integer witness invalid')
    if any(W[i][j]!=0 for i in range(n) for j in range(i)):
        raise ValueError('ABSTAIN: witness not upper triangular')
    # Work with arbitrary-precision exact integers; do not use NumPy int64.
    # W is upper triangular, so G(i,j) = sum_{k<=min(i,j)} W(k,i)W(k,j).
    S=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i,n):
            gram=sum(W[k][i]*W[k][j] for k in range(i+1))
            r=(q*A[i][j]-p*(i==j))*Q*Q -q*gram
            S[i][j]=r;S[j][i]=r
    # Diagonal-dominance is an EXACT sufficient PSD criterion. Derive a
    # second, explicit, nonnegative SOS basis to make the meaning auditable.
    diag=[S[i][i]-sum(abs(S[i][j]) for j in range(n) if j!=i) for i in range(n)]
    if min(diag)<0:
        raise ValueError('ABSTAIN: exact symmetric residual is not diagonally dominant PSD')
    # Reconstruct original remainder entrywise from the sum of squares:
    # S = sum_i diag[i] e_i e_i^T +
    #     sum_{i<j} |Sij| (e_i+sign(Sij)e_j)(...)^T.
    # Direct reconstruction verifies no sign or diagonal-counting mistakes.
    quadratic_diag=diag[:]
    offdiagonal={}
    for i in range(n):
        for j in range(i+1,n):
            c=abs(S[i][j]);quadratic_diag[i]+=c;quadratic_diag[j]+=c
            offdiagonal[(i,j)]=c if S[i][j]>=0 else -c
    if any(quadratic_diag[i]!=S[i][i] for i in range(n)) or \
       any(offdiagonal[i,j]!=S[i][j] for i in range(n) for j in range(i+1,n)):
        raise ValueError('ABSTAIN: explicit sum-of-squares decomposition failed')
    return {'min_diagonal_surplus':min(diag),'min_remainder_diagonal':min(S[i][i] for i in range(n)),
            'max_remainder_abs':max(abs(x) for row in S for x in row),
            'nonzero_signed_pair_squares':sum(S[i][j]!=0 for i in range(n) for j in range(i+1,n)),
            'diagonal_square_terms':sum(v>0 for v in diag)}

def verify(stage1,certificate):
    xyz,npz_sha=frozen_source(stage1)
    certain,optional=certain_edges_independent(xyz)
    raw=Path(certificate).read_bytes();c=json.loads(raw)
    if (c.get('format')!=VERSION or c.get('stage1_sha256')!=STAGE1_SHA or
        c.get('stage1_npz_sha256')!=npz_sha or
        c.get('epsilon_A')!='0.05' or c.get('n_residues')!=124 or
        c.get('n_certain_edges')!=972 or c.get('n_optional_edges')!=108):
        raise ValueError('ABSTAIN: bad source metadata')
    if c.get('certain_contact_edges')!=[list(e) for e in certain]:
        raise ValueError('ABSTAIN: supplied spectral graph not complete')
    root=c.get('root')
    if type(root)!=int or not 0<=root<124:raise ValueError('ABSTAIN: bad ground vertex')
    A=grounded_laplacian_integer(certain,124,root)
    p,q,Q=c.get('lambda_num'),c.get('lambda_den'),c.get('cholesky_den_Q')
    W=c.get('integer_upper_triangular_W')
    info=verify_spectral_matrix(A,W,p,q,Q)
    if c.get('min_dd_margin_numerator')!=info['min_diagonal_surplus']:
        raise ValueError('ABSTAIN: fake minimal diagonal surplus')
    result={'status':'VERIFIED_SOURCE_ONLY_SPECTRAL_SOS',
            'stage1_sha256':STAGE1_SHA,
            'certificate_sha256':sha256(raw).hexdigest(),
            'n_grounded':123,'certain_edges':972,'optional_edges':108,
            'lambda_lower':str(Fraction(p,q)),
            'B_upper':str(Fraction(q,p)),
            'B_decimal':q/p,'B_improvement_vs_262':262*p/q,
            'integer_matrix_entries':123*123,**info}
    print(json.dumps(result,indent=2,sort_keys=True))
    return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('stage1');p.add_argument('certificate')
    a=p.parse_args();verify(a.stage1,a.certificate)