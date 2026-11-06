#!/usr/bin/env python3
"""Construct an UNTRUSTED source-only Cholesky proposal for an exact SOS certificate.

The proposal never has authority. Independent integer checker must reconstruct
q Q² A - p Q² I - q WᵀW and verify symmetric diagonal dominance. This yields
λ_min(A) >= p/q by an explicit nonnegative sum-of-squares decomposition.
"""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction
import argparse, io, json, zipfile
import numpy as np

STAGE1_SHA='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
VERSION='FOLDBOUNDS_GNM_SOURCE_ONLY_SPECTRAL_SOS_V1'

def read_reference(path):
    raw=Path(path).read_bytes()
    if sha256(raw).hexdigest()!=STAGE1_SHA:raise ValueError('ABSTAIN: stage1 changed')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip():raise ValueError('ABSTAIN: broken reference archive')
        manifest=json.loads(z.read('SHA256SUMS.json'))
        if set(manifest)!=(set(z.namelist())-{'SHA256SUMS.json'}):raise ValueError('ABSTAIN: missing source')
        for k,v in manifest.items():
            if sha256(z.read(k)).hexdigest()!=v:raise ValueError('ABSTAIN: internal SHA mismatch')
        src=z.read('FROZEN_REFERENCE_ONLY.npz')
        with np.load(io.BytesIO(src),allow_pickle=False) as f:
            coords=f['reference_coords'].copy()
        return coords,sha256(src).hexdigest()

def certain_graph(coords):
    if coords.shape!=(124,3) or not np.isfinite(coords).all():
        raise ValueError('ABSTAIN: invalid reference coordinates')
    xyz=[[Fraction.from_float(float(t)) for t in row] for row in coords]
    fixed=[];optional=[];absent=[]
    eps=Fraction(1,20)
    for i in range(len(xyz)):
        for j in range(i+1,len(xyz)):
            d=[abs(xyz[i][k]-xyz[j][k]) for k in range(3)]
            lo=sum(max(Fraction(0), t-2*eps)**2 for t in d)
            hi=sum((t+2*eps)**2 for t in d)
            if lo==0:raise ValueError('ABSTAIN: possible coordinate collision')
            if hi<=100:fixed.append((i,j))
            elif lo>100:absent.append((i,j))
            else:optional.append((i,j))
    if (len(fixed),len(optional),len(absent))!=(972,108,6546):
        raise ValueError('ABSTAIN: source graph classification changed')
    return fixed,optional

def construct(stage1,out,p=169821,q=1000000,bits=30,delta='0.0000003'):
    if not (type(p)==int and type(q)==int and 0<p<q and 6<=bits<=32):
        raise ValueError('ABSTAIN: invalid spectral target')
    coords,npz_sha=read_reference(stage1);edges,optional=certain_graph(coords)
    root=45;n=123;ix={v:i for i,v in enumerate(v for v in range(124) if v!=root)}
    A=np.zeros((n,n),dtype=np.int64)
    for i,j in edges:
        if i!=root:A[ix[i],ix[i]]+=1
        if j!=root:A[ix[j],ix[j]]+=1
        if i!=root and j!=root:A[ix[i],ix[j]]-=1;A[ix[j],ix[i]]-=1
    # Numerical proposal only. Use a *fixed* source-only interior shift,
    # not an inferred bound. The exact certificate below decides acceptance.
    lam=float(Fraction(p,q));chosen_delta=float(Fraction(delta))
    if chosen_delta<=0:raise ValueError('ABSTAIN: δ must be positive')
    try:
        proposed=np.linalg.cholesky(A.astype(float)-(lam+chosen_delta)*np.eye(n)).T
    except np.linalg.LinAlgError as e:
        raise ValueError('ABSTAIN: numerical proposal was not PD') from e
    Q=1<<bits
    if not np.isfinite(proposed).all() or np.max(abs(proposed))*Q>=2**62:
        raise ValueError('ABSTAIN: unsafe proposal')
    W=np.rint(Q*proposed).astype(np.int64)
    if np.any(np.tril(W,k=-1)!=0):raise AssertionError('proposal not upper triangular')
    # CHECK exact integer identity and domination WITHOUT trusting numeric solver.
    # Must force object bigints: q Q² may exceed 64-bit.
    Wobj=np.array(W,dtype=object)
    gram=Wobj.T@Wobj
    S=[[(q*int(A[i,j])-p*int(i==j))*Q*Q - q*int(gram[i,j])
        for j in range(n)] for i in range(n)]
    if any(S[i][j]!=S[j][i] for i in range(n) for j in range(i)):
        raise ValueError('ABSTAIN: broken symmetric residual')
    margin=[S[i][i]-sum(abs(S[i][j]) for j in range(n) if j!=i) for i in range(n)]
    if min(margin)<0:raise ValueError('ABSTAIN: residual not diagonally dominant')
    cert={
      'format':VERSION,'stage1_sha256':STAGE1_SHA,'stage1_npz_sha256':npz_sha,
      'epsilon_A':'0.05','n_residues':124,'root':root,
      'n_certain_edges':len(edges),'n_optional_edges':len(optional),
      'certain_contact_edges':[list(x) for x in edges],
      'lambda_num':p,'lambda_den':q,'cholesky_den_Q':Q,
      'integer_upper_triangular_W':W.tolist(),
      'min_dd_margin_numerator':min(margin),
      'claim':'q*Q²*A - p*Q²*I = q*WᵀW + S. Symmetric diagonally dominant S with S_ii>=sum_j!=i |S_ij| is explicitly a sum of nonnegative diagonal squares and signed pair squares. Hence A >= (p/q) I for every grounded displacement on the fixed contacts, and on any graph with more nonnegative springs.',
      'scope':['Proposes Cholesky from Stage1 reference only; numerical proposal untrusted','All matrix identities and residual dominance require separate exact integer replay','The theorem is about the model, not real protein dynamics','Original blind RNase ANM result remains a FAIL','No 2AAS target, thioredoxin 1XOA target or NMR labels used']
    }
    dest=Path(out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(cert,sort_keys=True,separators=(',',':'))+'\n')
    summary={'status':'SOURCE_ONLY_SPECTRAL_SOS_CERTIFICATE','sha256':sha256(dest.read_bytes()).hexdigest(),
             'p':p,'q':q,'lambda_lower_exact':str(Fraction(p,q)),
             'B_upper_exact':str(Fraction(q,p)),'B_upper_decimal':q/p,'tree_B':262,
             'tightening_vs_tree':262*p/q,'tightening_vs_flow':(72886878439/4294967296)*p/q,
             'min_residual_dd_margin_integer':min(margin),'contact_count':len(edges),
             'matrix_entries':n*n,'witness_size_bytes':dest.stat().st_size}
    print(json.dumps(summary,indent=2,sort_keys=True))
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage1');ap.add_argument('output')
    ap.add_argument('--p',type=int,default=169821);ap.add_argument('--q',type=int,default=1000000)
    ap.add_argument('--bits',type=int,default=30);ap.add_argument('--delta',default='0.0000003')
    a=ap.parse_args();construct(a.stage1,a.output,a.p,a.q,a.bits,a.delta)