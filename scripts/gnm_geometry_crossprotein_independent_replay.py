#!/usr/bin/env python3
"""Independent exact integer audit for source-only CspB/Crambin no-goods.

No imports from either geometric certificate producer; never opens any
*EVALUATION_DATA* files. Requires pinned prospective archive SHA-256.
"""
import argparse,hashlib,io,itertools,json,math,zipfile
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from gnm_geometry_independent_replay import (
    certificate_contradiction,independent_count,strict_fraction)

SHA='666a5d734fbd61be0f182c954474841bb552b71b3a1e9156d6c21b8c2f095f95'
EXPECTED={'CspB':(67,509,543,561,52,7,2),'Crambin':(46,318,334,354,36,14,1)}

def replay(zip_path,certificate):
    cert=json.loads(Path(certificate).read_text());name=cert.get('protein')
    if name not in EXPECTED or cert.get('format')!='FOLDBOUNDS_FROZEN_CROSSPROTEIN_GEOMETRIC_NOGOOD_V1':
        raise ValueError('unsupported format/protein')
    if hashlib.sha256(Path(zip_path).read_bytes()).hexdigest()!=SHA:
        raise ValueError('source archive mismatch')
    with zipfile.ZipFile(zip_path) as z:
        # Only reference NPZ. Never read NMR target / evaluation members.
        data=z.read(f'{name}_FROZEN_BEFORE_TARGET.npz')
    if cert.get('original_frozen_archive_sha256')!=SHA or cert.get('reference_npz_sha256')!=hashlib.sha256(data).hexdigest():
        raise ValueError('source reference mismatch')
    with np.load(io.BytesIO(data),allow_pickle=False) as a:
        xyz=np.asarray(a['reference_coords'],dtype='float64')
    n,ns,nn,np_,no,k2,k3=EXPECTED[name]
    if xyz.shape!=(n,3) or not np.isfinite(xyz).all():raise ValueError('wrong source shape')
    r=[[float(v).as_integer_ratio() for v in row] for row in xyz]
    D=math.lcm(20,*(q for row in r for _,q in row));points=[[p*(D//q) for p,q in row] for row in r]
    sure=nom=poss=0;optional=[]
    for i,j in itertools.combinations(range(n),2):
        dist=[abs(points[i][k]-points[j][k]) for k in range(3)]
        lo=sum(max(0,t-D//10)**2 for t in dist)
        hi=sum((t+D//10)**2 for t in dist)
        now=sum(t*t for t in dist)
        sure+=hi<=100*D*D;nom+=now<=100*D*D;poss+=lo<=100*D*D
        if lo<=100*D*D and hi>100*D*D:optional.append((i,j))
    if (sure,nom,poss,len(optional))!=(ns,nn,np_,no):raise ValueError('unexpected graph topology')
    if (cert.get('optional_edges')!=[list(x) for x in optional] or
        cert.get('n_residues')!=n or cert.get('contact_counts')!=[ns,nn,np_] or
        cert.get('epsilon_A')!='1/20' or cert.get('cutoff_squared_A')!='100'):
        raise ValueError('contact classification changed')
    clauses=cert.get('forbidden_patterns')
    if not isinstance(clauses,list) or len(clauses)!=k2+k3:raise ValueError('missing clauses')
    records=[];last=None
    for c in clauses:
        ids=c.get('edge_ids');flags=c.get('forbidden_on_bits')
        if (not isinstance(ids,list) or len(ids) not in (2,3) or ids!=sorted(set(ids)) or
            any(type(j)!=int or j<0 or j>=no for j in ids) or
            not isinstance(flags,list) or len(flags)!=len(ids) or
            any(type(s)!=int or s not in (0,1) for s in flags)):
            raise ValueError('malformed clause')
        key=(len(ids),ids,flags)
        if last is not None and key<=last:raise ValueError('duplicate/misordered')
        last=key
        weights=[strict_fraction(t) for t in c.get('lambda',[])]
        gap=certificate_contradiction(points,D,[optional[j] for j in ids],flags,weights)
        if gap<=0 or str(gap)!=c.get('gap'):raise ValueError('false Farkas identity')
        records.append((ids,flags))
    if sum(len(ids)==2 for ids,st in records)!=k2 or sum(len(ids)==3 for ids,st in records)!=k3:
        raise ValueError('wrong clause sizes')
    allowed,components=independent_count(no,records);total=2**no
    if (cert.get('components')!=components or cert.get('abstract_binary_assignments')!=str(total)
        or cert.get('not_excluded_assignments')!=str(allowed) or
        cert.get('excluded_fraction')!=str(F(total-allowed,total)) or
        cert.get('pair_clauses')!=k2 or cert.get('triple_clauses')!=k3):
        raise ValueError('bad exact model count')
    return {'status':'PASS','protein':name,'source_only':True,'pair_nogoods':k2,
            'triple_nogoods':k3,'excluded_percent':float(100*(1-F(allowed,total))),
            'certificate_sha256':hashlib.sha256(Path(certificate).read_bytes()).hexdigest()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('archive');p.add_argument('certificate');x=p.parse_args()
    print(json.dumps(replay(x.archive,x.certificate),indent=2))