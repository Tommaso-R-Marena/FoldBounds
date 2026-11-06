#!/usr/bin/env python3
"""Independent exact-integer replay of every adaptive GNM case certificate.

NO imports from gnm_adaptive_partition, gnm_exact_certificates or the topology
producer. Coordinate boxes are rebuilt as scaled integers from frozen binary64
bytes; exact grounded graph residuals are recomputed from every untrusted
integer inverse matrix.  Target ranks are rederived from the pinned 2AAS file.

Verification of the Lean abstract Rayleigh/case-cover theorems is a separate
obligation: this is independent Python byte-and-arithmetic replay, not kernel
execution over the full 124-residue dataset.
"""
from __future__ import annotations
import argparse
from bisect import bisect_left,bisect_right,insort
from fractions import Fraction
from hashlib import sha256
from math import lcm
import json,io,zipfile
from pathlib import Path
import numpy as np

EXPECTED_STAGE1='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
EXPECTED_STAGE2='4f1ea752ada45453a44ea9128b6e7a43f4070e6824062b796a27144c340dfc6c'
VERSION='FOLDBOUNDS_GNM_ADAPTIVE_PARTITION_V1'


def sealed_zip(path, expected):
    raw=Path(path).read_bytes()
    if sha256(raw).hexdigest()!=expected:raise ValueError('ABSTAIN: source archive SHA mismatch')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip():raise ValueError('ABSTAIN: broken ZIP CRC')
        m=json.loads(z.read('SHA256SUMS.json'))
        if set(m)!=set(z.namelist())-{'SHA256SUMS.json'}:
            raise ValueError('ABSTAIN: incomplete archive manifest')
        if any(sha256(z.read(k)).hexdigest()!=v for k,v in m.items()):
            raise ValueError('ABSTAIN: archive content hash mismatch')
        return {name:z.read(name) for name in z.namelist()}


def exact_coord_boxes(x,epsilon):
    eps=Fraction(epsilon)
    if eps<0:raise ValueError('negative epsilon')
    raw=[[float(v).as_integer_ratio() for v in row] for row in x]
    den=lcm(eps.denominator,*(d for row in raw for _,d in row))
    xyz=[[n*(den//d) for n,d in row] for row in raw]
    a=eps.numerator*(den//eps.denominator)
    fixed=[];off=[];optional=[]
    for i in range(len(xyz)):
        for j in range(i+1,len(xyz)):
            distances=[abs(xyz[i][k]-xyz[j][k]) for k in range(3)]
            lower=sum(max(0,d-2*a)**2 for d in distances)
            upper=sum((d+2*a)**2 for d in distances)
            if lower==0:raise ValueError('ABSTAIN: possible coordinate collisions')
            if upper<=100*den*den:fixed.append((i,j))
            elif lower>100*den*den:off.append((i,j))
            else:optional.append((i,j))
    return fixed,off,optional


def tree_checked(n,root,parent,B,edges):
    if type(root)!=int or not 0<=root<n or len(parent)!=n or parent[root]!=-1:
        raise ValueError('bad anchor')
    edges_set=set(edges)
    depth=[None]*n;depth[root]=0
    for i in range(n):
        if i==root:continue
        seen=set();v=i
        while depth[v] is None:
            if v in seen or v==root:raise ValueError('tree cycle')
            seen.add(v)
            w=parent[v]
            if type(w)!=int or not 0<=w<n or (min(v,w),max(v,w)) not in edges_set:
                raise ValueError('invalid tree contact')
            v=w
        while seen:
            children=[u for u in seen if parent[u]==v]
            if len(children)!=1:raise ValueError('invalid tree parent structure')
            u=children[0];depth[u]=depth[v]+1;seen.remove(u);v=u
    if sum(depth)!=B or B<=0:raise ValueError('invalid tree coercivity constant')


def checked_scores(n,edges,pairs,model,root,B,Q):
    if model['n_contacts']!=len(edges):raise ValueError('graph count mismatch')
    N=model['integer_inverse_proposal']
    if len(N)!=n-1 or any(len(row)!=n-1 for row in N):
        raise ValueError('malformed witness matrix')
    if any(type(a)!=int or a.bit_length()>256 for row in N for a in row):
        raise ValueError('noninteger witness')
    if type(Q)!=int or not (1<=Q<=1<<128):raise ValueError('bad witness denominator')
    vertices=[i for i in range(n) if i!=root]
    ix={i:j for j,i in enumerate(vertices)}
    neighbors=[[] for _ in range(n)]
    for i,j in edges:
        neighbors[i].append(j);neighbors[j].append(i)
    # Recompute entire grounded residual matrix R = Q*I - A*N, using
    # adjacency differences, not the producer's grounded_laplacian matrix.
    R=[]
    for row,v in enumerate(vertices):
        rv=[]
        for col in range(n-1):
            Au=sum(N[row][col] - (N[ix[t]][col] if t!=root else 0) for t in neighbors[v])
            rv.append((Q if row==col else 0)-Au)
        R.append(rv)
    # Build all score bounds by independently expanding force, residual, energy.
    result=[]
    for i,j in pairs:
        a=ix.get(i);b=ix.get(j)
        u=[(N[k][a] if a is not None else 0)-(N[k][b] if b is not None else 0)
           for k in range(n-1)]
        rr=[(R[k][a] if a is not None else 0)-(R[k][b] if b is not None else 0)
            for k in range(n-1)]
        forcing=(u[a] if a is not None else 0)-(u[b] if b is not None else 0)
        lo=Q*forcing+sum(v*w for v,w in zip(u,rr))
        hi=lo+B*sum(w*w for w in rr)
        if not lo<=hi:raise AssertionError('bad nonnegative residual gap')
        result.append((lo,hi))
    # Cross-verify energy identity for independently chosen deterministic probe
    # pairs via direct graph-edge sum-of-squares, including grounded node.
    for pidx in sorted({0,len(pairs)//2,len(pairs)-1}):
        i,j=pairs[pidx];a=ix.get(i);b=ix.get(j)
        u=[(N[k][a] if a is not None else 0)-(N[k][b] if b is not None else 0)
           for k in range(n-1)]
        v=lambda k: u[ix[k]] if k!=root else 0
        energy=sum((v(s)-v(t))**2 for s,t in edges)
        forcing=v(i)-v(j)
        if result[pidx][0] != 2*Q*forcing-energy:
            raise ValueError('direct energy disagrees with residual dual')
    return result


def count_independent(lo,hi,observed):
    if len(lo)!=len(hi) or len(lo)!=len(observed):raise ValueError('mismatched vectors')
    n=len(lo)
    if any(type(a)!=int or type(b)!=int or a>b for a,b in zip(lo,hi)):
        raise ValueError('invalid integer interval')
    obs=np.asarray(observed,dtype=np.float64)
    if not np.isfinite(obs).all():raise ValueError('invalid labels')
    indices=sorted(range(n),key=lambda i:float(obs[i]))
    highs=[];lows=[];C=D=U=ties=0;k=0
    while k<n:
        e=k+1
        while e<n and obs[indices[e]]==obs[indices[k]]:e+=1
        for t in range(k,e):
            j=indices[t]
            cc=bisect_left(highs,lo[j]);dd=len(lows)-bisect_right(lows,hi[j])
            C+=cc;D+=dd;U+=len(highs)-cc-dd
        ties+=(e-k)*(e-k-1)//2
        for t in range(k,e):
            j=indices[t];insort(highs,hi[j]);insort(lows,lo[j])
        k=e
    denom=C+D+U
    if denom<1 or denom+ties!=n*(n-1)//2:raise ValueError('count mismatch')
    return {'n_pairs':n,'total_comparisons':n*(n-1)//2,'target_ties':ties,
            'distinct_target_comparisons':denom,'guaranteed_concordant':C,
            'guaranteed_discordant':D,'uncertain':U,
            'lower':str(Fraction(C-D-U,denom)),
            'upper':str(Fraction(C-D+U,denom)),
            'certain_fraction':str(Fraction(C+D,denom))}


def replay(stage1,stage2,certificate):
    a=sealed_zip(stage1,EXPECTED_STAGE1)
    b=sealed_zip(stage2,EXPECTED_STAGE2)
    with np.load(io.BytesIO(a['FROZEN_REFERENCE_ONLY.npz']),allow_pickle=False) as frozen,\
         np.load(io.BytesIO(b['NMR_COORDINATES_AND_FROZEN_PREDICTIONS.npz']),allow_pickle=False) as target:
        coords=frozen['reference_coords'].copy();n=len(coords)
        pairs=list(zip(map(int,frozen['pair_i']),map(int,frozen['pair_j'])))
        if n!=124 or len(pairs)!=2614:raise ValueError('unexpected experiment dimensions')
        if not np.array_equal(coords,target['original_reference_coords']):
            raise ValueError('tampered target reference')
        if not all(np.array_equal(frozen['prediction_'+name],target['frozen_prediction_'+name])
               for name in ['GNM_resistance','ANM_weighted_inverse_square','ANM_ridge',
                            'packing_inverse','reference_distance','sequence_distance']):
            raise ValueError('changed frozen predictions')
        if (not np.array_equal(frozen['pair_i'],target['evaluated_target_i']) or
            not np.array_equal(frozen['pair_j'],target['evaluated_target_j']) or
            not target['eligible_mask'].all() or
            target['target_coords'].shape!=(32,n,3)):
            raise ValueError('target pair / structure mismatch')
        obs=np.var(np.linalg.norm(target['target_coords'][:,target['evaluated_target_i'],:]-
                            target['target_coords'][:,target['evaluated_target_j'],:],axis=-1),
                   axis=0,ddof=1)
        if not np.array_equal(obs,target['observed_distance_variance']):
            raise ValueError('observed target labels changed')
    raw=Path(certificate).read_bytes();cert=json.loads(raw)
    if cert.get('format')!=VERSION or cert.get('stage1_sha256')!=EXPECTED_STAGE1 or cert.get('stage2_sha256')!=EXPECTED_STAGE2:
        raise ValueError('certificate format/source mismatch')
    if cert['epsilon_A']!='0.05' or cert['Q']!=1<<55:
        raise ValueError('unknown epsilon/denominator')
    if cert['n_residues']!=n or cert['n_score_pairs']!=len(pairs):
        raise ValueError('certificate dimension mismatch')
    fixed,off,optional=exact_coord_boxes(coords,cert['epsilon_A'])
    if (len(fixed),len(off),len(optional))!=(cert['always_on'],cert['always_off'],cert['optional_edges']):
        raise ValueError('contact geometry mismatch')
    selected=cert['selected_optional_indices']
    if (not isinstance(selected,list) or len(selected)>6 or
        len(selected)!=len(set(selected)) or
        any(type(j)!=int or not 0<=j<len(optional) for j in selected)):
        raise ValueError('invalid case split edge indices')
    if cert['selected_edges']!=[list(optional[k]) for k in selected]:
        raise ValueError('wrong selected-edge identifiers')
    # Critical: the solver is allowed to propose arbitrary split edges. We
    # verify coverage by enumerating ALL masks, not trust its chosen heuristic.
    root=cert['root'];B=cert['tree_B'];Q=cert['Q']
    tree_checked(n,root,cert['tree_parent'],B,fixed)
    rem=[edge for k,edge in enumerate(optional) if k not in selected]
    num_masks=1<<len(selected)
    cases=cert['cases'];models=cert['models']
    if len(cases)!=num_masks or len(models)!=2*num_masks:
        raise ValueError('missing/duplicate case coverage')
    lb_global=Fraction(1);ub_global=Fraction(-1)
    all_counts=[]
    for mask in range(num_masks):
        rec=cases[mask]
        if (rec['mask']!=mask or rec['all_remaining_on_model']!=2*mask or
            rec['all_remaining_off_model']!=2*mask+1):
            raise ValueError('bad case index/case coverage')
        forced=[optional[selected[k]] for k in range(len(selected)) if mask>>k&1]
        lo_graph=fixed+forced+rem
        hi_graph=fixed+forced
        mlo=models[2*mask];mhi=models[2*mask+1]
        if (mlo['mask']!=mask or mlo['kind']!='all_remaining_on' or
            mhi['mask']!=mask or mhi['kind']!='all_remaining_off'):
            raise ValueError('incorrect endpoint labels')
        lo_result=checked_scores(n,lo_graph,pairs,mlo,root,B,Q)
        hi_result=checked_scores(n,hi_graph,pairs,mhi,root,B,Q)
        lows=[p[0] for p in lo_result];highs=[p[1] for p in hi_result]
        if any(x>y for x,y in zip(lows,highs)):
            raise ValueError('model monotonicity certificate violated')
        verified=count_independent(lows,highs,obs)
        if rec['target_concordance_counts']!=verified:
            raise ValueError('false target-concordance claim in case '+str(mask))
        lb_global=min(lb_global,Fraction(verified['lower']))
        ub_global=max(ub_global,Fraction(verified['upper']))
        all_counts.append(verified)
    if cert['lower_concordance']!=str(lb_global) or cert['upper_concordance']!=str(ub_global):
        raise ValueError('false aggregate result')
    result={'status':'VERIFIED','n_cases':num_masks,'endpoint_graphs':len(models),
            'coverage':'All 2^k selected-edge masks enumerated; remaining optional edges sandwiched by graph monotonicity',
            'lower':str(lb_global),'upper':str(ub_global),
            'positive_guarantee':lb_global>0,'source_sha256':EXPECTED_STAGE1,
            'target_sha256':EXPECTED_STAGE2,'certificate_sha256':sha256(raw).hexdigest()}
    print('INDEPENDENT ADAPTIVE CERTIFICATE VERIFIED',json.dumps(result,sort_keys=True))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage1');p.add_argument('stage2');p.add_argument('certificate')
    args=p.parse_args();replay(args.stage1,args.stage2,args.certificate)