#!/usr/bin/env python3
"""Untrusted SOURCE-ONLY discovery of rational coupled-GNM rank witnesses.

Do not use the candidate floating filter as a proof: EVERY emitted pair is
checked again with arbitrary-precision integer arithmetic before writing.
No NMR data, Stage2, external network calls, or prospective targets.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from gnm_exact_case_independent_replay import (
    read_checked_base, independent_update, loaded_source, independent_edges,
    SELECTED, TOTAL, val,
)
from gnm_exact_case_source import count_nominal_pairwise_agreement
from fractions import Fraction

VERSION='FOLDBOUNDS_SOURCE_ONLY_COUPLED_RANK_V1'
ROOT=45

def scaled_pair_num(N,i,j,root=ROOT):
    a=None if i==root else (i if i<root else i-1)
    b=None if j==root else (j if j<root else j-1)
    v=(N[a][a] if a is not None else 0)+(N[b][b] if b is not None else 0)
    if a is not None and b is not None:v-=2*N[a][b]
    return v

def force_entries(pair,root=ROOT):
    i,j=pair
    return tuple((x if x<root else x-1,sgn) for x,sgn in ((i,1),(j,-1)) if x!=root)

def mixed_cross(M,f,g):
    return sum(si*sj*M[i][j] for i,si in f for j,sj in g)

def exact_coupled_sign(lo,hi,offdiag,f_i,f_j,ri,rj):
    """True iff the *exact rational cone bound* certifies R_j>R_i.

    Here all Rhi/Rlo have common positive denominator Dhi*Dlo and
    offdiag = Dhi*Nlo-Dlo*Nhi, encoding K=Alo^-1-Ahi^-1 >=0.
    Certificate: s>0 and s^2+4c^2>(di+dj)^2.
    """
    s=(lo[rj]+hi[rj])-(lo[ri]+hi[ri])
    if s<=0:return False
    drop_i=hi[ri]-lo[ri];drop_j=hi[rj]-lo[rj]
    if drop_i<0 or drop_j<0:raise ValueError('invalid endpoint monotonicity')
    c=mixed_cross(offdiag,f_i,f_j)
    # Needed for the Lean discriminant theorem; if the norm vanishes abstain.
    if drop_i+drop_j+2*c<=0:return False
    return s*s+4*c*c>(drop_i+drop_j)**2


def matrix_float(N,D):
    return np.array([[float(n/D) for n in row] for row in N],dtype=np.float64)


def generate(stage1,certain_path,possible_path,nominal_path,case_report_path,out):
    start=time.perf_counter()
    coords,pairs,_,npzsha=loaded_source(stage1)
    certain=set(independent_edges(coords,'certain'))
    possible=set(independent_edges(coords,'possible'))
    if len(certain)!=972 or len(possible)!=1080 or len(pairs)!=2614:raise ValueError('source integrity')
    baseline=json.loads(Path(case_report_path).read_text())
    if baseline['format']!='FOLDBOUNDS_EXACT_GNM_SOURCE_ONLY_64_CASE_RANK_STABILITY_V1':raise ValueError('wrong baseline')
    certfiles={'certain':certain_path,'possible':possible_path,'nominal':nominal_path}
    for mode,file in certfiles.items():
        if baseline[mode+'_certificate_sha256']!=hashlib.sha256(Path(file).read_bytes()).hexdigest():
            raise ValueError('baseline attachment hash mismatch')
    sure={0:read_checked_base(stage1,certain_path,'certain')}
    possible_state={0:read_checked_base(stage1,possible_path,'possible')}
    nominal_N,nominal_D=read_checked_base(stage1,nominal_path,'nominal')
    nominal=[scaled_pair_num(nominal_N,i,j) for i,j in pairs]
    if len(set(nominal))!=2614:raise ValueError('nominal ties')
    order=sorted(range(len(pairs)),key=lambda k:nominal[k])
    forces=[force_entries(pairs[k]) for k in order]
    for mask in range(1,64):
        step=(mask&-mask).bit_length()-1; prev=mask^(1<<step)
        sure[mask]=independent_update(*sure[prev],*SELECTED[step],1)
        possible_state[mask]=independent_update(*possible_state[prev],*SELECTED[step],-1)
    records=[]
    total_extra=0
    for mask in range(64):
        Nlo,Dlo=sure[mask];Nhi,Dhi=possible_state[(~mask)&63]
        # integer rational endpoints: Rlo<=Rtrue<=Rhi
        lows=[scaled_pair_num(Nhi,*pairs[k])*Dlo for k in order]
        highs=[scaled_pair_num(Nlo,*pairs[k])*Dhi for k in order]
        if any(x>y for x,y in zip(lows,highs)):raise ValueError('invalid rational endpoint order')
        M=[[Nlo[i][j]*Dhi-Nhi[i][j]*Dlo for j in range(123)] for i in range(123)]
        # only discovery uses floating point, followed by strict integer check
        a=matrix_float(Nlo,Dlo);b=matrix_float(Nhi,Dhi)
        F=np.zeros((2614,123),dtype=float)
        for k,f in enumerate(forces):
            for x,sgn in f:F[k,x]=sgn
        lo_float=np.einsum('ki,ij,kj->k',F,b,F)
        hi_float=np.einsum('ki,ij,kj->k',F,a,F)
        K=a-b;cross_basis=F@K
        delta=hi_float-lo_float; center=(hi_float+lo_float)/2
        extra=[]
        for i in range(2613):
            cross=cross_basis[i]@F[i+1:].T
            displ=center[i+1:]-center[i]
            radsq=np.maximum((delta[i]+delta[i+1:])**2-4*cross**2,0)
            gap=4*displ**2-radsq
            # Conservative numerical candidate selection; exact conditions
            # are the only accepted proof of each reported witness.
            margin=1e-8*np.maximum(4*displ**2+radsq,1e-25)
            candidate=np.nonzero((displ>0)&(gap>margin)&(lo_float[i+1:]<=hi_float[i]+1e-11))[0]
            for j0 in candidate:
                j=i+int(j0)+1
                # Already guaranteed by the old interval method -> skip;
                # old and new witnesses must be disjoint.
                if lows[j]>highs[i]:continue
                if exact_coupled_sign(lows,highs,M,forces[i],forces[j],i,j):
                    extra.append([i,j])
        old=baseline['cases'][mask]['definitely_correct']
        records.append({'mask':mask,'extra_rank_pairs':extra,'extra_count':len(extra),
                        'baseline_correct_count':old,
                        'certified_lower_numerator':2*(old+len(extra))-TOTAL})
        total_extra+=len(extra)
        if mask%8==0:print('case',mask,'extras',len(extra),'lower',records[-1]['certified_lower_numerator'],'sec',round(time.perf_counter()-start,2),flush=True)
    minimum=min(r['certified_lower_numerator'] for r in records)
    output={'format':VERSION,'scientific_scope':'Reference-only nominal GNM rank stability, NOT biological accuracy or NMR labels',
            'stage1_sha256':baseline['stage1_sha256'],'stage1_npz_sha256':npzsha,
            'base_certificate_sha256':{mode:hashlib.sha256(Path(file).read_bytes()).hexdigest() for mode,file in certfiles.items()},
            'baseline_case_report_sha256':hashlib.sha256(Path(case_report_path).read_bytes()).hexdigest(),
            'selected_contacts':[list(x) for x in SELECTED], 'n_cases':64,'n_pairs':2614,
            'n_pairwise_comparisons':TOTAL,'coupled_bound_type':'PSD_ORDER_CONTRACTION_EXACT_INTEGER',
            'universal_lower_numerator':minimum,'universal_lower_denominator':TOTAL,
            'cases':records}
    raw=(json.dumps(output,sort_keys=True,separators=(',',':'))+'\n').encode()
    Path(out).write_bytes(raw)
    result={'status':'PASS','n_cases':64,'new_witnesses':total_extra,'min_extra_case':min(len(r['extra_rank_pairs']) for r in records),
            'universal_lower_numerator':minimum,'baseline_lower_numerator':baseline['guaranteed_lower_numerator'],
            'case_certificate_bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'elapsed_seconds':round(time.perf_counter()-start,3)}
    print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for arg in ('stage1','certain','possible','nominal','oldcase','output'):parser.add_argument(arg)
    a=parser.parse_args()
    generate(a.stage1,a.certain,a.possible,a.nominal,a.oldcase,a.output)