#!/usr/bin/env python3
"""Independent all-integer replay of source-only *coupled* GNM order witnesses.

This program intentionally NEVER imports gnm_coupled_producer, numpy, sympy,
any NMR target data, or any floating-point calculation. It verifies the
exact graph+inverse identities at 128 endpoints, reproduces the PREVIOUS
independent interval counts, then checks every additional proposed pair by
squaring a rigorous PSD-contraction inequality over the inverse-order cone.
"""
from __future__ import annotations
import argparse,hashlib,json,math,time
from pathlib import Path
from fractions import Fraction
from gnm_exact_case_independent_replay import (
    loaded_source,independent_edges,read_checked_base,independent_update,
    validate_matrix,certified_count_offline,SELECTED,TOTAL
)

VERSION='FOLDBOUNDS_SOURCE_ONLY_COUPLED_RANK_V1'

def pair_terms(N,pair):
    a,b=pair
    def g(x):return None if x==45 else (x if x<45 else x-1)
    i,j=g(a),g(b)
    return ((N[i][i] if i is not None else 0)+
            (N[j][j] if j is not None else 0)-
            (2*N[i][j] if i is not None and j is not None else 0))


def signed_force(pair):
    a,b=pair
    out=[]
    for atom,sign in ((a,1),(b,-1)):
        if atom==45:continue
        out.append((atom if atom<45 else atom-1,sign))
    return out


def transfer(M,f,g):
    c=0
    for k,sgn in f:
        for j,sgnj in g:
            c+=sgn*sgnj*M[k][j]
    return c


def validate(stage1,lo_file,hi_file,nominal_file,previous_case_file,proposed_file):
    start=time.perf_counter()
    raw=Path(proposed_file).read_bytes()
    proposal=json.loads(raw)
    if proposal.get('format')!=VERSION:raise ValueError('ABSTAIN: wrong certificate format')
    base=json.loads(Path(previous_case_file).read_text())
    if base.get('format')!='FOLDBOUNDS_EXACT_GNM_SOURCE_ONLY_64_CASE_RANK_STABILITY_V1':
        raise ValueError('ABSTAIN: wrong predecessor certification')
    coords,pairs,_,npzsha=loaded_source(stage1)  # validates pinned original archive
    if (len(pairs)!=2614 or len(set(pairs))!=2614 or proposal.get('stage1_npz_sha256')!=npzsha or
        proposal.get('stage1_sha256')!=base.get('stage1_sha256') or
        proposal.get('baseline_case_report_sha256')!=hashlib.sha256(Path(previous_case_file).read_bytes()).hexdigest()):
        raise ValueError('ABSTAIN: changed original source or predecessor')
    files={'certain':lo_file,'possible':hi_file,'nominal':nominal_file}
    for mode,path in files.items():
        dig=hashlib.sha256(Path(path).read_bytes()).hexdigest()
        if proposal.get('base_certificate_sha256',{}).get(mode)!=dig or base.get(mode+'_certificate_sha256')!=dig:
            raise ValueError('ABSTAIN: changed original exact integer certificate')
    if (proposal.get('selected_contacts')!=[list(x) for x in SELECTED] or
        proposal.get('n_cases')!=64 or proposal.get('n_pairs')!=2614 or
        proposal.get('n_pairwise_comparisons')!=TOTAL or
        proposal.get('coupled_bound_type')!='PSD_ORDER_CONTRACTION_EXACT_INTEGER'):
        raise ValueError('ABSTAIN: case coverage metadata mismatch')
    sure=set(independent_edges(coords,'certain'));possibly=set(independent_edges(coords,'possible'))
    if len(sure)!=972 or len(possibly)!=1080 or len(possibly-sure)!=108:
        raise ValueError('ABSTAIN: wrong reference contact networks')
    initial={0:read_checked_base(stage1,lo_file,'certain')}
    final={0:read_checked_base(stage1,hi_file,'possible')}
    nomN,nomD=read_checked_base(stage1,nominal_file,'nominal')
    canonical_rank=sorted(range(2614),key=lambda k:pair_terms(nomN,pairs[k]))
    if len(set(pair_terms(nomN,p) for p in pairs))!=2614:
        raise ValueError('ABSTAIN: original exact scores have ties')
    sorted_pairs=[pairs[k] for k in canonical_rank]
    forces=[signed_force(p) for p in sorted_pairs]
    for mask in range(1,64):
        lowbit=mask&-mask;step=lowbit.bit_length()-1;prior=mask^lowbit
        initial[mask]=independent_update(*initial[prior],*SELECTED[step],1)
        final[mask]=independent_update(*final[prior],*SELECTED[step],-1)
    cases=proposal.get('cases')
    if not isinstance(cases,list) or len(cases)!=64 or [x.get('mask') for x in cases]!=list(range(64)):
        raise ValueError('ABSTAIN: case not exhaustive or duplicate')
    lowers=[];audited_extra=0
    for mask in range(64):
        loN,loD=initial[mask]
        complementary=(~mask)&63
        hiN,hiD=final[complementary]
        graph_lo=sure|{SELECTED[k] for k in range(6) if mask&(1<<k)}
        graph_hi=possibly-{SELECTED[k] for k in range(6) if complementary&(1<<k)}
        # Every integer entry of both 123 x 123 scaled inverse identities.
        validate_matrix(sorted(graph_lo),loN,loD,45)
        validate_matrix(sorted(graph_hi),hiN,hiD,45)
        # Uniform denominator D = loD*hiD. Naming here: lo[i] is the
        # minimal R from the dense graph, hi[i] maximal R from sparse.
        lower=[pair_terms(hiN,p)*loD for p in sorted_pairs]
        upper=[pair_terms(loN,p)*hiD for p in sorted_pairs]
        if any(l>u for l,u in zip(lower,upper)):
            raise ValueError('ABSTAIN: unphysical Rayleigh order')
        # Reconstruct old count using independent offline Fenwick ranking.
        yes,no=certified_count_offline(list(range(2614)),lower,upper)
        old=base['cases'][mask]
        if (old['mask']!=mask or old['definitely_correct']!=yes or
            old['definitely_wrong']!=no or old['ambiguous']!=TOTAL-yes-no):
            raise ValueError('ABSTAIN: old exact interval count failed independent replay')
        case=cases[mask]
        extras=case.get('extra_rank_pairs')
        if not isinstance(extras,list) or case.get('extra_count')!=len(extras):
            raise ValueError('ABSTAIN: wrong new witness length')
        # Compute exactly K=A_sparse^{-1}-A_dense^{-1}, with common D.
        # Subtraction is ordered in the graph Loewner cone K >= 0.
        K=[[loN[i][j]*hiD-hiN[i][j]*loD for j in range(123)] for i in range(123)]
        previous=(-1,-1)
        for indices in extras:
            if not (isinstance(indices,list) and len(indices)==2 and
                all(type(x) is int for x in indices)):
                raise ValueError('ABSTAIN: invalid rank pair representation')
            i,j=indices
            if not (0<=i<j<2614) or (i,j)<=previous:
                raise ValueError('ABSTAIN: noncanonical/duplicate rank witness')
            previous=(i,j)
            if lower[j]>upper[i]:
                raise ValueError('ABSTAIN: witness overlaps old certified comparison')
            diff=(upper[j]+lower[j])-(upper[i]+lower[i])
            if diff<=0:raise ValueError('ABSTAIN: no positive midpoint difference')
            di=upper[i]-lower[i];dj=upper[j]-lower[j]
            c=transfer(K,forces[i],forces[j])
            if di+dj+2*c<=0:
                raise ValueError('ABSTAIN: degenerate coupled direction')
            # Equivalent rigorous rational form, no square roots/floats.
            if diff*diff+4*c*c <= (di+dj)*(di+dj):
                raise ValueError('ABSTAIN: false coupled PSD-contraction witness')
        numerator=2*(yes+len(extras))-TOTAL
        if (case.get('baseline_correct_count')!=yes or
            case.get('certified_lower_numerator')!=numerator):
            raise ValueError('ABSTAIN: forged case conclusion')
        audited_extra+=len(extras);lowers.append(numerator)
        if mask%8==0:print('VERIFIED',mask,'new',len(extras),'lower numerator',numerator,'elapsed',round(time.perf_counter()-start,2),flush=True)
    if (proposal.get('universal_lower_numerator')!=min(lowers) or
        proposal.get('universal_lower_denominator')!=TOTAL):
        raise ValueError('ABSTAIN: forged universal guarantee')
    result={'status':'PASS','source_only':True,'NMR_target_accessed':False,
            'all_cases':64,'all_endpoint_inverse_checks':128,
            'witnessed_additional_pairs':audited_extra,
            'universal_signed_lower':str(Fraction(min(lowers),TOTAL)),
            'baseline_signed_lower':str(Fraction(base['guaranteed_lower_numerator'],TOTAL)),
            'witness_sha256':hashlib.sha256(raw).hexdigest(),
            'seconds':round(time.perf_counter()-start,3)}
    print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    a=argparse.ArgumentParser()
    for key in ('stage1','certain','possible','nominal','oldcase','newcase'):a.add_argument(key)
    x=a.parse_args();validate(x.stage1,x.certain,x.possible,x.nominal,x.oldcase,x.newcase)