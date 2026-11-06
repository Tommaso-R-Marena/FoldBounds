#!/usr/bin/env python3
"""Independently replay geometry-strengthened coupled GNM rank bounds.

The verifier does not import either NEW producer or the exploratory analysis.
It uses pre-existing exact-integer grounded-inverse verification components
from the already audited PR72/PR78. All observations are from the Stage1
source structure; no Stage2/NMR/thioredoxin target is opened.

On the full run, it validates old coupled witnesses first, then rederives all
64 contact-case upper/lower response envelopes and the proposed conclusion.
For development/performance only, --skip-existing-audit may be used after a
separate successful existing-certificate audit, but is NOT an end-to-end check.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from fractions import Fraction
from pathlib import Path
from gnm_geometry_independent_replay import replay as verify_geometry
from gnm_coupled_independent_replay import validate as verify_existing, pair_terms
from gnm_exact_case_independent_replay import (
    loaded_source,read_checked_base,independent_update,certified_count_offline,
    validate_matrix,SELECTED,TOTAL,val)

FORMAT='FOLDBOUNDS_GNM_GEOMETRY_AWARE_COUPLED_64CASE_V1'

def run(stage1,base_dir,geometry_file,report_file,lo=0,hi=64,skip_existing_audit=False):
    t=time.time();root=Path(base_dir)
    certain=root/'gnm_certain.json';possible=root/'gnm_possible.json';nominal=root/'gnm_nominal.json'
    baseline=root/'gnm_source_exact_64cases.json';coupled=root/'gnm_coupled_witnesses.json'
    references=[certain,possible,nominal,baseline,coupled]
    proposed=json.loads(Path(report_file).read_text())
    if proposed.get('format')!=FORMAT:raise ValueError('wrong signed certificate format')
    for path in references+[Path(geometry_file)]:
        if proposed['source_artifact_sha256'].get(path.name)!=hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError('changed source witness: '+path.name)
    verified_geometry=verify_geometry(stage1,geometry_file)
    if verified_geometry['at_least_excluded_fraction']!='32010617/33554432':raise ValueError('wrong geometric exclusion result')
    if not skip_existing_audit:
        result=verify_existing(stage1,certain,possible,nominal,baseline,coupled)
        if result['witnessed_additional_pairs']!=849112:raise ValueError('changed previous coupled guarantee')
    coords,pairs,_,_=loaded_source(stage1)
    original=json.loads(baseline.read_text());oldcoupled=json.loads(coupled.read_text())
    canonical=read_checked_base(stage1,nominal,'nominal')
    order=sorted(range(len(pairs)),key=lambda k:pair_terms(canonical[0],pairs[k]))
    if len(order)!=2614:raise ValueError('not all pair forces')
    forces=[pairs[k] for k in order]
    Lcache={0:read_checked_base(stage1,certain,'certain')}
    Ucache={0:read_checked_base(stage1,possible,'possible')}
    for mask in range(1,64):
        bit=(mask&-mask).bit_length()-1;prev=mask^(1<<bit)
        Lcache[mask]=independent_update(*Lcache[prev],*SELECTED[bit],1)
        Ucache[mask]=independent_update(*Ucache[prev],*SELECTED[bit],-1)
    geo=json.loads(Path(geometry_file).read_text())
    optional=[tuple(e) for e in geo['optional_edges']]
    # Only monotone SAT clauses are used here; all 20 are independently
    # checked and the full 95.4% excluded-fraction result is retained.
    monotone=[]
    for cl in geo['forbidden_patterns']:
        state=cl['forbidden_on_bits']
        if all(v==state[0] for v in state):
            monotone.append((bool(state[0]),[optional[i] for i in cl['edge_ids']]))
    if len(monotone)!=7:raise ValueError('missing geometry clauses')
    out=[]
    case_rows=proposed.get('cases')
    if not isinstance(case_rows,list) or len(case_rows)!=64:
        raise ValueError('case coverage must be exactly 64')
    if not 0<=lo<hi<=64:raise ValueError('invalid mask window')
    for mask in range(lo,hi):
        sparse=Lcache[mask];dense=Ucache[(~mask)&63]
        # Explicit rational response score interval R_dense<=R_world<=R_sparse.
        lower=[val(dense[0],dense[1],*pair) for pair in forces]
        upper=[val(sparse[0],sparse[1],*pair) for pair in forces]
        selected_on={SELECTED[k] for k in range(6) if mask&(1<<k)}
        selected_off={SELECTED[k] for k in range(6) if not mask&(1<<k)}
        use_upper=use_lower=0
        for forbidden_on,group in monotone:
            if not forbidden_on:
                # Not all group members can be off. Every geometrically
                # realizable world must have at least one group member on.
                if any(e in selected_on for e in group):continue
                choices=[e for e in group if e not in selected_off]
                if not choices:raise ValueError('geometrically impossible mask, not implemented')
                candidates=[]
                for e in choices:
                    N,D=independent_update(*sparse,*e,1)
                    candidates.append([val(N,D,*p) for p in forces])
                # For each choice e: G_world contains G_sparse+e, so
                # R_world<=R_sparse+e; OR of choices => max over them.
                upper=[min(u,max(scores[k] for scores in candidates)) for k,u in enumerate(upper)]
                use_upper+=1
            else:
                # Not all group members can be on. At least one missing edge
                # establishes a lower bound from the corresponding dense-e.
                if any(e in selected_off for e in group):continue
                choices=[e for e in group if e not in selected_on]
                if not choices:raise ValueError('geometrically impossible mask, not implemented')
                candidates=[]
                for e in choices:
                    N,D=independent_update(*dense,*e,-1)
                    candidates.append([val(N,D,*p) for p in forces])
                lower=[max(l,min(scores[k] for scores in candidates)) for k,l in enumerate(lower)]
                use_lower+=1
        if any(l>u for l,u in zip(lower,upper)):
            raise ValueError('inconsistent tightened resistance interval')
        definite,reverse=certified_count_offline(list(range(2614)),lower,upper)
        extras=oldcoupled['cases'][mask]['extra_rank_pairs']
        # Union: recompute whether every previous coupled witness has now
        # become a simple interval witness; count overlapping claims once.
        not_yet_interval=sum(upper[i]>=lower[j] for i,j in extras)
        olddef=original['cases'][mask]['definitely_correct']+len(extras)
        newdef=definite+not_yet_interval
        if newdef<olddef or newdef>TOTAL:raise ValueError('tightened proof weakened or exceeds all pairs')
        row={'mask':mask,'old_certified_orders':olddef,'new_certified_orders':newdef,
             'new_interval_orders':definite,'remaining_coupled_only_orders':not_yet_interval,
             'net_new_orders':newdef-olddef,'upper_nogoods_used':use_upper,
             'lower_nogoods_used':use_lower}
        if case_rows[mask]!=row:raise ValueError('forged case result mask '+str(mask))
        out.append(row)
    if lo==0 and hi==64:
        minnumer=min(2*c['new_certified_orders']-TOTAL for c in out)
        previous=min(2*c['old_certified_orders']-TOTAL for c in out)
        if (minnumer!=proposed.get('universal_lower_numerator') or
            previous!=proposed.get('previous_lower_numerator') or
            proposed.get('universal_lower_denominator')!=TOTAL or
            proposed.get('total_additional_casewise_orders')!=sum(c['net_new_orders'] for c in out)):
            raise ValueError('forged global guarantee')
    return {'status':'PASS','all_case_rows_in_window':hi-lo,'window':[lo,hi],
            'no_NMR_labels':True,'verified_exact_geometry':True,
            'new_case_order_witnesses':sum(c['net_new_orders'] for c in out),
            'elapsed_seconds':round(time.time()-t,3)}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage1');parser.add_argument('base_dir');parser.add_argument('geometry');parser.add_argument('report')
    parser.add_argument('--start',type=int,default=0);parser.add_argument('--end',type=int,default=64)
    parser.add_argument('--skip-existing-audit',action='store_true')
    a=parser.parse_args()
    print(json.dumps(run(a.stage1,a.base_dir,a.geometry,a.report,a.start,a.end,a.skip_existing_audit),indent=2))