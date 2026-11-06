#!/usr/bin/env python3
"""Independent exact integer replay of 64-case nominal-rank stability proof.

No imports from the case producer, edge updater, numerical inverse proposal or
SymPy. Uses a SEPARATE integer Sherman-Morrison implementation and an offline
Fenwick sweep for strictly classified nominal-order comparisons.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from fractions import Fraction
from pathlib import Path
from gnm_exact_inverse_independent_replay import check as check_base,loaded_source,independent_edges,validate_matrix,strict_bigint

FORMAT='FOLDBOUNDS_EXACT_GNM_SOURCE_ONLY_64_CASE_RANK_STABILITY_V1'
SELECTED=((51,115),(1,9),(6,32),(43,121),(1,7),(52,76))
TOTAL=3415191


def read_checked_base(stage,path,mode):
    r=check_base(stage,path)
    if r['mode']!=mode:raise ValueError('ABSTAIN: wrong endpoint base graph')
    c=json.loads(Path(path).read_text())
    return [[strict_bigint(v) for v in row] for row in c['N']],strict_bigint(c['D'])


def independent_update(N,D,i,j,sgn):
    """Rebuild rank-one graph update using integer matrix *columns*.

    Output is independently rechecked by full graph A*N=D*I afterward.
    """
    if sgn not in (-1,1) or not 0<=i<j<124 or D<=0:
        raise ValueError('ABSTAIN: invalid update')
    index=lambda x: None if x==45 else (x if x<45 else x-1)
    a,b=index(i),index(j)
    x=[(N[k][a] if a is not None else 0)-(N[k][b] if b is not None else 0) for k in range(123)]
    kappa=(x[a] if a is not None else 0)-(x[b] if b is not None else 0)
    factor=D+sgn*kappa
    if kappa<=0 or factor<=0:
        raise ValueError('ABSTAIN: invalid resistance or disconnected removal')
    newD=D*factor
    newN=[]
    for row_index,row in enumerate(N):
        newN.append([factor*element-sgn*x[row_index]*x[col] for col,element in enumerate(row)])
    divisor=math.gcd(newD,*(v for row in newN for v in row))
    newN=[[v//divisor for v in row] for row in newN]
    return newN,newD//divisor


def val(N,D,u,v):
    i=None if u==45 else (u if u<45 else u-1)
    j=None if v==45 else (v if v<45 else v-1)
    n=(N[i][i] if i is not None else 0)+(N[j][j] if j is not None else 0)
    if i is not None and j is not None:n-=2*N[i][j]
    return Fraction(n,D)


class BIT:
    def __init__(self,n):self.data=[0]*(n+1)
    def add(self,i):
        i+=1
        while i<len(self.data):self.data[i]+=1;i+=i&-i
    def prefix(self,i):
        answer=0;i+=1
        while i>0:answer+=self.data[i];i-=i&-i
        return answer


def certified_count_offline(order,lower,upper):
    """Two offline Fenwick passes, independent of producer's bisect/insort."""
    rank={index:k for k,index in enumerate(order)}
    # Strict correct: prior upper_i < later lower_j.
    events=[(v,0,rank[i]) for i,v in enumerate(lower)] + [(v,1,rank[i]) for i,v in enumerate(upper)]
    events.sort(key=lambda e:(e[0],e[1])) # query BEFORE add at equality
    fw=BIT(len(order));correct=0
    for _,kind,idx in events:
        if kind==0:correct+=fw.prefix(idx-1)
        else:fw.add(idx)
    # Strict wrong: prior lower_i > later upper_j.
    # Descending by score; query BEFORE add on exact ties.
    events=[(v,0,rank[i]) for i,v in enumerate(upper)] + [(v,1,rank[i]) for i,v in enumerate(lower)]
    events.sort(key=lambda e:(-e[0],e[1]))
    fw=BIT(len(order));wrong=0
    for _,kind,idx in events:
        if kind==0:wrong+=fw.prefix(idx-1)
        else:fw.add(idx)
    if correct+wrong>TOTAL:raise ValueError('ABSTAIN: incompatible rank-order evidence')
    return correct,wrong


def replay(stage,certain_path,possible_path,nominal_path,reported):
    report=json.loads(Path(reported).read_text())
    if report.get('format')!=FORMAT:
        raise ValueError('ABSTAIN: wrong 64-case format')
    coords,pairs,_,npz=loaded_source(stage)
    sourcehash='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
    if report.get('stage1_sha256')!=sourcehash or report.get('stage1_npz_sha256')!=npz:
        raise ValueError('ABSTAIN: original reference metadata changed')
    f={'certain':certain_path,'possible':possible_path,'nominal':nominal_path}
    for mode,path in f.items():
        key=mode+'_certificate_sha256'
        if report.get(key)!=hashlib.sha256(Path(path).read_bytes()).hexdigest():
            raise ValueError('ABSTAIN: base exact-inverse certificate hash mismatch')
    sure=set(independent_edges(coords,'certain'))
    possible=set(independent_edges(coords,'possible'))
    nominal=set(independent_edges(coords,'nominal'))
    if (len(possible-sure)!=108 or len(nominal-sure)!=48 or
        not set(SELECTED).issubset(possible-sure)):
        raise ValueError('ABSTAIN: invalid reference-only optional edge classification')
    if (report.get('selected_contacts')!=[list(e) for e in SELECTED] or
        report.get('n_optional_contacts')!=108 or report.get('n_unselected_contacts')!=102 or
        report.get('n_reference_pairs')!=2614 or report.get('n_rank_comparisons')!=TOTAL or
        report.get('n_cases')!=64 or report.get('n_endpoint_graphs')!=128):
        raise ValueError('ABSTAIN: incomplete case coverage metadata')
    on={0:read_checked_base(stage,certain_path,'certain')}
    off={0:read_checked_base(stage,possible_path,'possible')}
    numN,numD=read_checked_base(stage,nominal_path,'nominal')
    nominal_values=[val(numN,numD,i,j) for i,j in pairs]
    if len(set(nominal_values))!=2614:raise ValueError('ABSTAIN: nominal rank tie')
    rank_order=sorted(range(2614),key=lambda k:nominal_values[k])
    for mask in range(1,64):
        one=(mask&(-mask))
        bit=one.bit_length()-1
        prev=mask^one
        on[mask]=independent_update(*on[prev],*SELECTED[bit],+1)
        off[mask]=independent_update(*off[prev],*SELECTED[bit],-1)
    by_mask=report.get('cases')
    if not isinstance(by_mask,list) or len(by_mask)!=64 or [c.get('mask') for c in by_mask]!=list(range(64)):
        raise ValueError('ABSTAIN: missing, duplicate or unsorted case')
    low_bounds=[];upper_bounds=[]
    for mask in range(64):
        lowerN,lowerD=on[mask]
        subset=sure|{SELECTED[b] for b in range(6) if mask&(1<<b)}
        validate_matrix(sorted(subset),lowerN,lowerD,45)
        complement=(~mask)&63
        upperN,upperD=off[complement]
        superset=possible-{SELECTED[b] for b in range(6) if complement&(1<<b)}
        validate_matrix(sorted(superset),upperN,upperD,45)
        lower=[val(upperN,upperD,i,j) for i,j in pairs]
        upper=[val(lowerN,lowerD,i,j) for i,j in pairs]
        if any(a>b for a,b in zip(lower,upper)):
            raise ValueError('ABSTAIN: impossible monotone endpoints')
        correct,wrong=certified_count_offline(rank_order,lower,upper)
        c=by_mask[mask]
        if (c.get('definitely_correct')!=correct or c.get('definitely_wrong')!=wrong or
            c.get('ambiguous')!=TOTAL-correct-wrong or
            c.get('guaranteed_signed_lower_numerator')!=2*correct-TOTAL or
            c.get('guaranteed_signed_upper_numerator')!=TOTAL-2*wrong):
            raise ValueError('ABSTAIN: claimed case rank certificate incorrect')
        low_bounds.append(2*correct-TOTAL);upper_bounds.append(TOTAL-2*wrong)
    low=min(low_bounds);hi=max(upper_bounds)
    if (report.get('guaranteed_lower_numerator')!=low or
        report.get('guaranteed_upper_numerator')!=hi or
        report.get('lower_signed_concordance')!=str(Fraction(low,TOTAL)) or
        report.get('upper_signed_concordance')!=str(Fraction(hi,TOTAL))):
        raise ValueError('ABSTAIN: incorrect universal exact rank bound')
    return {'status':'PASS','n_endpoint_scaled_inverses_checked':128,
            'all_64_cases_checked':True,'all_contact_worlds_covered':True,
            'target_labels_accessed':False,
            'universal_lower':str(Fraction(low,TOTAL)),
            'universal_upper':str(Fraction(hi,TOTAL)),
            'case_certificate_sha256':hashlib.sha256(Path(reported).read_bytes()).hexdigest()}

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name in ('stage1','certain','possible','nominal','case_report'):ap.add_argument(name)
    args=ap.parse_args();print(json.dumps(replay(args.stage1,args.certain,args.possible,args.nominal,args.case_report),indent=2))