#!/usr/bin/env python3
"""Exact model-self-consistency across all contact graphs permitted by a box.

No NMR target or ML labels. All resistance interval endpoints are exact
fractions from A*N=D*I witnesses, instead of approximate residual enclosures.
Case partition is defined by six reference-only optional contacts; all other
102 contacts may independently be on or off. Bounds cover a superset of
physically realizable coordinate boxes, so guarantees are conservative.
"""
from __future__ import annotations
import argparse,bisect,hashlib,json,time
from pathlib import Path
from fractions import Fraction
from gnm_exact_inverse_producer import archive_source,exact_edges,grounded,check_scaled_inverse,pair_num
from gnm_exact_inverse_independent_replay import check as verify_base
from gnm_exact_edge_updates import update_contact_signed

SELECTED=((51,115),(1,9),(6,32),(43,121),(1,7),(52,76))
TOTAL=2614*2613//2
VERSION='FOLDBOUNDS_EXACT_GNM_SOURCE_ONLY_64_CASE_RANK_STABILITY_V1'

def load_matrix(path):
    c=json.loads(Path(path).read_text())
    return [[int(v) for v in row] for row in c['N']],int(c['D'])

def score_vectors(N,D,pairs,idx):
    return [Fraction(pair_num(i,j,N,idx),D) for i,j in pairs]

def count_nominal_pairwise_agreement(order,lowers,uppers):
    """Every world in this case satisfies lo[i]<=true[i]<=hi[i].

    Count strict proven same-order and strict proven reversed-order pairs
    relative to a nominal reference ranking, using O(P log P) bisection.
    """
    previous_uppers=[];previous_lowers=[]
    correct=reverse=0
    for i in order:
        correct+=bisect.bisect_left(previous_uppers,lowers[i])
        reverse+=len(previous_lowers)-bisect.bisect_right(previous_lowers,uppers[i])
        bisect.insort(previous_uppers,uppers[i]);bisect.insort(previous_lowers,lowers[i])
    if not 0<=correct+reverse<=TOTAL:
        raise ValueError('ABSTAIN: impossible rank concordance certificate')
    return correct,reverse

def run(stage1,certain_path,possible_path,nominal_path,out):
    t=time.perf_counter()
    coords,pairs,_,npzsha=archive_source(stage1)
    for file,mode in ((certain_path,'certain'),(possible_path,'possible'),(nominal_path,'nominal')):
        if verify_base(stage1,file)['mode']!=mode:raise ValueError('ABSTAIN: source-only witness mode mismatch')
    sure=set(exact_edges(coords,'certain'))
    possible=set(exact_edges(coords,'possible'))
    nominal=set(exact_edges(coords,'nominal'))
    if any(e not in possible-sure for e in SELECTED):
        raise ValueError('ABSTAIN: selected contact is not reference-optional')
    if len(possible-sure)!=108 or len(nominal-sure)!=48:
        raise ValueError('ABSTAIN: frozen source changed')
    low_states={0:load_matrix(certain_path)}
    upper_states={0:load_matrix(possible_path)}
    for mask in range(1,64):
        lsb=mask&-mask
        k=lsb.bit_length()-1
        prev=mask^lsb
        low_states[mask]=update_contact_signed(*low_states[prev],*SELECTED[k],45,+1)
        upper_states[mask]=update_contact_signed(*upper_states[prev],*SELECTED[k],45,-1)
    _,idx=grounded(sorted(sure))
    nominalN,nominalD=load_matrix(nominal_path)
    nominal_scores=score_vectors(nominalN,nominalD,pairs,idx)
    if len(set(nominal_scores))!=2614:raise ValueError('ABSTAIN: nominal ranks include exact ties')
    nominal_order=sorted(range(2614),key=lambda k:nominal_scores[k])
    result=[]
    for mask in range(64):
        # Case includes selected-on bits; the complement is forced-off.
        # Sparse lower endpoint has certain-on + chosen bits only.
        ln,ld=low_states[mask]
        active=sure | {SELECTED[i] for i in range(6) if mask>>i&1}
        A,_=grounded(sorted(active))
        check_scaled_inverse(A,ln,ld)
        # Complete upper endpoint has all possible edges except forced-offs.
        off=(~mask)&63
        un,ud=upper_states[off]
        active_upper=possible-{SELECTED[i] for i in range(6) if off>>i&1}
        Au,_=grounded(sorted(active_upper))
        check_scaled_inverse(Au,un,ud)
        lo=score_vectors(un,ud,pairs,idx)
        hi=score_vectors(ln,ld,pairs,idx)
        if any(a>b for a,b in zip(lo,hi)):
            raise ValueError('ABSTAIN: monotone resistance bounds reversed')
        yes,no=count_nominal_pairwise_agreement(nominal_order,lo,hi)
        result.append({'mask':mask,'definitely_correct':yes,'definitely_wrong':no,
                       'ambiguous':TOTAL-yes-no,
                       'guaranteed_signed_lower_numerator':2*yes-TOTAL,
                       'guaranteed_signed_upper_numerator':TOTAL-2*no})
    lower=min(r['guaranteed_signed_lower_numerator'] for r in result)
    upper=max(r['guaranteed_signed_upper_numerator'] for r in result)
    output={
        'format':VERSION,'stage1_sha256':'ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff',
        'stage1_npz_sha256':npzsha,
        'certain_certificate_sha256':hashlib.sha256(Path(certain_path).read_bytes()).hexdigest(),
        'possible_certificate_sha256':hashlib.sha256(Path(possible_path).read_bytes()).hexdigest(),
        'nominal_certificate_sha256':hashlib.sha256(Path(nominal_path).read_bytes()).hexdigest(),
        'selected_contacts':[list(e) for e in SELECTED],
        'n_optional_contacts':108,'n_unselected_contacts':102,
        'n_reference_pairs':2614,'n_rank_comparisons':TOTAL,
        'n_cases':64,'n_endpoint_graphs':128,
        'lower_signed_concordance':str(Fraction(lower,TOTAL)),
        'upper_signed_concordance':str(Fraction(upper,TOTAL)),
        'guaranteed_lower_numerator':lower,'guaranteed_upper_numerator':upper,
        'cases':result,
        'scientific_scope':'Source-only exact-model rank stability; NOT NMR accuracy, dynamics, or a blind prediction.'
    }
    raw=(json.dumps(output,sort_keys=True,separators=(',',':'))+'\n').encode()
    Path(out).write_bytes(raw)
    return {'status':'PASS','elapsed_seconds':round(time.perf_counter()-t,3),
            'n_cases':64,'n_endpoint_graphs_checked':128,
            'lower':output['lower_signed_concordance'],
            'upper':output['upper_signed_concordance'],
            'case_certificate_bytes':len(raw),
            'case_certificate_sha256':hashlib.sha256(raw).hexdigest(),
            'first_case':result[0],'last_case':result[-1]}

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('stage1');ap.add_argument('certain');ap.add_argument('possible');ap.add_argument('nominal');ap.add_argument('out')
    args=ap.parse_args();print(json.dumps(run(args.stage1,args.certain,args.possible,args.nominal,args.out),indent=2))