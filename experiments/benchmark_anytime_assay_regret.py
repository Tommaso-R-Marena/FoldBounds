"""Reproduce joint, anytime-valid S669 assay-budget regret results.

No raw third-party data is included in this PR; --download fetches the exact
published GeoStab Git commit, authenticated by the original Git blob hash.
All evidence/rank arithmetic is exact Fraction and published Decimal.
"""
from __future__ import annotations

import argparse
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from urllib.request import urlopen

from foldbounds.assay_regret_certificates import (
    SOURCE_REPO, SOURCE_COMMIT, SOURCE_PATH, SOURCE_GIT_BLOB,
    git_blob_sha1, parse_s669, _shuffle,
)
from foldbounds.anytime_assay_regret import (
    DEFAULT_BUDGETS, assay_group_evidence, propose_certificate,
    signed_conformal_cutoff, retrospectively_audit,
)
from foldbounds.anytime_assay_regret_verifier import verify_anytime_certificate


def _rounded(x: Fraction | None) -> float | None:
    return None if x is None else round(float(x),9)


def study(rows, *, seed=20261105, motif_length=None,
          lam=Fraction(1), epsilon=Fraction(0),
          check_independent=True, alpha=Fraction(1,10)):
    group_data=assay_group_evidence(rows,DEFAULT_BUDGETS,motif_length=motif_length)
    keys=tuple(sorted(group_data))
    if len(keys)<5:
        raise ValueError('insufficient eligible complete protein groups')
    assignment={g:i%5 for i,g in enumerate(_shuffle(keys,seed))}
    budgets=DEFAULT_BUDGETS
    bins={str(b):{'mean_joint':Fraction(0),'mean_bonf':Fraction(0),
                  'mean_single':Fraction(0),'joint_groups_covered':0,
                  'bonf_groups_covered':0,'single_groups_covered':0}
          for b in budgets}
    covered_joint=0
    stop_stats={t:{'counts':{str(b):0 for b in budgets},'met_target':0,
                   'groups_covered':0,'total_assays':0,'total_candidates':0}
                for t in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(5,2))}
    fold_reports=[]
    for f in range(5):
        cal=tuple(sorted(g for g in keys if assignment[g]!=f))
        test=tuple(sorted(g for g in keys if assignment[g]==f))
        cert=propose_certificate(rows,cal,budgets=budgets,
                                 lambda_margin=lam,alpha=alpha,motif_length=motif_length)
        if check_independent:
            verify_anytime_certificate(cert,rows)
        n=len(cal)
        rank_b={}
        for b in budgets:
            scores=[max(x.actual_regret-lam*x.observed_margin
                        for x in group_data[g] if x.budget==b)
                    for g in cal]
            rank_b[b]={
                'single': signed_conformal_cutoff(scores,alpha),
                'bonf': signed_conformal_cutoff(scores,alpha/Fraction(len(budgets))),
            }
        group_hits=0
        for g in test:
            records=group_data[g]
            budget_bounds={}
            any_failed=False
            for b in budgets:
                items=tuple(x for x in records if x.budget==b)
                if len(items)==0:
                    raise AssertionError('missing calibration or held-out horizon')
                if cert.cutoff is None:
                    raise ValueError('unexpected abstention: specify more calibration groups')
                eta=max((Fraction(2)+lam)*epsilon/w for w in cert.weights)
                weight=cert.weights[budgets.index(b)]
                b_joint=[];b_solo=[];b_bonf=[]
                for item in items:
                    if item.actual_regret is None:
                        raise AssertionError('held-out truth required for retrospective assessment')
                    q1=rank_b[b]['single']
                    qb=rank_b[b]['bonf']
                    # For benchmark epsilon=0, these are exact. With epsilon>0
                    # the additive extension is a hypothetical worst-case
                    # sensitivity, not an estimate of laboratory noise.
                    b_joint.append(max(Fraction(0),lam*(item.observed_margin+epsilon)
                                       +weight*(cert.cutoff+eta)))
                    b_solo.append(None if q1 is None else max(
                        Fraction(0),lam*(item.observed_margin+epsilon)
                        +q1+(Fraction(2)+lam)*epsilon))
                    b_bonf.append(None if qb is None else max(
                        Fraction(0),lam*(item.observed_margin+epsilon)
                        +qb+(Fraction(2)+lam)*epsilon))
                t_truth=max(x.actual_regret for x in items)
                good_joint=all(x.actual_regret<=val for x,val in zip(items,b_joint))
                good_single=all(val is None or x.actual_regret<=val for x,val in zip(items,b_solo))
                good_bonf=all(val is None or x.actual_regret<=val for x,val in zip(items,b_bonf))
                any_failed |= not good_joint
                key=str(b)
                bins[key]['mean_joint']+=sum(b_joint,Fraction(0))/len(b_joint)
                bins[key]['mean_bonf']=(None if any(v is None for v in b_bonf) or bins[key]['mean_bonf'] is None
                                         else bins[key]['mean_bonf']+sum(b_bonf,Fraction(0))/len(b_bonf))
                bins[key]['mean_single']=(None if any(v is None for v in b_solo) or bins[key]['mean_single'] is None
                                           else bins[key]['mean_single']+sum(b_solo,Fraction(0))/len(b_solo))
                bins[key]['joint_groups_covered']+=int(good_joint)
                bins[key]['bonf_groups_covered']+=int(good_bonf)
                bins[key]['single_groups_covered']+=int(good_single)
                budget_bounds[b]=max(b_joint)
            covered_joint+=int(not any_failed)
            group_hits+=int(not any_failed)
            for target,record in stop_stats.items():
                meeting=next((b for b in budgets if budget_bounds[b]<=target),None)
                stop=meeting if meeting is not None else budgets[-1]
                record['met_target']+=int(meeting is not None)
                record['counts'][str(stop)]+=1
                stop_items=[x for x in records if x.budget==stop]
                record['groups_covered']+=int(all(x.actual_regret<=budget_bounds[stop]
                                                 for x in stop_items))
                record['total_assays']+=sum(x.assayed for x in stop_items)
                record['total_candidates']+=sum(x.total for x in stop_items)
        fold_reports.append({'fold':f,'calibration':len(cal),'test':len(test),
                             'q_joint':_rounded(cert.cutoff),'joint_groups_covered':group_hits})
    for b in budgets:
        for field in ('mean_joint','mean_bonf','mean_single'):
            bins[str(b)][field]=_rounded(bins[str(b)][field]/len(keys)) if bins[str(b)][field] is not None else None
    return {'seed':seed,'group_unit':'motif-12' if motif_length==12 else 'exact-wildtype',
            'protein_or_motif_groups':len(keys),'nominal_simultaneous_alpha':str(alpha),
            'joint_groups_covered':covered_joint,
            'joint_groups_coverage':covered_joint/len(keys),
            'budget_results':bins,
            'adaptive_stops':{str(t):v for t,v in stop_stats.items()},
            'folds':fold_reports}


def main():
    cli=argparse.ArgumentParser()
    src=cli.add_mutually_exclusive_group(required=True)
    src.add_argument('--csv',type=Path)
    src.add_argument('--download',action='store_true')
    cli.add_argument('--repeat',type=int,default=60)
    cli.add_argument('--seed',type=int,default=20261105)
    cli.add_argument('--epsilon',default='0',help='hypothetical maximum absolute assay error, kcal/mol; Fraction')
    cli.add_argument('--no-independent-audit',action='store_true')
    cli.add_argument('--output',type=Path)
    opt=cli.parse_args()
    if not 1<=opt.repeat<=60:
        cli.error('--repeat must be in [1,60]')
    epsilon=Fraction(opt.epsilon)
    if epsilon<0:
        cli.error('--epsilon must be nonnegative')
    if opt.download:
        uri=f'https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{SOURCE_PATH}'
        with urlopen(uri,timeout=60) as response:
            data=response.read(2_000_000)
    else:
        data=opt.csv.read_bytes()
    if git_blob_sha1(data)!=SOURCE_GIT_BLOB:
        raise ValueError('unrecognized S669 source; exact pinned Git blob required')
    records=parse_s669(data)
    report={
        'source':{'repository':SOURCE_REPO,'commit':SOURCE_COMMIT,
                  'path':SOURCE_PATH,'git_blob_sha1':SOURCE_GIT_BLOB},
        'n_mutations':len(records),
        'assumptions':['complete protein/motif group exchangeability',
                        'published 3D predictions fixed before evaluation',
                        'matched reported pH and temperature within panels',
                        'calibration truth never used to choose candidates',
                        'epsilon is hypothetical deterministic maximum error, not fitted from S669'],
        'budgets':[str(x) for x in DEFAULT_BUDGETS],
        'lambda_penalty':'1',
        'uniform_horizon_weights':True,
        'epsilon_kcal_mol':str(epsilon),
        'results':[],
    }
    for seed in range(opt.seed,opt.seed+opt.repeat):
        record={'seed':seed}
        for mode,mlen in (('protein',None),('motif12',12)):
            record[mode]=study(records,seed=seed,motif_length=mlen,
                     epsilon=epsilon,check_independent=not opt.no_independent_audit)
        report['results'].append(record)
    output=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if opt.output:
        opt.output.write_text(output,encoding='utf-8')
    print(output)


if __name__=='__main__':
    main()
