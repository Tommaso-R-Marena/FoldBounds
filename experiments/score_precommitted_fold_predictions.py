"""Leakage-aware scoring interface for REAL fold-state population predictions.

Predictions CSV: id,alpha,beta,unfolded,model_name,model_git_sha,trained_until
Labels CSV: id,alpha,beta,unfolded. Independent experimental NMR estimates.
No data used for model training; this is an evaluation harness only.
The workflow must record a pre-prediction commit hash and train/test policy.
"""
import argparse,csv,json,math
from pathlib import Path

def load(path,need_extra=False):
    out={}
    with path.open(newline='') as f:
        rows=csv.DictReader(f)
        required={'id','alpha','beta','unfolded'}|({'model_name','model_git_sha','trained_until'} if need_extra else set())
        if not rows.fieldnames or not required.issubset(rows.fieldnames):raise ValueError('missing CSV headers')
        for row in rows:
            k=row['id'].strip()
            if not k or k in out:raise ValueError('blank/duplicate variant ID')
            vals=[float(row[x]) for x in ('alpha','beta','unfolded')]
            if any(not math.isfinite(v) or v<0 or v>1 for v in vals) or abs(sum(vals)-1)>1e-5:
                raise ValueError(f'invalid state distribution for {k}')
            if need_extra and any(not row[x].strip() for x in ('model_name','model_git_sha','trained_until')):
                raise ValueError('model provenance missing')
            out[k]=dict(zip(('alpha','beta','unfolded'),vals))
    return out

def score(pred,truth):
    if pred.keys()!=truth.keys():raise ValueError('predictions and labels refer to different variants')
    ids=sorted(pred)
    if not ids:raise ValueError('empty dataset')
    brier=0.0;ce=0.0;correct=0
    for k in ids:
        p=pred[k];t=truth[k]
        brier+=sum((p[c]-t[c])**2 for c in ('alpha','beta','unfolded'))
        ce-=sum(t[c]*math.log(max(p[c],1e-15)) for c in ('alpha','beta','unfolded'))
        correct+=max(p,key=p.get)==max(t,key=t.get)
    return {'variants':len(ids),'mean_multiclass_brier':brier/len(ids),
            'mean_cross_entropy_nats':ce/len(ids),'dominant_fold_accuracy':correct/len(ids),
            'explicit_warning':'Only valid as a held-out benchmark if predictions were frozen before labels used. Scoring alone does not prove that.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--predictions',type=Path,required=True);p.add_argument('--labels',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(score(load(a.predictions,need_extra=True),load(a.labels)),sort_keys=True,indent=2))
if __name__=='__main__':main()