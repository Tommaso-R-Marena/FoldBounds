#!/usr/bin/env python3
"""Exploratory dual-model consensus on previously revealed protein ensembles.

No new scientific validation or selective-inference significance claim.
Tests if ANM and GNM agree on the sign of each pairwise ranking of mobility
scores, and whether that mutually agreed sign tracks the target NMR variance.
"""
import io,zipfile,json,hashlib
from pathlib import Path
import numpy as np

STAGE1_SHA='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
STAGE2_SHA='4f1ea752ada45453a44ea9128b6e7a43f4070e6824062b796a27144c340dfc6c'
PRIOR_SHA='400dd9a242c1f34c71781f50079ad87f8a394eb03b3fc309aa71736aa2595923'


def pair_sign(a,b):
    return (a>b).astype(np.int8)-(a<b).astype(np.int8)


def concordance(gnm,anm,obs):
    n=len(obs);ii,jj=np.triu_indices(n,1)
    g=pair_sign(gnm[ii],gnm[jj]);a=pair_sign(anm[ii],anm[jj]);t=pair_sign(obs[ii],obs[jj])
    mask=t!=0
    g=g[mask];a=a[mask];t=t[mask]
    cons=(g!=0)&(a!=0)&(g==a)
    disagree=(g!=0)&(a!=0)&(g!=a)
    out={'n_predictions':n,'n_target_distinct_pairs':len(t),'consensus_coverage':float(cons.mean()),'model_disagreement':float(disagree.mean()),
         'gnm_signed_concordance':float(np.mean(g*t)), 'anm_signed_concordance':float(np.mean(a*t)),
         'consensus_signed_concordance':float(np.mean(g[cons]*t[cons])) if cons.any() else None,
         'consensus_accuracy':float(np.mean(g[cons]==t[cons])) if cons.any() else None,
         'gnm_nominal_accuracy':float(np.mean(g==t)), 'anm_nominal_accuracy':float(np.mean(a==t)),
         'n_consensus_pairs':int(cons.sum())}
    return out


def run(stage1,stage2,archive):
    if hashlib.sha256(Path(stage1).read_bytes()).hexdigest()!=STAGE1_SHA:raise ValueError('Stage1 SHA')
    if hashlib.sha256(Path(stage2).read_bytes()).hexdigest()!=STAGE2_SHA:raise ValueError('Stage2 SHA')
    if hashlib.sha256(Path(archive).read_bytes()).hexdigest()!=PRIOR_SHA:raise ValueError('Prior SHA')
    with zipfile.ZipFile(stage1) as a,zipfile.ZipFile(stage2) as b:
        with np.load(io.BytesIO(a.read('FROZEN_REFERENCE_ONLY.npz')),allow_pickle=False) as x, \
             np.load(io.BytesIO(b.read('NMR_COORDINATES_AND_FROZEN_PREDICTIONS.npz')),allow_pickle=False) as y:
            out={'RNase A':concordance(x['prediction_GNM_resistance'],x['prediction_ANM_weighted_inverse_square'],y['observed_distance_variance'])}
    with zipfile.ZipFile(archive) as a:
        inner=a.read('FoldBounds_Bulk_Rational_Audit_20261104/FoldBounds_Prospective_SEND_BACK.zip')
        with zipfile.ZipFile(io.BytesIO(inner)) as z:
            for p in ['CspB','Crambin']:
                with np.load(io.BytesIO(z.read(p+'_FROZEN_BEFORE_TARGET.npz')),allow_pickle=False) as x, \
                     np.load(io.BytesIO(z.read(p+'_EVALUATION_DATA.npz')),allow_pickle=False) as y:
                    ix={(int(i),int(j)):k for k,(i,j) in enumerate(zip(x['pair_i'],x['pair_j']))}
                    keep=[ix[(int(i),int(j))] for i,j in zip(y['eligible_pair_i'],y['eligible_pair_j'])]
                    out[p]=concordance(x['prediction_GNM_resistance'][keep],x['prediction_ANM_weighted_inverse_square'][keep],y['observed_pair_variance'])
    return out

if __name__=='__main__':
    import sys
    print(json.dumps(run(*sys.argv[1:]),indent=2))