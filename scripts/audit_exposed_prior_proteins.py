#!/usr/bin/env python3
"""Exploratory (NOT blind) ablation on previously revealed CspB/crambin labels.

Requires archive FoldBounds_Bulk_Rational_Audit_COMPLETE_20261104.zip. Must
not be used as independent validation or to tune any frozen stage 1.
"""
from pathlib import Path
from hashlib import sha256
import zipfile,io,json
import numpy as np
from hybrid_reference_model import rankdata
from hybrid_reference_model import frozen_methods

ARCHIVE_SHA='400dd9a242c1f34c71781f50079ad87f8a394eb03b3fc309aa71736aa2595923'

def evaluate(archive):
    raw=Path(archive).read_bytes()
    if sha256(raw).hexdigest()!=ARCHIVE_SHA:raise ValueError('Archive SHA mismatch')
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
        z=outer.read('FoldBounds_Bulk_Rational_Audit_20261104/FoldBounds_Prospective_SEND_BACK.zip')
    results={}
    with zipfile.ZipFile(io.BytesIO(z)) as fold:
        hashes=json.loads(fold.read('SHA256SUMS.json'))
        if set(hashes)!=set(fold.namelist())-{'SHA256SUMS.json'}:
            raise ValueError('Manifest missing or extra entries')
        for name,digest in hashes.items():
            if sha256(fold.read(name)).hexdigest()!=digest:
                raise ValueError('Integrity failure for '+name)
        for protein in ('CspB','Crambin'):
            with np.load(io.BytesIO(fold.read(protein+'_FROZEN_BEFORE_TARGET.npz')),allow_pickle=False) as f, \
                 np.load(io.BytesIO(fold.read(protein+'_EVALUATION_DATA.npz')),allow_pickle=False) as e:
                ii,jj=f['pair_i'],f['pair_j']
                xyz=f['reference_coords']; n=len(xyz)
                observed=e['observed_pair_variance']
                subset={(int(i),int(j)):k for k,(i,j) in enumerate(zip(ii,jj))}
                keep=[subset[(int(i),int(j))] for i,j in zip(e['eligible_pair_i'],e['eligible_pair_j'])]
                methods,_=frozen_methods(xyz,ii,jj)
                # The original stage1 classical controls must reproduce exactly.
                for m in ['GNM_resistance','ANM_weighted_inverse_square','ANM_ridge',
                          'packing_inverse','reference_distance','sequence_distance']:
                    x=methods[m]; target=f['prediction_'+m]
                    if np.max(abs(x-target))>1e-9:raise AssertionError('Old baseline changed: '+m)
                scores={m:float(np.corrcoef(rankdata(x[keep]),rankdata(observed))[0,1])
                        for m,x in methods.items()}
                best_old=max(scores[x] for x in ['GNM_resistance','ANM_weighted_inverse_square',
                          'ANM_ridge','packing_inverse','reference_distance','sequence_distance'])
                results[protein]={'n_reference_pairs':len(ii),'n_nmr_pairs':len(keep),
                     'scores':scores,'best_frozen_baseline':best_old,
                     'primary_margin':scores['GANM_half_isotropic']-best_old,
                     'frozen_gate_pass':scores['GANM_half_isotropic']>=.5 and scores['GANM_half_isotropic']-best_old>=.15}
    return results

if __name__=='__main__':
    import sys
    print(json.dumps(evaluate(sys.argv[1]),indent=2))