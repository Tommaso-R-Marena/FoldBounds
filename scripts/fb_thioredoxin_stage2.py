#!/usr/bin/env python3
"""Strict Stage2 runner for pre-committed 2TRX -> 1XOA blind benchmark.

DO NOT RUN until Stage 1 archive has been returned and SHA256 has been
independently confirmed. Then call `fb_thio_evaluate <SHA256>` in PyMOL.
A SHA256 string is required, NEVER inferred from the untrusted archive.
No target access occurs until the external SHA, manifest, protocol and
score arrays have been validated. No model coefficients are recalculated.
"""
from pathlib import Path
import hashlib
import io
import json
import math
import re
import traceback
import zipfile
import numpy as np

STAGE2_VERSION = 'FOLDBOUNDS_THIO_2TRX_1XOA_STAGE2_20261105_V1'
STAGE1_PREFIX = 'FoldBounds_THIO_STAGE1_20261105'
EXPECTED_STAGE1_VERSION = 'FOLDBOUNDS_THIO_2TRX_STAGE1_20261105_V1'
EXPECTED_TARGET = '1XOA'
EXPECTED_REFERENCE = '2TRX crystal chain A state 1'
BASELINES = ('GNM_resistance','ANM_weighted_inverse_square','ANM_ridge',
             'packing_inverse','reference_distance','sequence_distance')
PRIMARY = 'GANM_half_isotropic'
SECONDARY = 'GNM_ANM_rank_fusion_75_25'
MODELS = BASELINES + (PRIMARY,SECONDARY)


def sha(b):return hashlib.sha256(b).hexdigest()


def _parse_stage1(expected_sha, path):
    if not isinstance(expected_sha,str) or not re.fullmatch(r'[0-9a-f]{64}',expected_sha):
        raise ValueError('ABSTAIN: pass the exact 64-character PRETARGET Stage1 SHA256')
    blob=Path(path).read_bytes()
    if sha(blob)!=expected_sha:raise ValueError('ABSTAIN: Stage1 SHA256 does not match the externally pinned value')
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        if z.testzip():raise ValueError('ABSTAIN: input ZIP corrupt')
        m=json.loads(z.read('SHA256SUMS.json'))
        if set(m)!=set(z.namelist())-{'SHA256SUMS.json'}:
            raise ValueError('ABSTAIN: bad Stage1 ZIP manifest')
        for filename,digest in m.items():
            if sha(z.read(filename))!=digest:raise ValueError('ABSTAIN: Stage1 manifest mismatch '+filename)
        p=json.loads(z.read('FROZEN_PROTOCOL.json'))
        if (p.get('version')!=EXPECTED_STAGE1_VERSION or
            p.get('reference')!=EXPECTED_REFERENCE or
            p.get('sealed_target')!='1XOA NMR all conformers' or
            p.get('primary')!=PRIMARY or
            p.get('secondary')!=SECONDARY or
            p.get('alpha_isotropic')!=0.5 or
            p.get('rank_fusion_GNM_weight')!=0.75 or
            p.get('baselines')!=list(BASELINES)):
            raise ValueError('ABSTAIN: Stage1 parameters differ from frozen protocol')
        source=z.read('fb_thioredoxin_stage1.py')
        if sha(source)!=p.get('code_sha256'):
            raise ValueError('ABSTAIN: Stage1 script provenance mismatch')
        frozen=z.read('FROZEN_THIO_REFERENCE.npz')
        with np.load(io.BytesIO(frozen),allow_pickle=False) as obj:
            arrays={key:obj[key].copy() for key in obj.files}
    if not set('prediction_'+x for x in MODELS).issubset(arrays):
        raise ValueError('ABSTAIN: missing a frozen predictor')
    res=arrays['residues'];names=arrays['names'];xyz=arrays['reference_coords']
    i,j=arrays['pair_i'],arrays['pair_j']
    if (res.ndim!=1 or len(res)<45 or len(set(map(int,res)))!=len(res) or
        names.shape!=res.shape or xyz.shape!=(len(res),3) or
        len(i)!=len(j) or len(i)<100 or not np.isfinite(xyz).all() or
        np.any(i<0) or np.any(j>=len(res)) or np.any(i>=j)):
        raise ValueError('ABSTAIN: invalid Stage1 reference arrays')
    for name in MODELS:
        x=arrays['prediction_'+name]
        if x.shape != i.shape or not np.isfinite(x).all():
            raise ValueError('ABSTAIN: invalid predictor '+name)
    return arrays,p,sha(frozen)


def _read_target(cmd,obj):
    states=int(cmd.count_states(obj))
    if states<10:raise ValueError('ABSTAIN: fewer than 10 target NMR conformers')
    frames=[];types=[]
    for state in range(1,states+1):
        model=cmd.get_model('(%s) and chain A and name CA'%obj,state=state)
        current={};names={}
        for atom in model.atom:
            label=str(atom.resi).strip()
            if not re.fullmatch(r'[0-9]+',label):continue
            r=int(label)
            if r in current:raise ValueError('ABSTAIN: duplicate NMR C-alpha residue')
            current[r]=tuple(map(float,atom.coord))
            names[r]=str(atom.resn).upper()
        if len(current)<45:raise ValueError('ABSTAIN: too few target C-alphas')
        frames.append(current);types.append(names)
    common=sorted(set.intersection(*(set(frame) for frame in frames)))
    if len(common)<45:raise ValueError('ABSTAIN: NMR models have few common residues')
    for r in common:
        if len({names[r] for names in types})!=1:
            raise ValueError('ABSTAIN: target sequence changes between NMR models')
    xyz=np.asarray([[frame[k] for k in common] for frame in frames],float)
    if not np.isfinite(xyz).all():raise ValueError('ABSTAIN: nonfinite NMR coordinates')
    return np.asarray(common,np.int64),np.asarray([types[0][r] for r in common],dtype='U4'),xyz


def rankdata(v):
    x=np.asarray(v,float)
    order=np.argsort(x,kind='stable');result=np.empty(len(x),float)
    i=0
    while i<len(x):
        j=i+1
        while j<len(x) and x[order[j]]==x[order[i]]:j+=1
        result[order[i:j]]=.5*(i+j-1)+1
        i=j
    return result


def spearman(v,observed):
    x=rankdata(v);y=rankdata(observed)
    x-=x.mean();y-=y.mean()
    d=np.linalg.norm(x)*np.linalg.norm(y)
    return float((x@y)/d) if d else float('nan')


def auroc(v,observed):
    v=np.asarray(v,float); observed=np.asarray(observed,float)
    k=max(1,math.ceil(.1*len(v)))
    hot=np.zeros(len(v),bool)
    hot[np.argsort(observed,kind='stable')[-k:]]=True
    r=rankdata(v);p=int(hot.sum());n=len(v)-p
    return float((r[hot].sum()-p*(p+1)/2)/(p*n))


def _evaluate(arr,target_res,target_names,target_coord):
    ref=arr['residues'];rn=arr['names'];ii=arr['pair_i'];jj=arr['pair_j']
    tmap={int(r):k for k,r in enumerate(target_res)}
    match=np.asarray([int(r) in tmap and str(rn[k])==str(target_names[tmap[int(r)]])
                      for k,r in enumerate(ref)],bool)
    allowed=match[ii]&match[jj]
    if allowed.sum()<100 or allowed.mean()<.85:
        raise ValueError('ABSTAIN: fewer than 85% reference pairs retain residue identity')
    ti=np.asarray([tmap[int(ref[k])] for k in ii[allowed]],dtype=int)
    tj=np.asarray([tmap[int(ref[k])] for k in jj[allowed]],dtype=int)
    diff=target_coord[:,ti,:]-target_coord[:,tj,:]
    obs=np.linalg.norm(diff,axis=2).var(axis=0,ddof=1)
    if not np.isfinite(obs).all() or not np.isfinite(target_coord).all():
        raise ValueError('ABSTAIN: invalid observed NMR variances')
    metrics={}
    for name in MODELS:
        pr=arr['prediction_'+name][allowed]
        metrics[name]={'spearman':spearman(pr,obs),'auroc_top10':auroc(pr,obs)}
        if not all(np.isfinite(list(metrics[name].values()))):
            raise ValueError('ABSTAIN: undefined correlation or AUROC')
    best=max(metrics[m]['spearman'] for m in BASELINES)
    rho=metrics[PRIMARY]['spearman']
    result=dict(stage2_version=STAGE2_VERSION,reference='2TRX',target='1XOA',
        conformers=len(target_coord),matched_residues=int(match.sum()),
        retained_pairs=int(allowed.sum()),reference_pairs=len(ii),
        retained_fraction=float(allowed.mean()),method_metrics=metrics,
        primary=PRIMARY,primary_spearman=rho,best_baseline_spearman=best,
        primary_margin=rho-best,
        strict_primary_pass=bool(rho>=.50 and rho-best>=.15),
        status='PASS' if rho>=.50 and rho-best>=.15 else 'FAIL',
        interpretation='Association with deposited NMR distance-variance, not true equilibrium or kinetics')
    arrays=dict(target_residues=target_res,target_names=target_names,target_coords=target_coord,
        retained_pair_mask=allowed,reference_residues=ref,reference_names=rn,
        reference_coords=arr['reference_coords'],pair_i=ii,pair_j=jj,
        target_i=ti,target_j=tj,observed_distance_variance=obs,
        **{'frozen_prediction_'+name:arr['prediction_'+name] for name in MODELS})
    return result,arrays


def fb_thio_evaluate(expected_sha='',stage1_zip=''):
    from pymol import cmd
    print('\n=== FOLDBOUNDS 2TRX->1XOA: PRECOMMITTED STAGE2 ===')
    folder=Path.home()/'Downloads'
    if not folder.is_dir():folder=Path.home()
    if stage1_zip:
        src=Path(stage1_zip)
    else:
        candidates=sorted(folder.glob('FoldBounds_THIO_STAGE1_20261105*SEND_BACK.zip'))
        if len(candidates)!=1:raise ValueError('ABSTAIN: Stage1 ZIP missing/ambiguous; provide exact stage1_zip path')
        src=candidates[0]
    frozen,protocol,npz_sha=_parse_stage1(expected_sha,str(src))
    print('BEFORE target fetch, Stage1 SHA verified:',expected_sha)
    print('Frozen primary G-ANM alpha 0.5, all six baselines, secondary rank fusion')
    obj='fb_thio_1xoa_target'
    if obj in cmd.get_names('objects'):cmd.delete(obj)
    cmd.fetch(EXPECTED_TARGET.lower(),name=obj,type='pdb',async_=0)
    if obj not in cmd.get_names('objects') or cmd.count_atoms(obj)<100:
        raise RuntimeError('ABSTAIN: could not load 1XOA')
    target_res,target_names,target_coord=_read_target(cmd,obj)
    result,arrays=_evaluate(frozen,target_res,target_names,target_coord)
    result['stage1_archive_sha256']=expected_sha
    result['stage1_npz_sha256']=npz_sha
    out=folder/('FoldBounds_THIO_STAGE2_20261105_'+expected_sha[:12])
    if out.exists():
        for k in range(2,100):
            proposal=Path(str(out)+'_%02d'%k)
            if not proposal.exists():out=proposal;break
    out.mkdir(parents=True)
    np.savez_compressed(out/'TARGET_NMR_EVALUATION.npz',**arrays)
    (out/'FROZEN_EVALUATION.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    (out/'READ_ME_FIRST.txt').write_text('Return whole ZIP for independent audit. Do not tune after target.\n')
    (out/'SHA256SUMS.json').write_text(json.dumps({p.name:sha(p.read_bytes()) for p in sorted(out.iterdir()) if p.is_file()},indent=2)+'\n')
    archive=folder/('FoldBounds_THIO_STAGE2_20261105_'+expected_sha[:12]+'_SEND_BACK.zip')
    if archive.exists():raise FileExistsError('ABSTAIN: output archive would overwrite an earlier result')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:
        for p in sorted(out.iterdir()):z.write(p,p.name)
    print('\n=== PRIMARY RESULTS (prespecified) ===')
    for key,m in result['method_metrics'].items():print('  %-32s rho %+0.6f AUROC %0.6f'%(key,m['spearman'],m['auroc_top10']))
    print('PRIMARY:',result['status'],'margin:',round(result['primary_margin'],6))
    print('Results SHA256:',sha(archive.read_bytes()))
    print('RETURN THIS ZIP:',archive)
    return str(archive)


try:
    from pymol import cmd
    cmd.extend('fb_thio_evaluate',fb_thio_evaluate)
    print('FoldBounds Stage2 loaded. DO NOT RUN until Stage1 SHA externally pinned.')
except ImportError:
    pass

if __name__=='__main__':print('Import inside PyMOL only; Stage1 must be independently verified first.')