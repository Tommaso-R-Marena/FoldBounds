#!/usr/bin/env python3
"""Exact GNM structural-uncertainty concordance bound against revealed 2AAS.

Input BLOBs are pinned by SHA256; no target optimization. Given certified GNM
score envelopes (integer numerators / common Q^2 denominator), and the
binary64 target-ensemble distance variance, compute a signed rank agreement
interval for ALL allowable reference-coordinate contact networks.

All score-comparison certification is exact integer arithmetic. Sorting the
already revealed target's fixed floating-point variances is an exact ordering
of the archived IEEE 754 values, but those samples are NOT equilibrium dynamics.
No independence assumption on residue-pair comparisons.

This returns Kendall-style average directional agreement over non-tied target
pairs, NOT conventional tau-b, a statistical confidence interval, or biological
out-of-sample proof. It may be conservative when interval signs correlate.
"""
from __future__ import annotations
from bisect import bisect_left, bisect_right, insort
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import argparse
import io
import json
import zipfile

import numpy as np

STAGE1_SHA='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
STAGE2_SHA='4f1ea752ada45453a44ea9128b6e7a43f4070e6824062b796a27144c340dfc6c'


def verified_zip(path, expected):
    data=Path(path).read_bytes()
    if sha256(data).hexdigest()!=expected:
        raise ValueError('ABSTAIN: original archive SHA256 mismatch')
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        if z.testzip() is not None:raise ValueError('ABSTAIN: ZIP CRC failure')
        manifest=json.loads(z.read('SHA256SUMS.json'))
        if set(manifest)!=set(z.namelist())-{'SHA256SUMS.json'}:
            raise ValueError('ABSTAIN: incomplete archive manifest')
        for k,v in manifest.items():
            if sha256(z.read(k)).hexdigest()!=v:
                raise ValueError('ABSTAIN: internal archive digest mismatch')
        return {k:z.read(k) for k in z.namelist()}


def counts_fast(lo,hi,observed):
    """Correct comparison count with O(n^2) total ordered insertion worst-case.

    Sorting target groups before update handles tied target values without
    arbitrarily declaring model concordance. Interval endpoints are integers.
    """
    n=len(lo)
    if len(hi)!=n or len(observed)!=n or n<2:
        raise ValueError('ABSTAIN: input dimension invalid')
    if any(type(x) is not int for x in lo+hi):
        raise ValueError('ABSTAIN: integer interval endpoints required')
    if any(a>b for a,b in zip(lo,hi)):
        raise ValueError('ABSTAIN: inverted interval')
    obs=np.asarray(observed,dtype=np.float64)
    if not np.isfinite(obs).all():raise ValueError('ABSTAIN: invalid target')
    order=sorted(range(n),key=lambda i:float(obs[i]))
    high=[];low=[]
    C=D=U=tied=0
    ix=0
    while ix<n:
        end=ix+1
        while end<n and obs[order[end]]==obs[order[ix]]:end+=1
        # Prior groups have strictly lower target observations.
        for k in range(ix,end):
            i=order[k]
            c=bisect_left(high,lo[i])
            d=len(low)-bisect_right(low,hi[i])
            if c+d>len(high):raise AssertionError('inconsistent interval comparison')
            C+=c;D+=d;U+=len(high)-c-d
        tied+=(end-ix)*(end-ix-1)//2
        for k in range(ix,end):
            i=order[k]
            insort(low,lo[i]);insort(high,hi[i])
        ix=end
    m=C+D+U
    if m<1 or m+tied!=n*(n-1)//2:
        raise ValueError('ABSTAIN: no target-distinguishable pairs')
    return {'n_pairs':n,'total_comparisons':n*(n-1)//2,
            'target_ties':tied,'distinct_target_comparisons':m,
            'guaranteed_concordant':C,'guaranteed_discordant':D,
            'uncertain':U,'lower':str(Fraction(C-D-U,m)),
            'upper':str(Fraction(C-D+U,m)),
            'certain_fraction':str(Fraction(C+D,m))}


def counts_slow(lo,hi,observed):
    n=len(lo);C=D=U=t=0
    for i in range(n):
        for j in range(i+1,n):
            if observed[i]==observed[j]:t+=1;continue
            a,b=(i,j) if observed[i]>observed[j] else (j,i)
            if lo[a]>hi[b]:C+=1
            elif hi[a]<lo[b]:D+=1
            else:U+=1
    return C,D,U,t


def evaluate(stage1,stage2,certificate):
    a=verified_zip(stage1,STAGE1_SHA)
    b=verified_zip(stage2,STAGE2_SHA)
    with np.load(io.BytesIO(a['FROZEN_REFERENCE_ONLY.npz']),allow_pickle=False) as x, \
         np.load(io.BytesIO(b['NMR_COORDINATES_AND_FROZEN_PREDICTIONS.npz']),allow_pickle=False) as y:
        source_i,source_j=x['pair_i'],x['pair_j']
        if len(source_i)!=2614 or not (
            np.array_equal(source_i,y['original_pair_i']) and
            np.array_equal(source_j,y['original_pair_j']) and
            np.array_equal(x['reference_coords'],y['original_reference_coords']) and
            np.array_equal(x['prediction_GNM_resistance'],y['frozen_prediction_GNM_resistance'])):
            raise ValueError('ABSTAIN: frozen predictions/coordinates changed')
        if not y['eligible_mask'].all() or y['target_coords'].shape!=(32,124,3):
            raise ValueError('ABSTAIN: target eligibility or conformers invalid')
        ii,jj=y['evaluated_target_i'],y['evaluated_target_j']
        if not np.array_equal(ii,source_i) or not np.array_equal(jj,source_j):
            raise ValueError('ABSTAIN: target/ref alignment changed')
        obs=np.var(np.linalg.norm(y['target_coords'][:,ii,:]-
                                     y['target_coords'][:,jj,:],axis=-1),axis=0,ddof=1)
        if not np.array_equal(obs,y['observed_distance_variance']):
            raise ValueError('ABSTAIN: observed distances do not reproduce')
        scores=x['prediction_GNM_resistance'].copy()
    # Before examining the target-ranked signs, independently recheck EVERY
    # coordinate-box contact status, endpoint witness, residual, and rank count.
    # The verifier has a separate coordinate decoder and direct edge-energy
    # implementation, avoiding self-certification by the envelope producer.
    from gnm_independent_box_replay import replay as independently_verify_boxes
    independently_verify_boxes(stage1,certificate)
    record=json.loads(Path(certificate).read_bytes())
    if record['format']!='GNM_TOPOLOGY_BOX_HULL_V1' or record['source_archive_sha256']!=STAGE1_SHA:
        raise ValueError('ABSTAIN: unrecognized source of score envelopes')
    if len(record['scores'])<1:raise ValueError('ABSTAIN: no coordinate uncertainty envelopes')
    output={'format':'RNASE_REVEALED_GNM_AGREEMENT_V1',
            'source_stage1_sha256':STAGE1_SHA,'source_stage2_sha256':STAGE2_SHA,
            'structural_envelope_sha256':sha256(Path(certificate).read_bytes()).hexdigest(),
            'description':'signed direction agreement on all target-distinct pair comparisons; exact interval endpoint comparisons but only conditional on externally sound structural envelopes',
            'coordinate_intervals':[]}
    # Raw intervals reconstructed from exact integer numerators; soundness
    # of the source certificate must be validated by separate replay script.
    for item in record['scores']:
        lo=item['lower_score_numerators'];hi=item['upper_score_numerators']
        if not(len(lo)==len(hi)==len(scores)):
            raise ValueError('ABSTAIN: score vector lengths incompatible')
        counts=counts_fast(lo,hi,obs)
        counts['epsilon_A']=item['epsilon'];counts['source_endpoint_denominator']=item['denominator']
        output['coordinate_intervals'].append(counts)
    # Nominal GNM score ordering against exposed target (exploratory).
    nominal=counts_fast([int.from_bytes(np.asarray(float(v),dtype='<f8').tobytes(),'little') for v in scores],
                        [int.from_bytes(np.asarray(float(v),dtype='<f8').tobytes(),'little') for v in scores],obs)
    # Using IEEE754 bit patterns as monotonically ordered ints is only valid for
    # *nonnegative* floats, which effective resistances must be. This diagnostic
    # is not used to establish the certified envelope bounds.
    if not np.isfinite(scores).all() or (scores<=0).any():
        raise ValueError('ABSTAIN: invalid nominal resistance scores')
    output['nominal_concordance']=nominal
    return output


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('stage1');ap.add_argument('stage2');ap.add_argument('structural_json');ap.add_argument('output_json')
    args=ap.parse_args()
    out=evaluate(args.stage1,args.stage2,args.structural_json)
    Path(args.output_json).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    for r in out['coordinate_intervals']:
        print('epsilon',r['epsilon_A'],'certain',r['certain_fraction'],
              'concordance lower',r['lower'],'upper',r['upper'])
    print('nominal',out['nominal_concordance']['lower'])


if __name__=='__main__':main()