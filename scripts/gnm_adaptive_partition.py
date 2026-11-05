#!/usr/bin/env python3
"""Fail-closed adaptive-contact GNM concordance certificates.

Novelty is the *application/workflow*, NOT Rayleigh monotonicity, case
splitting, or Sherman-Morrison identities.  For independent coordinate boxes,
exact integer geometry yields certain/optional/off edges; each contact mask
is contained in exactly one selected-edge case. Two GNM endpoint models
sandwich all remaining optional edges by the existing Rayleigh theorem.

A source-only, label-independent sensitivity heuristic chooses the split edges;
we do NOT use the NMR target to select these edges.  The independent replay
checks every endpoint integer witness, all score enclosures and target signs.

The sample variance across NMR conformer deposits is not molecular dynamics.
Full dataset-specific code is Python-checked; the Lean case-cover theorem
certifies the conditional mathematical composition, not byte decoding.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from hashlib import sha256
import io
import json
from pathlib import Path
from typing import Any
import zipfile
import numpy as np

import gnm_exact_certificates as g
from gnm_topology_boxes import classify_coord_box
from gnm_target_concordance import verified_zip, counts_fast, STAGE2_SHA

VERSION='FOLDBOUNDS_GNM_ADAPTIVE_PARTITION_V1'
EPSILON='0.05'
Q=1<<55
# The branch policy is deliberately chosen using reference data ONLY.  The
# selected indices are rechecked from the source-only algebraic sensitivity.
N_SPLITS=6

def source_only_edge_policy(coords, fixed, optional, root, count=N_SPLITS):
    n=len(coords)
    A, idx=g.grounded_laplacian(fixed,n,root)
    H=np.linalg.inv(np.asarray(A,float))  # untrusted proposal for branch ORDER only
    if not np.isfinite(H).all():raise ValueError('Invalid selection inverse')
    frozen_i=[];frozen_j=[]
    # Selection is pair-list-dependent, not target-label-dependent.
    # Caller provides exact stage1 eligible pairs below.
    return H,idx

def select_split_edges(coords, fixed, optional, root, pairs, count=N_SPLITS):
    H,idx=source_only_edge_policy(coords,fixed,optional,root,count)
    n=len(coords)
    ai=np.array([i for i,j in pairs],dtype=int);aj=np.array([j for i,j in pairs],dtype=int)
    choice=[]
    for k,(i,j) in enumerate(optional):
        u=np.zeros(n-1,float)
        if i!=root:u[idx[i]]+=1
        if j!=root:u[idx[j]]-=1
        v=H@u
        vi=np.zeros(n,float)
        for a in range(n):
            if a!=root:vi[a]=v[idx[a]]
        denom=float(1+u@v)
        if not denom>1:raise ValueError('Invalid grounded effective-resistance sensitivity')
        importance=float(np.sum((vi[ai]-vi[aj])**2/denom))
        if not np.isfinite(importance):raise ValueError('Nonfinite contact sensitivity')
        choice.append((importance,k))
    ranks=sorted(choice,key=lambda x:(-x[0],x[1]))[:count]
    return [int(k) for _,k in ranks], [float(v) for v,_ in ranks]

def frozen_target(stage1_arrays,stage2):
    files=verified_zip(stage2,STAGE2_SHA)
    with np.load(io.BytesIO(files['NMR_COORDINATES_AND_FROZEN_PREDICTIONS.npz']),allow_pickle=False) as y:
        if not np.array_equal(y['original_reference_coords'],stage1_arrays['reference_coords']):
            raise ValueError('ABSTAIN: stage1 coords changed')
        for first,second in [('pair_i','original_pair_i'),('pair_j','original_pair_j')]:
            if not np.array_equal(stage1_arrays[first],y[second]):raise ValueError('ABSTAIN: pair arrays changed')
        for method in ['GNM_resistance','ANM_weighted_inverse_square','ANM_ridge',
                       'packing_inverse','reference_distance','sequence_distance']:
            if not np.array_equal(stage1_arrays['prediction_'+method],y['frozen_prediction_'+method]):
                raise ValueError('ABSTAIN: frozen score mutated: '+method)
        if (y['target_coords'].shape!=(32,124,3) or
            not y['eligible_mask'].all() or
            not np.array_equal(y['evaluated_target_i'],stage1_arrays['pair_i']) or
            not np.array_equal(y['evaluated_target_j'],stage1_arrays['pair_j'])):
            raise ValueError('ABSTAIN: target mismatch')
        observed=np.var(np.linalg.norm(y['target_coords'][:,y['evaluated_target_i'],:]-
                                     y['target_coords'][:,y['evaluated_target_j'],:],axis=-1),
                        axis=0,ddof=1)
        if not np.isfinite(observed).all() or not np.array_equal(observed,y['observed_distance_variance']):
            raise ValueError('ABSTAIN: target variance changed')
    return observed

def make_graph_witness(edges,n,root,tree_B,pairs):
    A,idx=g.grounded_laplacian(edges,n,root)
    N=g.untrusted_inverse_proposal(A,Q)
    R=g.verify_integer_witness(A,N,Q)
    scores=[g.score_numerators(p,N,R,Q,tree_B,idx) for p in pairs]
    return N,scores

def make_certificate(stage1_zip,stage2_zip,output_path,epsilon=EPSILON,split_count=N_SPLITS):
    if epsilon != EPSILON:raise ValueError('ABSTAIN: this trial has exactly the frozen 0.05 A radius')
    if type(split_count) is not int or not 0<=split_count<=N_SPLITS:
        raise ValueError('ABSTAIN: invalid split count')
    arrays,sha_npz=g.open_sealed_stage1(stage1_zip)
    obs=frozen_target(arrays,stage2_zip)
    coords=arrays['reference_coords'];n=len(coords)
    xi,xj=arrays['pair_i'],arrays['pair_j'];x,pairs=g._validate_inputs(coords,xi,xj)
    if len(pairs)!=2614 or n!=124:raise ValueError('ABSTAIN: unexpected source dimensions')
    fixed,off,optional=classify_coord_box(coords,Fraction(epsilon))
    graph=[set() for _ in range(n)]
    for i,j in fixed: graph[i].add(j);graph[j].add(i)
    root,B,depth,parent=g.choose_tree(graph)
    selected,importance=select_split_edges(coords,fixed,optional,root,pairs,split_count)
    if len(set(selected))!=split_count:raise AssertionError('Bad split uniqueness')
    remaining=[e for k,e in enumerate(optional) if k not in selected]
    # 2^split_count disjoint cases, each with two extreme graphs.
    records=[];models=[]
    lower_global=Fraction(1,1);upper_global=Fraction(-1,1)
    for mask in range(1<<split_count):
        forced=[optional[selected[k]] for k in range(split_count) if (mask>>k)&1]
        # The endpoint with all remaining optional edges ON has the SMALLER response.
        graphs=[fixed+forced+remaining,fixed+forced]
        scores=[]
        graph_ids=[]
        for kind,edges in [('all_remaining_on',graphs[0]),('all_remaining_off',graphs[1])]:
            N,val=make_graph_witness(edges,n,root,B,pairs)
            model={'mask':mask,'kind':kind,'n_contacts':len(edges),
                   'integer_inverse_proposal':N}
            models.append(model)
            graph_ids.append(len(models)-1)
            scores.append(val)
        lower=[int(v[0]) for v in scores[0]]
        upper=[int(v[1]) for v in scores[1]]
        if any(a>b for a,b in zip(lower,upper)):
            raise AssertionError('Bad Rayleigh sandwich or numerical bound')
        counts=counts_fast(lower,upper,obs)
        this_lo=Fraction(counts['lower']);this_hi=Fraction(counts['upper'])
        lower_global=min(lower_global,this_lo)
        upper_global=max(upper_global,this_hi)
        records.append({'mask':mask,'all_remaining_on_model':graph_ids[0],
                        'all_remaining_off_model':graph_ids[1],
                        'target_concordance_counts':counts})
        print(f'case {mask+1}/{1<<split_count}: lower {float(this_lo):+.6f}',flush=True)
    ret={'format':VERSION,'stage1_sha256':g.STAGE1_SHA,'stage1_npz_sha256':sha_npz,
         'stage2_sha256':STAGE2_SHA,'epsilon_A':epsilon,
         'n_residues':n,'n_score_pairs':len(pairs),
         'always_on':len(fixed),'always_off':len(off),'optional_edges':len(optional),
         'root':root,'tree_parent':parent,'tree_B':B,
         'Q':Q,'selected_optional_indices':selected,'selected_edges':[list(optional[k]) for k in selected],
         'selection_importances':importance,
         'method':'all 2^k selected-contact cases, each with exact-integer endpoint certificates; every other optional contact is universally quantified via Rayleigh monotonicity',
         'lower_concordance':str(lower_global),'upper_concordance':str(upper_global),
         'cases':records,'models':models,
         'limitations':['Original blind 7RSA-to-2AAS primary outcome was FAIL; no rescue',
            'Target-aware concordance evaluated retrospectively; split policy source-only',
            'Coordinate boxes are hypothetical, not measurement-error distributions',
            'The target NMR conformers do not sample a physical equilibrium trajectory',
            'Lean proof of general inequalities is separate from Python-specific byte decoding']}
    path=Path(output_path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(ret,separators=(',',':'),sort_keys=True)+'\n')
    print('SOURCE-ONLY SPLIT EDGES',ret['selected_edges'])
    print('FINAL EXACT BOUND',ret['lower_concordance'],ret['upper_concordance'])
    return ret

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage1_zip');parser.add_argument('stage2_zip');parser.add_argument('output_json')
    parser.add_argument('--splits',type=int,default=N_SPLITS)
    args=parser.parse_args()
    make_certificate(args.stage1_zip,args.stage2_zip,args.output_json,split_count=args.splits)