#!/usr/bin/env python3
"""Sound heterogeneous-depth adaptive contact-case GNM witness producer.

All split contacts and order are frozen from reference 7RSA alone BEFORE 2AAS
is opened.  NMR labels are used only to stop traversing retrospective branches,
not to propose spring strengths, contacts, or a new protein motion model.

Each terminal leaf represents ALL completions of its unselected contact bits.
Graph endpoint witnesses are reused between parent and child branches, and
all eventual numerical claims are rechecked by a separate verifier.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from pathlib import Path
import json, hashlib
import numpy as np

import gnm_exact_certificates as g
from gnm_topology_boxes import classify_coord_box
from gnm_target_concordance import counts_fast, STAGE2_SHA
from gnm_adaptive_partition import select_split_edges, make_graph_witness, frozen_target

VERSION='FOLDBOUNDS_GNM_PREFIX_TREE_V1'
EPSILON='0.05'
Q=1<<55
MAX_SPLITS=6


def make_tree(stage1, stage2, output, threshold=Fraction(1,10000), max_splits=MAX_SPLITS):
    if not isinstance(threshold,Fraction) or not (0<threshold<1):
        raise ValueError('ABSTAIN: rational positive threshold required')
    if type(max_splits) is not int or not 0<=max_splits<=6:
        raise ValueError('ABSTAIN: max splits out of bound')
    arrays,npz_sha=g.open_sealed_stage1(stage1)
    coords=arrays['reference_coords'];n=len(coords)
    xyz,pairs=g._validate_inputs(coords,arrays['pair_i'],arrays['pair_j'])
    if n!=124 or len(pairs)!=2614:
        raise ValueError('ABSTAIN: changed reference dimensions')
    fixed,off,optional=classify_coord_box(coords,Fraction(EPSILON))
    if (len(fixed),len(off),len(optional))!=(972,6546,108):
        raise ValueError('ABSTAIN: coordinate cutoff classification changed')
    neighbors=[set() for _ in range(n)]
    for i,j in fixed:
        neighbors[i].add(j);neighbors[j].add(i)
    root,B,depth,parent=g.choose_tree(neighbors)
    selected,importance=select_split_edges(coords,fixed,optional,root,pairs,MAX_SPLITS)
    # A sealed, source-only selection reproduced from the existing research PR.
    if selected != [65,3,12,55,2,66]:
        raise ValueError('ABSTAIN: original source-only edge policy changed')
    print('FROZEN SOURCE-ONLY EDGE ORDER', [optional[i] for i in selected],flush=True)
    print('FROZEN STOP THRESHOLD',str(threshold), 'MAX SPLITS',max_splits,flush=True)
    observed=frozen_target(arrays,stage2)

    models=[]
    cached={}
    visited=0
    leaves=[]

    def endpoint(edges):
        key=tuple(sorted(edges))
        if len(set(key))!=len(key):
            raise ValueError('ABSTAIN: duplicate edge in graph')
        if key in cached:return cached[key]
        N,scores=make_graph_witness(edges,n,root,B,pairs)
        mid=len(models)
        models.append({'n_contacts':len(edges),'integer_inverse_proposal':N})
        cached[key]=(mid,scores)
        return cached[key]

    def descend(d,mask):
        nonlocal visited
        visited+=1
        forced_on={optional[selected[k]] for k in range(d) if (mask>>k)&1}
        forced_off={optional[selected[k]] for k in range(d) if not (mask>>k)&1}
        low_graph=fixed+[e for e in optional if e not in forced_off]
        high_graph=fixed+[e for e in optional if e in forced_on]
        mid_low, scores_lo = endpoint(low_graph)
        mid_high, scores_hi = endpoint(high_graph)
        lows=[v[0] for v in scores_lo]
        highs=[v[1] for v in scores_hi]
        if any(lo>hi for lo,hi in zip(lows,highs)):
            raise ValueError('ABSTAIN: Rayleigh endpoint ordering failed')
        counts=counts_fast(lows,highs,observed)
        score_lo=Fraction(counts['lower'])
        if score_lo>=threshold or d==max_splits:
            leaves.append({'depth':d,'prefix_mask':mask,
                'all_unfixed_on_model':mid_low,
                'all_unfixed_off_model':mid_high,
                'target_concordance_counts':counts})
            return
        descend(d+1,mask)
        descend(d+1,mask|(1<<d))

    descend(0,0)
    lows=[Fraction(x['target_concordance_counts']['lower']) for x in leaves]
    highs=[Fraction(x['target_concordance_counts']['upper']) for x in leaves]
    certified=min(lows)
    if not (certified>=threshold and certified>0):
        raise ValueError('ABSTAIN: requested threshold not certified')
    if len(models)!=2*len(leaves):
        raise AssertionError('Cache/tree graph accounting invalid')
    all_mask_leaves=[]
    for full in range(1<<MAX_SPLITS):
        matching=[x for x in leaves if full&((1<<x['depth'])-1)==x['prefix_mask']]
        if len(matching)!=1:raise AssertionError('Incomplete or duplicate tree coverage')
        all_mask_leaves.append((matching[0]['depth'],matching[0]['prefix_mask']))
    report={
      'format':VERSION,'stage1_sha256':g.STAGE1_SHA,'stage1_npz_sha256':npz_sha,
      'stage2_sha256':STAGE2_SHA,'epsilon_A':EPSILON,'n_residues':n,'n_score_pairs':len(pairs),
      'always_on':len(fixed),'always_off':len(off),'optional_edges':len(optional),
      'root':root,'tree_parent':parent,'tree_B':B,'Q':Q,
      'source_only_selected_indices':selected,
      'source_only_selected_edges':[list(optional[k]) for k in selected],
      'threshold':str(threshold),'max_splits':max_splits,
      'n_nodes_visited':visited,'n_terminal_leaves':len(leaves),
      'lower_concordance':str(certified),'upper_concordance':str(max(highs)),
      'leaves':leaves,'models':models,
      'guarantee':'Every binary optional contact assignment extends exactly one prefix leaf; each leaf bounds all other contacts using the two extreme Laplacians. Monotone Rayleigh theorem gives all 2614 pair score enclosures and signed lower concordance.',
      'limitations':['Original blind RNase A trial failed; this is an exploratory retrospective certification',
        'Leaf stopping uses already revealed NMR target labels',
        'The selected contact edge order is frozen from the 7RSA reference only',
        'Uncalibrated independent coordinate box and deposited NMR ensemble',
        'Large protein integer replay is Python checked; universal tree theorem separately Lean checked']}
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,sort_keys=True,separators=(',',':'))+'\n',encoding='utf8')
    digest=hashlib.sha256(output.read_bytes()).hexdigest()
    print('PREFIX-TREE CERTIFICATE',json.dumps({
       'sha256':digest,'threshold':str(threshold),'lower':str(certified),
       'upper':str(max(highs)),'leaf_cases':len(leaves),'graph_witnesses':len(models),
       'tree_nodes_visited':visited,'full_mask_count':1<<MAX_SPLITS,
       'file':str(output)},sort_keys=True))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('stage1');p.add_argument('stage2');p.add_argument('output')
    p.add_argument('--threshold',default='1/10000')
    p.add_argument('--max-splits',type=int,default=6)
    args=p.parse_args()
    make_tree(args.stage1,args.stage2,args.output,Fraction(args.threshold),args.max_splits)