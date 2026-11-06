#!/usr/bin/env python3
"""Independent exact-integer replay of variable-depth GNM prefix-tree bounds.

Does NOT import the certificate producer, the floating inverse routines, or
source-only selection algorithm. It rechecks the sealed source ZIPs, exact
binary64 interval geometry, exhaustive prefix-code coverage, every reused
endpoint graph/witness, and all 3.4 million signed target comparisons per leaf.
The generic graph-response and prefix-cover theorems are separately in Lean.
"""
from __future__ import annotations
import argparse
import io,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction
import numpy as np
import gnm_adaptive_independent_replay as check

VERSION='FOLDBOUNDS_GNM_PREFIX_TREE_V1'
SELECTED=[65,3,12,55,2,66]


def verify_prefix_cover(prefixes, max_depth):
    """Cover every selected-bit assignment exactly once, with no overlap.

    This is a byte-level executable check of the finite prefix-code syntax.
    The general tree covering theorem lives separately in Lean.
    """
    if type(max_depth)!=int or not 0<=max_depth<=16:
        raise ValueError('ABSTAIN: invalid maximum mask depth')
    seen=set()
    for d,m in prefixes:
        if type(d)!=int or type(m)!=int or not 0<=d<=max_depth or not 0<=m<1<<d:
            raise ValueError('ABSTAIN: invalid prefix length or bits')
        if (d,m) in seen:raise ValueError('ABSTAIN: duplicate prefix')
        seen.add((d,m))
    owners=[]
    for world in range(1<<max_depth):
        matches=[i for i,(d,m) in enumerate(prefixes)
                 if world&((1<<d)-1)==m]
        if len(matches)!=1:
            raise ValueError('ABSTAIN: missing or overlapping Boolean contact worlds')
        owners.append(matches[0])
    if set(owners)!=set(range(len(prefixes))):
        raise ValueError('ABSTAIN: unused prefix')
    return owners


def verify(stage1, stage2, certificate):
    s1=check.sealed_zip(stage1,check.EXPECTED_STAGE1)
    s2=check.sealed_zip(stage2,check.EXPECTED_STAGE2)
    with np.load(io.BytesIO(s1['FROZEN_REFERENCE_ONLY.npz']),allow_pickle=False) as f,\
         np.load(io.BytesIO(s2['NMR_COORDINATES_AND_FROZEN_PREDICTIONS.npz']),allow_pickle=False) as t:
        coords=f['reference_coords'].copy();n=len(coords)
        pairs=list(zip(map(int,f['pair_i']),map(int,f['pair_j'])))
        if n!=124 or len(pairs)!=2614 or len(set(pairs))!=len(pairs):
            raise ValueError('ABSTAIN: changed Stage1 dimensions or repeated pair')
        if (not np.array_equal(coords,t['original_reference_coords']) or
            not np.array_equal(f['pair_i'],t['evaluated_target_i']) or
            not np.array_equal(f['pair_j'],t['evaluated_target_j']) or
            not t['eligible_mask'].all() or
            not t['matched_residues'].all() or
            t['target_coords'].shape!=(32,n,3) or
            not np.array_equal(f['residues'],t['target_residues'])):
            raise ValueError('ABSTAIN: target provenance or structure dimensions changed')
        if not all(np.array_equal(f['prediction_'+name],t['frozen_prediction_'+name])
           for name in ['GNM_resistance','ANM_weighted_inverse_square','ANM_ridge',
                        'packing_inverse','reference_distance','sequence_distance']):
            raise ValueError('ABSTAIN: original frozen model predictions altered')
        obs=np.var(np.linalg.norm(t['target_coords'][:,t['evaluated_target_i'],:]-
                         t['target_coords'][:,t['evaluated_target_j'],:],axis=-1),axis=0,ddof=1)
        if not np.array_equal(obs,t['observed_distance_variance']):
            raise ValueError('ABSTAIN: modified target pair distances/variance')
    raw=Path(certificate).read_bytes()
    cert=json.loads(raw)
    if (cert.get('format')!=VERSION or
        cert.get('stage1_sha256')!=check.EXPECTED_STAGE1 or
        cert.get('stage1_npz_sha256')!=sha256(s1['FROZEN_REFERENCE_ONLY.npz']).hexdigest() or
        cert.get('stage2_sha256')!=check.EXPECTED_STAGE2 or
        cert.get('epsilon_A')!='0.05' or cert.get('Q')!=1<<55 or
        cert.get('n_residues')!=n or cert.get('n_score_pairs')!=len(pairs)):
        raise ValueError('ABSTAIN: invalid certificate format or sealed input provenance')
    fixed,off,optional=check.exact_coord_boxes(coords,cert['epsilon_A'])
    if (len(fixed),len(off),len(optional))!=(972,6546,108):
        raise ValueError('ABSTAIN: wrong exact coordinate cutoff classification')
    if (cert.get('always_on'),cert.get('always_off'),cert.get('optional_edges'))!=(972,6546,108):
        raise ValueError('ABSTAIN: wrong declared contact classification')
    selected=cert['source_only_selected_indices']
    if selected!=SELECTED or cert['source_only_selected_edges']!=[list(optional[i]) for i in SELECTED]:
        raise ValueError('ABSTAIN: source-only contact split ordering changed')
    threshold=Fraction(cert['threshold'])
    if not (0<threshold<1) or cert['max_splits']!=6:
        raise ValueError('ABSTAIN: invalid positivity claim or depth limit')
    root=cert['root'];B=cert['tree_B'];Q=cert['Q']
    check.tree_checked(n,root,cert['tree_parent'],B,fixed)
    leaves=cert['leaves'];models=cert['models']
    if not isinstance(leaves,list) or not 1<=len(leaves)<=64:
        raise ValueError('ABSTAIN: invalid prefix leaves')
    if (not isinstance(models,list) or len(models)!=2*len(leaves) or
        cert['n_terminal_leaves']!=len(leaves) or
        cert['n_nodes_visited']!=2*len(leaves)-1):
        raise ValueError('ABSTAIN: incomplete cached graph witness set')
    prefixes=[(leaf['depth'],leaf['prefix_mask']) for leaf in leaves]
    coverage=verify_prefix_cover(prefixes,6)

    cache={}  # model ID -> (exact graph, per-pair rational numerators)
    graph_to_id={}
    used=set()
    def check_endpoint(mid,edges):
        if type(mid)!=int or not 0<=mid<len(models):
            raise ValueError('ABSTAIN: invalid endpoint model reference')
        key=tuple(sorted(edges))
        if len(set(key))!=len(key):raise ValueError('ABSTAIN: duplicate spring')
        if mid in cache:
            if cache[mid][0]!=key:raise ValueError('ABSTAIN: reused witness for different graph')
            used.add(mid)
            return cache[mid][1]
        if key in graph_to_id:
            raise ValueError('ABSTAIN: duplicate numerical witnesses for one canonical graph')
        scores=check.checked_scores(n,edges,pairs,models[mid],root,B,Q)
        cache[mid]=(key,scores)
        graph_to_id[key]=mid
        used.add(mid)
        return scores

    lower=Fraction(1);upper=Fraction(-1)
    exact_comparisons=0
    for leaf in leaves:
        d=leaf['depth'];mask=leaf['prefix_mask']
        forced_on={optional[selected[k]] for k in range(d) if (mask>>k)&1}
        forced_off={optional[selected[k]] for k in range(d) if not (mask>>k)&1}
        lo_graph=fixed+[e for e in optional if e not in forced_off]
        hi_graph=fixed+[e for e in optional if e in forced_on]
        lo=check_endpoint(leaf['all_unfixed_on_model'],lo_graph)
        hi=check_endpoint(leaf['all_unfixed_off_model'],hi_graph)
        lows=[x[0] for x in lo]
        highs=[x[1] for x in hi]
        if any(x>y for x,y in zip(lows,highs)):
            raise ValueError('ABSTAIN: graph response interval ordering violated')
        count=check.count_independent(lows,highs,obs)
        if leaf['target_concordance_counts']!=count:
            raise ValueError('ABSTAIN: fabricated leaf concordance')
        val=Fraction(count['lower']);up=Fraction(count['upper'])
        if val<threshold:
            raise ValueError('ABSTAIN: one terminal leaf does not certify requested threshold')
        lower=min(lower,val);upper=max(upper,up)
        exact_comparisons+=count['distinct_target_comparisons']
    if used!=set(range(len(models))):
        raise ValueError('ABSTAIN: unused/uncounted graph witnesses')
    if (cert['lower_concordance']!=str(lower) or
        cert['upper_concordance']!=str(upper) or lower<=0):
        raise ValueError('ABSTAIN: aggregate lower/upper claim not supported')
    result={'status':'VERIFIED','lower':str(lower),'upper':str(upper),
            'threshold':str(threshold),'leaf_cases':len(leaves),
            'endpoint_graphs':len(models),'complete_masks':len(coverage),
            'exact_pair_comparisons':exact_comparisons,
            'sealed_source':check.EXPECTED_STAGE1,
            'sealed_target':check.EXPECTED_STAGE2,
            'certificate_sha256':sha256(raw).hexdigest()}
    print('INDEPENDENT PREFIX-TREE CERTIFICATE VERIFIED',json.dumps(result,sort_keys=True))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('stage1');p.add_argument('stage2');p.add_argument('certificate')
    args=p.parse_args()
    verify(args.stage1,args.stage2,args.certificate)