#!/usr/bin/env python3
"""Uniform all-coordinate-box uncertainty bounds for GNM graph resistance.

For each reference coordinate, allow INDEPENDENT Cartesian deviation +/-epsilon.
Exact rational interval geometry classifies every potential C-alpha contact:
ALWAYS_ON (even at maximal separation), ALWAYS_OFF (even at minimal
separation), or OPTIONAL. The true graph in every coordinate world lies
between G_on and G_off. Rayleigh monotonicity of electrical resistance
then bounds each pair response between the scores of G_on (minimum) and
G_off (maximum), WITHOUT enumerating 2^k contact-switch combinations.

Exact score intervals for both endpoint graphs are built with untrusted
floating inverse proposals but verified by all-integer residual arithmetic.
This is topology-only structural robustness, NOT biological validation, and
it does not bound ANM's directional changes within a fixed contact graph.
"""
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import numpy as np
import gnm_exact_certificates as gc


def classify_coord_box(coords, epsilon, cutoff_sq=Fraction(100)):
    if epsilon < 0:raise ValueError('negative epsilon')
    x=np.asarray(coords,np.float64)
    if not np.isfinite(x).all():raise ValueError('nonfinite coordinates')
    X=[[Fraction.from_float(float(q)) for q in row] for row in x]
    always_on=[];always_off=[];maybe=[]
    for i in range(len(x)):
        for j in range(i+1,len(x)):
            delta=[abs(X[i][k]-X[j][k]) for k in range(3)]
            min_d2=sum(max(Fraction(0),d-2*epsilon)**2 for d in delta)
            max_d2=sum((d+2*epsilon)**2 for d in delta)
            if min_d2==0:raise ValueError('ABSTAIN: coincident residues possible in uncertainty box')
            if max_d2<=cutoff_sq:always_on.append((i,j))
            elif min_d2>cutoff_sq:always_off.append((i,j))
            else:maybe.append((i,j))
    return always_on,always_off,maybe


def _endpoint_score(edges,n,pairs,Q):
    neighbors=[set() for _ in range(n)]
    for i,j in edges:
        neighbors[i].add(j);neighbors[j].add(i)
    root,B,depth,parent=gc.choose_tree(neighbors)
    A,idx=gc.grounded_laplacian(edges,n,root)
    N=gc.untrusted_inverse_proposal(A,Q)
    R=gc.verify_integer_witness(A,N,Q)
    bounds=[gc.score_numerators(p,N,R,Q,B,idx) for p in pairs]
    return bounds,{'root':root,'B':B,'tree_parent':parent,'n_edges':len(edges),
                   'integer_inverse_proposal':N,'denominator_Q':Q}


def robust_scores(coords,pairs_i,pairs_j,epsilon,Q=gc.Q_DEFAULT):
    x,pairs=gc._validate_inputs(coords,pairs_i,pairs_j)
    fixed,off,optional=classify_coord_box(x,epsilon)
    # Superset of *all physically realizable* networks. Optional edges need
    # not be jointly realizable; extrema remain mathematically conservative.
    lower,_meta_lo=_endpoint_score(fixed+optional,len(x),pairs,Q)
    upper,_meta_hi=_endpoint_score(fixed,len(x),pairs,Q)
    hull=[]
    for i in range(len(pairs)):
        lo=lower[i][0]
        hi=upper[i][1]
        if lo>hi:
            raise AssertionError('Rayleigh extreme graph monotonicity violated')
        hull.append((lo,hi))
    # All inequalities share common denominator Q^2.
    old=gc.verify_and_score(x,pairs_i,pairs_j,Q,
        gc.untrusted_inverse_proposal(
            gc.grounded_laplacian(gc.exact_contact_graph(x)[0],len(x),
                                  gc.choose_tree(gc.exact_contact_graph(x)[1])[0])[0],Q))
    nominal_order=old['ordered_indices']
    lo=[v[0] for v in hull];hi=[v[1] for v in hull];n=len(pairs)
    pairs_certified=0
    for a in range(n):
        va=nominal_order[a]
        pairs_certified+=sum(hi[va]<lo[nominal_order[b]] for b in range(a+1,n))
    adjacent=sum(hi[nominal_order[i]]<lo[nominal_order[i+1]] for i in range(n-1))
    k=max(1,int(np.ceil(.10*n)))
    # For descending (high response) ranking, a guaranteed top-k item must
    # outrank enough entries by strict interval separation.
    definitely_above=[sum(j!=i and lo[j]>hi[i] for j in range(n)) for i in range(n)]
    definitely_below=[sum(j!=i and hi[j]<lo[i] for j in range(n)) for i in range(n)]
    guaranteed_top=[i for i in range(n) if n-definitely_below[i]<=k]
    possible_top=[i for i in range(n) if 1+definitely_above[i]<=k]
    return dict(epsilon=str(epsilon),n_pairs=n,always_on=len(fixed),
                always_off=len(off),optional_edges=len(optional),
                certified_pair_order=pairs_certified,total_pair_order=n*(n-1)//2,
                certified_adjacent=adjacent,total_adjacent=n-1,
                guaranteed_top_k=len(guaranteed_top),possible_top_k=len(possible_top),
                top_k=k,root_allon=_meta_lo['root'],root_alloff=_meta_hi['root'],
                lower_score_numerators=lo,upper_score_numerators=hi,
                denominator=Q*Q,
                endpoint_meta={'all_on':_meta_lo,'all_off':_meta_hi},
                note='Rayleigh monotonicity for every possible contact subset; structural boxes only')


def main():
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('stage1_zip');p.add_argument('output_json');p.add_argument('--eps',nargs='+',default=['0.0005','0.005','0.05'])
    x=p.parse_args()
    arrays,_=gc.open_sealed_stage1(x.stage1_zip)
    out={'format':'GNM_TOPOLOGY_BOX_HULL_V1',
         'source_archive_sha256':gc.STAGE1_SHA,
         'scores':[]}
    for eps in x.eps:
        print('VERIFY coordinate epsilon',eps,flush=True)
        score=robust_scores(arrays['reference_coords'],arrays['pair_i'],arrays['pair_j'],Fraction(eps))
        out['scores'].append(score)
        print('CERTIFIED:',score['certified_pair_order'],'/',score['total_pair_order'],
              'guaranteed top-k',score['guaranteed_top_k'],'/',score['top_k'],flush=True)
    path=Path(x.output_json);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print('saved',path)


if __name__=='__main__':main()