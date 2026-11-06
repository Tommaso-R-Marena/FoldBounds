#!/usr/bin/env python3
"""Find exact, SOURCE-ONLY Farkas witnesses for impossible GNM contact states.

The untrusted optimizer is SciPy. The producer certifies every reported
inequality again with fractions; a completely separate verifier rechecks it
without importing SciPy, optimizer, or this module. Never opens an NMR target.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from gnm_exact_inverse_producer import archive_source,exact_edges

EPS=F(1,20); R2=F(100); SOURCE='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
FORMAT='FOLDBOUNDS_REFERENCE_ONLY_GEOMETRIC_NOGOOD_V1'

def halfspace(edge,is_on,xyz):
    i,j=edge
    r=[xyz[i][k]-xyz[j][k] for k in range(3)]
    diff_sq=sum(z*z for z in r)
    sign=1 if is_on else -1
    coefficients={(i,k):2*sign*r[k] for k in range(3)}
    coefficients.update({(j,k):-2*sign*r[k] for k in range(3)})
    # ON: |r+delta_i-delta_j|^2 <= R2 -> linear <= R2-|r|^2
    # OFF: >R2 and residual quadratic <=12eps² -> -linear <= |r|²-R2+12eps²
    limit=R2-diff_sq if is_on else diff_sq-R2+12*EPS*EPS
    return coefficients,limit

def verify_witness(edges,flags,weights,xyz):
    if not (len(edges)==len(flags)==len(weights) and all(w>=0 for w in weights)
            and any(w>0 for w in weights)):raise ValueError('invalid multipliers')
    source=[halfspace(e,on,xyz) for e,on in zip(edges,flags)]
    variables=sorted(set().union(*(set(a) for a,b in source)))
    lhs=-EPS*sum(abs(sum(weights[j]*source[j][0].get(v,F(0)) for j in range(len(edges)))) for v in variables)
    rhs=sum(weights[j]*source[j][1] for j in range(len(edges)))
    gap=lhs-rhs
    if gap<=0:raise ValueError('Farkas witness is not a contradiction')
    return gap

def pair_witness(e,f,on1,on2,xyz):
    src=[halfspace(e,on1,xyz),halfspace(f,on2,xyz)]
    keys=sorted(set(src[0][0])|set(src[1][0]))
    coeff1=[src[0][0].get(k,F(0)) for k in keys]
    coeff2=[src[1][0].get(k,F(0)) for k in keys]
    candidates=[F(0)]+[-a/b for a,b in zip(coeff1,coeff2) if b and -a/b>0]
    # Piecewise-affine cone separation; maximum occurs at one of these exact breakpoints.
    weights=[[F(1),t] for t in candidates]+[[F(0),F(1)]]
    return [(ws,verify_witness([e,f],[on1,on2],ws,xyz)) for ws in weights
            if candidate_is_valid([e,f],[on1,on2],ws,xyz)]

def candidate_is_valid(edges,flags,weights,xyz):
    try:verify_witness(edges,flags,weights,xyz);return True
    except ValueError:return False

def higher_order_witness(edges,flags,xyz):
    src=[halfspace(e,on,xyz) for e,on in zip(edges,flags)]
    variables=sorted(set().union(*(set(a) for a,b in src)))
    n=len(edges);m=len(variables)
    A=np.asarray([[float(a.get(k,0)) for k in variables] for a,b in src])
    b=np.asarray([float(rhs) for a,rhs in src])
    # Minimize sum lambda_i*b_i + eps*sum_j |sum_i lambda_i a_ij|
    # for lambda>=0 and sum lambda=1. A negative solution is a proposed
    # separating hyperplane. Only the exact check is authoritative.
    obj=np.r_[b,np.full(m,float(EPS))]
    ub=np.zeros((m*2,n+m))
    ub[:m,:n]=A.T;ub[:m,n:]=-np.eye(m)
    ub[m:,:n]=-A.T;ub[m:,n:]=-np.eye(m)
    equal=np.r_[np.ones(n),np.zeros(m)].reshape(1,-1)
    result=linprog(obj,A_ub=ub,b_ub=np.zeros(2*m),A_eq=equal,b_eq=[1],
                   bounds=[(0,None)]*(n+m),method='highs')
    if not result.success:raise ValueError('optimizer failed: '+result.message)
    if result.fun>=-1e-9:return None
    weights=[F(float(v)).limit_denominator(1_000_000) for v in result.x[:n]]
    if not candidate_is_valid(edges,flags,weights,xyz):return None
    return weights

def exact_component_count(n_edges,clauses):
    parent=list(range(n_edges))
    def root(x):
        while parent[x]!=x: parent[x]=parent[parent[x]];x=parent[x]
        return x
    for indices,_ in clauses:
        for v in indices[1:]:parent[root(v)]=root(indices[0])
    grouped={}
    for i in range(n_edges):grouped.setdefault(root(i),[]).append(i)
    product=1;components=[]
    for vertices in sorted(grouped.values(),key=lambda v:(-len(v),v[0])):
        contained=[(ids,flags) for ids,flags in clauses if ids[0] in vertices]
        if len(vertices)>20:raise ValueError('component too large; use #SAT')
        valid=0
        for mask in range(1<<len(vertices)):
            choice={node:bool(mask>>j&1) for j,node in enumerate(vertices)}
            if not any(all(choice[e]==on for e,on in zip(ids,flags)) for ids,flags in contained):valid+=1
        components.append({'edge_ids':vertices,'forbidden_clauses':len(contained),
                           'valid_assignments':valid,'total_assignments':1<<len(vertices)})
        product*=valid
    return product,components

def produce(stage1,output):
    start=time.time()
    coords,pairs,_,npzsha=archive_source(stage1)
    xyz=[[F.from_float(float(v)) for v in row] for row in coords]
    certain=set(exact_edges(coords,'certain'));possible=set(exact_edges(coords,'possible'))
    edges=sorted(possible-certain);assert len(edges)==108
    clauses=[];has_pair=set()
    for i,j in itertools.combinations(range(108),2):
        if not (set(edges[i])&set(edges[j])):continue
        for state in itertools.product((False,True),repeat=2):
            cert=pair_witness(edges[i],edges[j],*state,xyz)
            if not cert:continue
            weights,_=max(cert,key=lambda entry:entry[1])
            gap=verify_witness([edges[i],edges[j]],state,weights,xyz)
            clauses.append({'edge_ids':[i,j],'forbidden_on_bits':[int(x) for x in state],
                            'lambda':[str(w) for w in weights],'gap':str(gap)})
            has_pair.add(((i,j),state))
    n_pair=len(clauses)
    for ids in itertools.combinations(range(108),3):
        if len(set(edges[ids[0]])|set(edges[ids[1]])|set(edges[ids[2]]))>4:continue
        for state in itertools.product((False,True),repeat=3):
            if any(((ids[a],ids[b]),(state[a],state[b])) in has_pair
                   for a,b in ((0,1),(0,2),(1,2))):continue
            weights=higher_order_witness([edges[i] for i in ids],state,xyz)
            if weights is None:continue
            gap=verify_witness([edges[i] for i in ids],state,weights,xyz)
            clauses.append({'edge_ids':list(ids),'forbidden_on_bits':[int(x) for x in state],
                            'lambda':[str(w) for w in weights],'gap':str(gap)})
    clauses.sort(key=lambda c:(len(c['edge_ids']),c['edge_ids'],c['forbidden_on_bits']))
    allowed,comps=exact_component_count(len(edges),[(c['edge_ids'],[bool(v) for v in c['forbidden_on_bits']]) for c in clauses])
    total=1<<len(edges)
    output_data={'format':FORMAT,'source_stage1_sha256':SOURCE,'source_npz_sha256':npzsha,
                 'box_epsilon_A':'1/20','contact_cutoff_squared_A':'100',
                 'n_optional_edges':108,'optional_edges':[list(e) for e in edges],
                 'n_pairs_impossible':n_pair,'n_triples_impossible':len(clauses)-n_pair,
                 'forbidden_patterns':clauses,'components':comps,
                 'independent_binary_assignments':str(total),
                 'not_yet_excluded_assignments':str(allowed),
                 'excluded_assignments':str(total-allowed),
                 'excluded_fraction':str(F(total-allowed,total)),
                 'scope':'Necessary contact-geometry constraints; an upper bound on realizable network count; not a distribution or biological prediction.'}
    encoded=(json.dumps(output_data,sort_keys=True,separators=(',',':'))+'\n').encode()
    Path(output).write_bytes(encoded)
    return {'source_only':True,'pair_conflicts':n_pair,'new_triple_conflicts':len(clauses)-n_pair,
            'allowed_network_assignments':allowed,'possible_binary_assignments':total,
            'at_least_fraction_ruled_out':str(F(total-allowed,total)),
            'certificate_sha256':hashlib.sha256(encoded).hexdigest(),'bytes':len(encoded),'elapsed_s':round(time.time()-start,3)}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage1');ap.add_argument('output');args=ap.parse_args()
    print(json.dumps(produce(args.stage1,args.output),indent=2))