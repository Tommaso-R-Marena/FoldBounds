"""UNTRUSTED numerical search for a compact complete GNM contact-world proof.

Outputs only a *proposal*. Always run replay_1d3z_complete_worlds.py
afterward: no numerical comparison in this file is an accepted certificate.
No NMR mobility labels, held-out protein data, or target optimization used.

Regenerate:
    python -m pip install numpy sympy
    python experiments/propose_1d3z_complete_worlds.py
    python experiments/replay_1d3z_complete_worlds.py
"""
from __future__ import annotations
import argparse,json
from itertools import combinations
from pathlib import Path
from fractions import Fraction
import numpy as np

from foldbounds.gnm_joint_rank import classify_contact_boxes
from foldbounds.gnm_reduced_woodbury import ReducedWorldCache,exact_joint_sign

if __package__:
    from .real_1d3z_joint_rank import DEFAULT_FIXTURE,load_ca
else:
    from real_1d3z_joint_rank import DEFAULT_FIXTURE,load_ca

ROOT=Path(__file__).resolve().parents[1]
PROPOSAL=ROOT/'research/gnm_1d3z_complete_worlds.json'
PRIOR=((41,49),(37,40),(5,11))


def make_proposal(fixture:Path=DEFAULT_FIXTURE):
    coords=load_ca(fixture);n=len(coords)
    certain,possible=classify_contact_boxes(coords,Fraction(8),Fraction(1,20))
    opt=tuple(sorted(set(possible)-set(certain)))
    if n!=69 or len(opt)!=15:raise ValueError('expected frozen 1D3Z 69×15 graph')

    def numeric_inverse(edges):
        mat=np.ones((n,n))/n
        for i,j in edges:
            mat[i,i]+=1;mat[j,j]+=1;mat[i,j]-=1;mat[j,i]-=1
        return np.linalg.inv(mat)-np.ones((n,n))/n

    def update(C,edge,adding):
        i,j=edge;v=C[:,i]-C[:,j];r=v[i]-v[j]
        den=(1+r if adding else 1-r)
        if den<=0:raise ValueError('non-PD optional contact update')
        return C+((-1 if adding else 1)/den)*np.outer(v,v)

    # Numerical/full 2^15-world explorer. NEVER accepted by itself.
    n_worlds=1<<len(opt)
    D=np.empty((n_worlds,n),dtype=np.float64)
    lo=numeric_inverse(certain);hi=numeric_inverse(possible)
    C=lo.copy()
    for t in range(n_worlds):
        gray=t^(t>>1)
        if t:
            changed=(gray^previous).bit_length()-1
            C=update(C,opt[changed],bool(gray&(1<<changed)))
        D[gray]=np.diag(C)
        previous=gray
    pairs=list(combinations(range(n),2))
    all_stable=[];flips=[]
    for i,j in pairs:
        difference=D[:,i]-D[:,j]
        low=float(difference.min());high=float(difference.max())
        if low>1e-8:all_stable.append([i,j,1])
        elif high< -1e-8:all_stable.append([i,j,-1])
        elif low< -1e-8 and high>1e-8:flips.append([i,j])
        else:raise ArithmeticError('inconclusive floating proposal; abstaining')

    # Re-certify prior signs using exact small-rank inverses, not floats.
    reduced=ReducedWorldCache(n,certain,possible)
    prior_idx=tuple(opt.index(e) for e in PRIOR)
    selected=sum(1<<k for k in prior_idx);full=n_worlds-1
    baseline=None
    for bits in range(1<<len(PRIOR)):
        lower=sum(1<<k for p,k in enumerate(prior_idx) if bits&(1<<p))
        upper=full^(selected^lower)
        L=reduced.get(lower);H=reduced.get(upper)
        signs={(i,j):s for i,j in pairs
               if (s:=exact_joint_sign(L,H,i,j))}
        baseline=(signs if baseline is None else
                  {p:s for p,s in baseline.items() if signs.get(p)==s})
    assert baseline is not None
    assert all((i,j) not in baseline or baseline[i,j]==s for i,j,s in all_stable)
    extras=[(i,j,s) for i,j,s in all_stable if (i,j) not in baseline]
    if len(all_stable)!=1980 or len(flips)!=366 or len(extras)!=33:
        raise ArithmeticError('unexpected preliminary model classification')

    # A short counterexample set covers both signs for all variable pairs.
    ai=np.array([i for i,j in flips]);aj=np.array([j for i,j in flips])
    diff=D[:,ai]-D[:,aj]
    pos=diff>1e-7;neg=diff< -1e-7
    needp=np.ones(len(flips),dtype=bool);needn=needp.copy();masks=[]
    while needp.any() or needn.any():
        scores=np.count_nonzero(pos[:,needp],axis=1)+np.count_nonzero(neg[:,needn],axis=1)
        mask=int(np.argmax(scores))
        if scores[mask]==0:raise ArithmeticError('unable to find counterexamples')
        masks.append(mask)
        needp &= ~pos[mask];needn &= ~neg[mask]

    # A dynamic rank-two bound-driven tree is a proposal. The independent
    # exact checker reconstructs EVERY leaf and all branch coverage.
    ii=np.array([i for i,j,s in extras]);jj=np.array([j for i,j,s in extras]);ss=np.array([s for i,j,s in extras])
    nodes=0;leaves=0;max_depth=0
    def unresolved(L,H,active):
        if not len(active):return active
        i=ii[active];j=jj[active];sgn=ss[active]
        a=L[i,i]-H[i,i];b=L[j,j]-H[j,j];c=L[i,j]-H[i,j]
        s=2*(H[i,i]-H[j,j])+a-b
        delta=(a+b)**2-4*c*c
        certified=(s*sgn>1e-9)&(s*s>delta+1e-9)
        return active[~certified]
    def branch(L,H,remaining,active,depth):
        nonlocal nodes,leaves,max_depth
        nodes+=1;max_depth=max(max_depth,depth)
        active=unresolved(L,H,active)
        if not len(active):
            leaves+=1
            return {'leaf':True}
        if not remaining:raise ArithmeticError('full leaf remains numerically unresolved')
        best=None
        for edge in remaining:
            L1=update(L,opt[edge],True)
            H0=update(H,opt[edge],False)
            a=unresolved(L1,H,active);b=unresolved(L,H0,active)
            rank=(len(a)+len(b),max(len(a),len(b)),edge)
            if best is None or rank<best[0]:best=(rank,edge,L1,H0)
        _,index,L1,H0=best
        rest=[e for e in remaining if e!=index]
        return {'split':index,'zero':branch(L,H0,rest,active,depth+1),
                'one':branch(L1,H,rest,active,depth+1)}
    tree=branch(lo,hi,list(range(len(opt))),np.arange(len(extras)),0)
    return {'schema':'foldbounds-gnm-full-contact-worlds-v1',
            'protein':'PDB 1D3Z MODEL 1 CA residues 2-70',
            'coordinate_error_angstrom':'1/20',
            'radius_angstrom':'8',
            'optional_edges':[list(e) for e in opt],
            'extra_stable_orders':[list(x) for x in extras],
            'changing_pairs':flips,'witness_masks':masks,'tree':tree,
            'proposer_statistics':{'full_boolean_worlds':n_worlds,
                 'stably_signed':len(all_stable),'variable':len(flips),
                 'additional_over_prior_3split':len(extras),
                 'tree_nodes':nodes,'tree_leaves':leaves,
                 'tree_max_depth':max_depth,'witness_masks':len(masks)}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture',type=Path,default=DEFAULT_FIXTURE)
    parser.add_argument('--output',type=Path,default=PROPOSAL)
    args=parser.parse_args()
    result=make_proposal(args.fixture)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,separators=(',',':'))+'\n',encoding='utf-8')
    print('UNTRUSTED proposal saved at',args.output)
    print(json.dumps(result['proposer_statistics'],sort_keys=True))
    print('Run the independent exact verifier before trusting ANY claims.')

if __name__=='__main__':main()
