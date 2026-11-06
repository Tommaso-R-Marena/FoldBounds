#!/usr/bin/env python3
"""Source-only dyadic precision/certification comparison (NO NMR target)."""
from pathlib import Path
import sys,json
import numpy as np
from fractions import Fraction
from gnm_exact_certificates import open_sealed_stage1,grounded_laplacian,verify_integer_witness
from gnm_topology_boxes import classify_coord_box

def experiment(path):
    a,_=open_sealed_stage1(path)
    coords=a['reference_coords']; edges,_,opt=classify_coord_box(coords,Fraction('.05'))
    n=len(coords);anchor=45;A,idx=grounded_laplacian(edges,n,anchor)
    matrix=np.asarray(A,dtype=np.float64)
    inv=np.linalg.inv(matrix)
    pairs=list(zip(map(int,a['pair_i']),map(int,a['pair_j'])))
    def check(Bn,Bd,Q,N,R):
        intervals=[]
        for i,j in pairs:
            mi,mj=idx.get(i),idx.get(j)
            u=[(N[k][mi] if mi is not None else 0)-(N[k][mj] if mj is not None else 0) for k in range(n-1)]
            r=[(R[k][mi] if mi is not None else 0)-(R[k][mj] if mj is not None else 0) for k in range(n-1)]
            bx=(u[mi] if mi is not None else 0)-(u[mj] if mj is not None else 0)
            lo=Q*bx+sum(x*y for x,y in zip(u,r))
            hi=lo*Bd+Bn*sum(x*x for x in r)
            intervals.append((lo*Bd,hi))
        order=sorted(range(len(pairs)),key=lambda i:(intervals[i][0],i))
        certified=sum(intervals[order[i]][1]<intervals[order[i+1]][0] for i in range(len(pairs)-1))
        maxwidth=max((hi-lo)/(Bd*Q**2) for lo,hi in intervals)
        return certified,maxwidth
    out=[]
    for bits in (10,12,14,16,18,20,22,24,26,28,30,32,36,40,44,48,52,55):
        Q=1<<bits;N=[[int(round(float(x)*Q)) for x in row] for row in inv]
        R=verify_integer_witness(A,N,Q)
        old=check(262,1,Q,N,R)
        new=check(1000000,169821,Q,N,R)
        out.append({'q_bits':bits,'old_adjacent_orders':old[0],'new_adjacent_orders':new[0],
                    'old_max_score_width':old[1],'new_max_score_width':new[1]})
        print('bits',bits,'old_rank',old[0],'new_rank',new[0],'maxWidths',old[1],new[1],flush=True)
    return out
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser()
    pa.add_argument('source_stage1')
    pa.add_argument('output_json')
    args=pa.parse_args()
    rows=experiment(args.source_stage1)
    result={'format':'GNM_SOURCE_ONLY_PRECISION_ABLATION_V1',
       'source_sha256':'ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff',
       'coercivity_old':262, 'coercivity_spectral':'1000000/169821',
       'pairs':2614,'n_certain_on_edges':972,'n_optional_edges':108,
       'results':rows,
       'interpretation':'Adjacency counts are over the sorted lower-bound interval order at each resolution; all intervals are exact-rational residual certificates. No target NMR data accessed.'}
    Path(args.output_json).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('WROTE',args.output_json)