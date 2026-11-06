"""Opt-in source-only greedy contact selection; costly but deterministic.

It reads the 1D3Z MODEL-1 C-alpha coordinates and exact network
certificate coverage only. It NEVER reads or tunes to any NMR target.
Selected source-only pairs for depth two are PDB (43,51) and (39,42).
"""
import argparse
from fractions import Fraction
# Resolve the sibling module under both direct script and package execution.
# Direct 'python experiments/script.py' puts experiments/ (not repo root)
# on sys.path; 'python -m experiments.script' uses the package import.
if __package__:
    from .real_1d3z_joint_rank import load_ca, DEFAULT_FIXTURE
else:
    from real_1d3z_joint_rank import load_ca, DEFAULT_FIXTURE
from foldbounds.gnm_joint_rank import classify_contact_boxes
from foldbounds.gnm_joint_cases import certify_split_cases
from foldbounds.gnm_integer_updates import ExactEndpointCache
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--depth',type=int,default=2)
args=parser.parse_args()
if not 1<=args.depth<=3:raise ValueError('source-only exhaustive heuristic depth must be 1..3')
coords=load_ca(DEFAULT_FIXTURE)
certain,possible=classify_contact_boxes(coords,Fraction(8),Fraction(1,20))
optional=sorted(set(possible)-set(certain))
cache=ExactEndpointCache(len(coords),certain)
chosen=[]
for depth in range(1,args.depth+1):
 results=[]
 for e in optional:
  if e in chosen:continue
  selected=tuple(chosen+[e])
  result=certify_split_cases(len(coords),certain,possible,selected_count=depth,selected_override=selected,update_cache=cache)
  results.append((result.case_certified_pairs,e,result))
 winner=max(results,key=lambda x:(x[0],tuple(-v for v in x[1])))
 chosen.append(winner[1])
 print('cached inverses',cache.count,'depth',depth,'certified',winner[0],'selected',chosen,'other candidates',[(x[0],x[1]) for x in sorted(results,key=lambda x:-x[0])[:4]],flush=True)
