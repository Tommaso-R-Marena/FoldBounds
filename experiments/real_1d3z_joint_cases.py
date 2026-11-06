"""Source-only, exhaustive exact contact case splits on public PDB 1D3Z.

Recomputes every endpoint integer inverse and checks all B*N=dI identities.
No NMR conformer spread or target labels are ever read by this script.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import json
from foldbounds.gnm_joint_rank import classify_contact_boxes
from foldbounds.gnm_joint_cases import certify_split_cases
# Resolve the sibling module under both direct script and package execution.
# Direct 'python experiments/script.py' puts experiments/ (not repo root)
# on sys.path; 'python -m experiments.script' uses the package import.
if __package__:
    from .real_1d3z_joint_rank import load_ca, DEFAULT_FIXTURE
else:
    from real_1d3z_joint_rank import load_ca, DEFAULT_FIXTURE


def evaluate(epsilon: Fraction=Fraction(1,20), count:int=3, *, fixed_two:bool=False):
    coords=load_ca(DEFAULT_FIXTURE)
    certain,possible=classify_contact_boxes(coords,Fraction(8),epsilon)
    if fixed_two:
        if count!=2:
            raise ValueError('fixed_two benchmark requires selected_count=2')
        # PDB residue pairs 43-51 and 39-42, selected by source-only
        # exact rank-coverage greedy search without observing NMR labels.
        override=((41,49),(37,40))
    else:
        override=None
    r=certify_split_cases(len(coords),certain,possible,selected_count=count,
                          selected_override=override)
    return {'pdb':'1D3Z','coordinate_box_angstrom':str(epsilon),
            'optional_nonlocal_edges':r.optional_edge_count,
            'selected_source_only_optional_edges':[[u+2,v+2] for u,v in r.selected_edges],
            'selected_count':len(r.selected_edges),
            'exhaustive_cases':len(r.conclusions),
            'verified_integer_endpoint_matrices':r.inverse_witness_count,
            'single_case_joint_certified':r.baseline_pairs,
            'all_case_joint_certified':r.case_certified_pairs,
            'newly_certified':r.added_pairs,
            'all_residue_pairs':len(coords)*(len(coords)-1)//2}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epsilon',default='0.05')
    p.add_argument('--selected-count',type=int,default=3)
    p.add_argument('--fixed-two',action='store_true',help='exact 4-case source-only selected baseline')
    a=p.parse_args()
    print(json.dumps(evaluate(Fraction(a.epsilon),a.selected_count,fixed_two=a.fixed_two),indent=2))


if __name__=='__main__': main()
