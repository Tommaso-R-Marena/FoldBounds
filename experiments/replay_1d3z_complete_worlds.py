"""Independent exact verification of a proposed COMPLETE 1D3Z GNM rank atlas.

No floating-point calculations and no PDB motion labels in acceptance.
Reproduces all 1980 invariant / 366 reversible Boolean-graph rankings.
"""
from __future__ import annotations
import argparse,json,sys
from fractions import Fraction
from pathlib import Path
from foldbounds.gnm_joint_rank import classify_contact_boxes
from foldbounds.gnm_complete_worlds import verify_complete_contact_worlds


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixture',type=Path,default=Path(__file__).resolve().parents[1]/'tests/fixtures/1d3z_model1_ca_milli.csv')
    p.add_argument('--proposal',type=Path,default=Path(__file__).resolve().parents[1]/'research/gnm_1d3z_complete_worlds.json')
    p.add_argument('--max-tree-nodes',type=int,default=10000)
    args=p.parse_args()
    # Reuse the separately tested PDB coordinate fixture parser from #76.
    if __package__:
        from .real_1d3z_joint_rank import load_ca
    else:
        from real_1d3z_joint_rank import load_ca
    coords=load_ca(args.fixture)
    certain,possible=classify_contact_boxes(coords,Fraction(8),Fraction(1,20))
    certificate=json.loads(args.proposal.read_text(encoding='utf-8'))
    audit=verify_complete_contact_worlds(len(coords),certain,possible,certificate,
              selected_prior=((41,49),(37,40),(5,11)),
              max_tree_nodes=args.max_tree_nodes)
    print(json.dumps(audit.__dict__,indent=2,sort_keys=True))
    if (audit.residues,audit.optional_contacts,audit.invariant_pairs,
        audit.reversing_pairs,audit.extra_certified)!=(69,15,1980,366,33):
        raise ArithmeticError('unexpected pinned dataset result')

if __name__=='__main__':main()
