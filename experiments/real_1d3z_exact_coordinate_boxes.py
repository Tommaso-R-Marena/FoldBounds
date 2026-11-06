"""Reproduce the tighter exact ±epsilon Cartesian-box GNM graph atlas.

Public PDB 1D3Z MODEL 1 source CA positions; no NMR target labels.
Not a biologically calibrated protein motion prediction.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from fractions import Fraction
from foldbounds.gnm_joint_rank import classify_contact_boxes
from foldbounds.gnm_exact_coordinate_boxes import (
    exact_coordinate_box_graphs,audit_exact_box_gnm_ranks
)
if __package__:
    from .real_1d3z_joint_rank import load_ca,DEFAULT_FIXTURE
else:
    from real_1d3z_joint_rank import load_ca,DEFAULT_FIXTURE


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixture',type=Path,default=DEFAULT_FIXTURE)
    p.add_argument('--epsilon',default='0.05')
    p.add_argument('--max-worlds',type=int,default=8192)
    args=p.parse_args()
    coords=load_ca(args.fixture)
    eps=Fraction(args.epsilon)
    tight=exact_coordinate_box_graphs(coords,Fraction(8),eps)
    coarse_certain,coarse_possible=classify_contact_boxes(coords,Fraction(8),eps)
    if not set(coarse_certain)<=set(tight.certainly_on) or not set(tight.possibly_on)<=set(coarse_possible):
        raise ArithmeticError('tight box geometry not nested inside conservative graph envelope')
    audit=audit_exact_box_gnm_ranks(tight,max_worlds=args.max_worlds)
    excluded=sorted(set(coarse_possible)-set(tight.possibly_on))
    print(json.dumps({
        'source':'public PDB 1D3Z model 1; source-only, uncalibrated geometry',
        'epsilon_angstrom':str(eps),
        'coarse_optional_count':len(set(coarse_possible)-set(coarse_certain)),
        'exact_pairwise_optional_count':audit.optional_contacts,
        'excluded_impossible_contacts_pdb':[[i+2,j+2] for i,j in excluded],
        'exact_boolean_worlds_checked':audit.exact_graph_worlds,
        'invariant_residue_pairs':audit.globally_invariant,
        'opposite_sign_pairs':audit.opposite_sign_witnessed,
        'tie_only':audit.noninvariant_with_ties_only,
        'always_tied':audit.always_tied,
        'small_exact_inverses_verified':audit.reduced_inverses_checked,
    },indent=2,sort_keys=True))
    if eps==Fraction(1,20):
        assert (audit.optional_contacts,audit.exact_graph_worlds,
                audit.globally_invariant,audit.opposite_sign_witnessed,
                audit.noninvariant_with_ties_only,audit.always_tied)==(
                 12,4096,2096,250,0,0)

if __name__=='__main__':main()
