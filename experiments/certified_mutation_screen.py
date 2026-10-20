"""Full exact mutation, adversarial uncertainty and nonlocal epistasis screen.

Run: PYTHONPATH=src python experiments/certified_mutation_screen.py
Only standard-library arithmetic; finite 3D lattice state space (8 residues).
"""
from __future__ import annotations
import json
from fractions import Fraction as F
from math import comb
from foldbounds.mutation_switch import (
    _sequence, compact_spectrum, geometric_fingerprints, nominal_probability,
    certify_mutation_switch,verify_mutation_switch,
    certify_uniform_mutation_landscape,verify_uniform_mutation_landscape)
from foldbounds.ensemble_uncertainty import lower_probability,upper_probability


def exact(v):return str(v)


def run():
    q0,q1,b0,b1=F(8),F(12),F(3,4),F(5,4)
    flo,fhi=F(5,6),F(6,5)
    x=certify_mutation_switch('PHPHPHPH','HHPHPHPH',q0,q1,b0,b1,flo,fhi)
    if not verify_mutation_switch(x):raise AssertionError('mutation failed independent checker')
    landscape=certify_uniform_mutation_landscape(q0,q1,b0,b1,flo,fhi)
    if not verify_uniform_mutation_landscape(landscape):raise AssertionError('landscape failed independent checker')
    assert len(landscape.robust_switch_edges)==36
    assert len(geometric_fingerprints())==171
    def rev(mask):return int(f'{mask:08b}'[::-1],2)
    assert all(compact_spectrum(m)==compact_spectrum(rev(m)) for m in range(256))
    orbit_keys={min((a,b,p),(rev(a),rev(b),9-p)) for a,b,p in landscape.robust_switch_edges}
    assert len(orbit_keys)==18
    # Full double-mutant cycles at q=8,b=1; geometry event is mutation invariant.
    pop=[nominal_probability(m,F(8),F(1)) for m in range(256)]
    cycles=0;even_gap=0;even_gap_epistasis=0;double_flip=0;double_flip_robust_10=0
    example=None
    for m in range(256):
        for i in range(8):
            if m&(1<<i):continue
            for j in range(i+1,8):
                if m&(1<<j):continue
                a=m|(1<<i);b=m|(1<<j);z=a|(1<<j)
                cycles+=1
                if (j-i)%2==0:
                    even_gap+=1
                    # Compare *odds* multiplicatively with exact integer arithmetic.
                    def odds(x):return x/(1-x)
                    if odds(pop[z])*odds(pop[m]) != odds(pop[a])*odds(pop[b]):
                        even_gap_epistasis+=1
                    if max(pop[m],pop[a],pop[b])<F(1,2)<pop[z]:
                        double_flip+=1
                        if all(upper_probability(pop[y],F(10,11),F(11,10))<F(1,2) for y in (m,a,b)) and lower_probability(pop[z],F(10,11),F(11,10))>F(1,2):
                            double_flip_robust_10+=1
                        if m==1 and i==3 and j==5:
                            example=(_sequence(m),_sequence(a),_sequence(b),_sequence(z),*[exact(pop[y]) for y in (m,a,b,z)])
    assert (cycles,even_gap,even_gap_epistasis,double_flip,double_flip_robust_10)==(1792,768,504,64,30)
    return {
      'model':'anchored 3D cubic-lattice 8-site binary HP chain (not atomistic)',
      'event':'at least two nonconsecutive geometric lattice-neighbor contacts; event invariant under mutations',
      'anchor_oriented_conformations':13565,
      'lossless_geometric_fingerprints':len(geometric_fingerprints()),
      'sequences':256, 'undirected_one_site_edges_checked':landscape.checked_edges,
      'robust_one_site_switches_parameter_rectangle':len(landscape.robust_switch_edges),
      'robust_switch_reversal_orbits':len(orbit_keys),
      'bounds':dict(contact_weight=[exact(q0),exact(q1)],bend_weight=[exact(b0),exact(b1)],
                    per_state_correction=[exact(flo),exact(fhi)]),
      'example':{
         'parent':x.parent.sequence,'mutant':x.mutant.sequence,'mutated_residue':x.mutated_residue,
         'nominal_parent_probability':[exact(x.parent.event_lower),exact(x.parent.event_upper)],
         'nominal_mutant_probability':[exact(x.mutant.event_lower),exact(x.mutant.event_upper)],
         'robust_parent_upper':exact(x.parent_robust_upper),
         'robust_mutant_lower':exact(x.mutant_robust_lower),
         'absolute_population_separation_lower_bound':exact(x.mutant_robust_lower-x.parent_robust_upper)
      },
      'all_robust_switch_edges': [dict(from_sequence=_sequence(a),to_sequence=_sequence(b),residue=p) for a,b,p in landscape.robust_switch_edges],
      'double_mutation_epistasis':{
          'total_2_site_mutational_squares':cycles,
          'same_parity_site_squares_no_possible_direct_contact':even_gap,
          'nonzero_exact_odds_epistasis_despite_no_direct_contact':even_gap_epistasis,
          'nominal_double_only_majority_crossings_in_same_parity_squares_at_q8_b1':double_flip,
          'same_crossings_robust_to_statewise_factors_10_over_11_to_11_over_10_at_fixed_q8_b1':double_flip_robust_10,
          'example_parent_single_single_double':example
      },
      'assurance_boundary':'Exact within enumerated toy model; Bernstein/weight proofs are mathematical but histogram construction Python-tested rather than Lean kernel verified; no physical prediction.'
    }


if __name__=='__main__':
    print(json.dumps(run(),indent=2,sort_keys=True))