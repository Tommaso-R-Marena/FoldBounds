"""Reproduce full robust same-parity 3D HP mutation landscapes exactly."""
import json
from fractions import Fraction as F
from foldbounds.composition_swap import (
    DEFAULT,certify_landscape,verify_landscape,certify_swap,
    certify_error_budget,sequence)

def run():
    out={'model':'anchored 3D cubic HP model, not real protein folding',
         'parameter_box':[str(x) for x in DEFAULT],'cases':[]}
    for n,parent,mutant,states in (
        (8,'PPPHPHHH','HPPHPHPH',13565),
        (10,'PPPHPHHHPH','HPPHPHPHPH',308981)):
        landscape=certify_landscape(n)
        assert verify_landscape(landscape)
        c=certify_swap(parent,mutant)
        threshold=certify_error_budget(parent,mutant)
        expected=(768,12,6) if n==8 else (5120,4,2)
        assert (landscape.allowed_edges,len(landscape.certified_switches),
                landscape.reversal_orbits)==expected
        assert c.conformations==states
        assert c.parent_robust_upper<F(1,2)<c.mutant_robust_lower
        row=dict(n=n,full_conformations=states,complete_sequence_universe=2**n,
                 same_parity_H_P_swap_edges=landscape.allowed_edges,
                 robust_switches=len(landscape.certified_switches),
                 chain_reversal_orbits=landscape.reversal_orbits,
                 focal_original=parent,focal_mutant=mutant,
                 original_robust_upper=str(c.parent_robust_upper),
                 mutant_robust_lower=str(c.mutant_robust_lower),
                 guaranteed_separation=str(c.separation),
                 maximum_certificate_error_factor_squared=str(threshold.max_symmetric_multiplier_squared),
                 all_certified_switches=[
                   dict(original=sequence(a,n),mutant=sequence(b,n),H_to_P=i,P_to_H=j)
                   for a,b,i,j in landscape.certified_switches])
        out['cases'].append(row)
    return out

if __name__=='__main__':
    print(json.dumps(run(),sort_keys=True,indent=2))
