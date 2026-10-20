"""Exact no-direct-contact epistasis & real-data thermodynamic sanity checks."""
from __future__ import annotations
import json,sys,time,math
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from foldbounds.counterfactual_reweighting import (
    reweight_cycle,critical_energy_error_kcal,population_with_factor_error)
from foldbounds.mutation_switch import certify_mutation_switch
from foldbounds.mutation_switch_10 import certify_ten_mutation


def main():
    result={}
    t=time.monotonic();all_effect=0;no_effect=0;positive=0;negative=0;largest=(F(0),None)
    # At baseline q=8,b=1, all 768 same-sublattice two-site P->H cycles.
    for i,j in combinations(range(1,9),2):
        if (i-j)%2:continue
        for others in range(1<<6):
            k=0
            slots=[t for t in range(8) if t not in (i-1,j-1)]
            for bit,idx in enumerate(slots):
                if others&(1<<bit):k|=1<<idx
            seq=''.join('H' if k&(1<<s) else 'P' for s in range(8))
            x=reweight_cycle(seq,i,j,F(8),F(1))
            assert not x.direct_contact_possible
            if x.odds_epistasis_factor==1:no_effect+=1
            else:
                all_effect+=1
                if x.odds_epistasis_factor>1:positive+=1
                else:negative+=1
            # Prioritize deviation in relative odds in both directions
            deviation=max(x.odds_epistasis_factor,1/x.odds_epistasis_factor)
            if deviation>largest[0]:largest=deviation,(seq,i,j)
    assert all_effect+no_effect==768
    result['eight_mer_pairwise_nondirect']={
        'same_parity_two_site_cycles':768,'nonzero_exact_odds_epistasis':all_effect,
        'zero_exact_odds_epistasis':no_effect,'positive_sign_count':positive,
        'negative_sign_count':negative,'max_absolute_odds_interaction_factor':str(largest[0]),
        'max_case':list(largest[1]),'elapsed_seconds_not_reproducibility_key':round(time.monotonic()-t,3)}
    for n,seq in ((8,'PHPHPHPH'),(10,'PHPHPHPHPH')):
        x=reweight_cycle(seq,1,3,F(8),F(1))
        result[str(n)+'_mer']={
            'initial':seq,'mutated_sites':[1,3], 'direct_contact_possible':x.direct_contact_possible,
            'four_event_probabilities':[str(x.baseline_event_probability),str(x.single_i_event_probability),str(x.single_j_event_probability),str(x.double_event_probability)],
            'event_conditional_covariance':str(x.event_conditional_covariance),
            'complement_conditional_covariance':str(x.complement_conditional_covariance),
            'event_interaction':str(x.event_interaction_factor),
            'complement_interaction':str(x.complement_interaction_factor),
            'odds_epistasis_factor':str(x.odds_epistasis_factor),
            'microstates':x.enumerated_microstates}
    c8=certify_mutation_switch('PHPHPHPH','HHPHPHPH',F(8),F(12),F(3,4),F(5,4))
    c10=certify_ten_mutation()
    RT=.00198720425864083*293.15
    result['energy_calibration']={
        'temperature_K':293.15, 'RT_kcal_per_mol':RT,
        'energy_error_of_weight_multiplier_1p2_kcal_mol':RT*math.log(1.2),
        'critical_abs_independent_per_conformation_energy_error_kcal_mol_8mer':critical_energy_error_kcal(c8.parent.event_upper,c8.mutant.event_lower),
        'critical_abs_independent_per_conformation_energy_error_kcal_mol_10mer':critical_energy_error_kcal(c10.original_nominal_upper,c10.mutant_nominal_lower),
        'published_GA95_GB95_opposite_90_percent_population_minimum_mutation_free_energy_shift_kcal_mol':2*RT*math.log(9),
        'scope':'Thermodynamic implication of literature populations; no predicted GA95/GB95 relative free energy'}
    return result

if __name__=='__main__':
    d=main()
    d['eight_mer_pairwise_nondirect'].pop('elapsed_seconds_not_reproducibility_key')
    print(json.dumps(d,sort_keys=True,indent=2))