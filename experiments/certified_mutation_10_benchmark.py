"""Independent 308,981-conformation 3D length-10 phase certificate reproduction."""
import json
from fractions import Fraction as F
from foldbounds.mutation_switch_10 import (
    geometry_fingerprints10, certify_ten_mutation,verify_ten_mutation,
    population10, ORIGINAL, MUTANT)

def run():
    cert=certify_ten_mutation()
    assert verify_ten_mutation(cert)
    assert cert.original_robust_upper < F(1,2) < cert.mutant_robust_lower
    assert cert.fingerprints==1852 and cert.conformations==308981
    return {
        'model':'3D cubic HP length-10, >=3 nonconsecutive geometric contacts',
        'original':ORIGINAL, 'mutant':MUTANT, 'mutation':'P1H',
        'q_range':['8','12'], 'b_range':['3/4','5/4'],
        'independent_per_state_multiplier_range':['5/6','6/5'],
        'complete_SAW_conformations':cert.conformations,
        'lossless_contact_graph_fingerprints':cert.fingerprints,
        'original_nominal_range':[str(cert.original_nominal_lower),str(cert.original_nominal_upper)],
        'mutant_nominal_range':[str(cert.mutant_nominal_lower),str(cert.mutant_nominal_upper)],
        'original_robust_upper':str(cert.original_robust_upper),
        'mutant_robust_lower':str(cert.mutant_robust_lower),
        'minimum_population_separation':str(cert.mutant_robust_lower-cert.original_robust_upper),
        'at_q8_b1':[str(population10(ORIGINAL,F(8),F(1))),
                   str(population10(MUTANT,F(8),F(1)))],
        'qualification':'exact finite-model result, not a real protein folding probability'
    }

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))