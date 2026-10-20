"""Recompute rigorous positive-root brackets of robust mutation-majority thresholds."""
import json
from fractions import Fraction as F
from foldbounds.mutation_root import certify_contact_root,verify_contact_root


def run():
 out=[]
 for n,seq,a,b in [(8,'HHPHPHPH',F(6),F(7)),(10,'HHPHPHPHPH',F(7),F(8))]:
  c=certify_contact_root(seq,F(3,4),F(5,6),F(6,5),a,b,bits=48)
  assert verify_contact_root(c) and c.positive_root_sign_changes==1
  out.append({'length':n,'sequence':seq,'bend_weight':'3/4',
              'statewise_factor_range':['5/6','6/5'],
              'robust_majority_polynomial_coefficients':[str(v) for v in c.polynomial_coefficients],
              'descartes_positive_root_sign_changes':c.positive_root_sign_changes,
              'exact_low_bracket':str(c.bracket_low),
              'exact_high_bracket':str(c.bracket_high),
              'approx_critical_q':(float(c.bracket_low)+float(c.bracket_high))/2,
              'isolation_rounds':c.bisection_rounds})
 return {'results':out,'qualification':'contact/bend toy Hamiltonians, not experimental folding transitions'}

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))