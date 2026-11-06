"""Exhaustive n=8 H/P sequence scan of Bend monotonicity certificate.

Enumerates all 256 binary sequences at 8 residues, each in both 2D and 3D.
A failure of the sufficient coefficient condition is NOT proof that the
probability is nonmonotone; we separately search for exact counterexamples.
"""
import json
from fractions import Fraction as F
from itertools import product
from foldbounds.ensemble_crossover import spectrum,probability,robust_bend_sign_certificate
from foldbounds.ensemble_crossover_3d import spectrum_3d

def run(dim):
    total,nontrivial,certificate,negative,nonmonotone,examples=0,0,0,0,0,[]
    for word in product('HP',repeat=8):
        seq=''.join(word)
        h=spectrum(seq) if dim==2 else spectrum_3d(seq)
        total+=1
        # Threshold one actual H-H contact; nontrivial means both event and
        # complement classes occur among the microstates.
        isevent=any(k>=1 and count>0 for (k,r),count in h.items())
        iscomp=any(k<1 and count>0 for (k,r),count in h.items())
        if not (isevent and iscomp):continue
        nontrivial+=1
        terms,fail,positive=robust_bend_sign_certificate(h,1,F(3,2),F(1,4))
        if not fail:certificate+=1
        else:
            negative+=1
            if len(examples)<3:examples.append({'sequence':seq,'failed_coefficients':len(fail)})
        # Deterministic exact-rational *witness* of nonmonotonicity if any;
        # no claim of monotonicity from a finite positive sample.
        for q in (F(3,2),F(2),F(3)):
            pp=[probability(h,q,b,1) for b in (F(1,4),F(1,2),F(1),F(2),F(4),F(8))]
            if any(pp[i]>pp[i+1] for i in range(len(pp)-1)):
                nonmonotone+=1
                if len(examples)<8:examples.append({'counterexample_sequence':seq,'contact_q':str(q),'p_at_b_1_4':str(pp[0]),'p_at_b_1_2':str(pp[1]),'p_at_b_1':str(pp[2]),'p_at_b_2':str(pp[3]),'p_at_b_4':str(pp[4]),'p_at_b_8':str(pp[5])})
                break
    return {'dimension':dim,'sequences_total':total,'nontrivial_event_sequences':nontrivial,
            'sufficient_quadrant_coefficient_certificates':certificate,
            'sufficient_certificate_failed':negative,
            'demonstrably_nonmonotone_sequences_sampled_q_b':nonmonotone,
            'examples':examples}

if __name__=='__main__':print(json.dumps([run(2),run(3)],indent=2,sort_keys=True))