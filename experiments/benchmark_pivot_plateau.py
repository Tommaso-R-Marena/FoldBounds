"""Deterministic all-rational nonenumerative kinetic/equilibrium certificate benchmark."""
from fractions import Fraction
from time import perf_counter

from foldbounds.pivot_plateau_producer import propose_rectangular_plateau
from foldbounds.pivot_plateau_verifier import verify_pivot_plateau


def main():
    q, steps = Fraction(10**9), 10**6
    print('n,plateau_states,witness_gap,boundary_slots,exit_hit_upper,'
          'equilibrium_bad_upper,full_ensemble_enumerated,elapsed_seconds')
    for h in (10, 20, 50, 100, 200, 300):
        t0 = perf_counter()
        cert = propose_rectangular_plateau(h)
        proof = verify_pivot_plateau(cert,q,steps)
        assert proof.total_enumerated_conformations == 0
        bad_upper = proof.equilibrium_bad_event_upper
        # Precise decimal order, without conversion to floating point.
        exponent = max(0, len(str(bad_upper.denominator))-len(str(bad_upper.numerator))-1)
        while bad_upper < Fraction(1,10**(exponent+1)):
            exponent+=1
        while exponent and bad_upper >= Fraction(1,10**exponent):
            exponent-=1
        tag = f'<1e-{exponent}' if exponent else '<1'
        print(f'{proof.n},{proof.plateau_states},{proof.witness_gap},'
              f'{proof.boundary_proposal_slots},{float(proof.fixed_seed_hit_by_steps_upper):.9g},'
              f'{tag},0,{perf_counter()-t0:.3f}',flush=True)
    print('This benchmark is a 2D lattice HP toy-model certificate, NOT physical folding kinetics.')


if __name__ == '__main__':
    main()
