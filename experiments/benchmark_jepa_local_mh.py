"""Check local JEPA/Hastings sampling versus exhaustive small-world Boltzmann oracle.

Large-chain sampling deliberately does NOT enumerate the conformation space.
It reports no physical times, folding rates, or validated mixing advantage.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import time

from foldbounds.jepa_local_hastings import (
    StabilityGuard,audit_finite_mh,neural_policy,sample_pivot_chain,
)
from foldbounds.jepa_reversible_kinetics import (
    build_kernel,move_classes,uniform_classes,
)


def structural_policy(path):
    classes=move_classes(len(path))
    weights=[3+abs(path[k][0])+len(path)*int(family=="half")
             for k,family in classes]
    z=sum(weights)
    return tuple(Fraction(value,z) for value in weights)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--steps",type=int,default=8)
    parser.add_argument("--train",action="store_true")
    args=parser.parse_args()
    if args.steps<0:
        parser.error("steps must be nonnegative")
    for sequence in ("HHHHHH","HPHHPPHH"):
        start=time.perf_counter()
        report=audit_finite_mh(sequence,Fraction(2),
                               StabilityGuard(structural_policy))
        original=build_kernel(sequence,Fraction(2))
        assert report.partition==original.partition
        print(json.dumps({
            "stage":"exhaustive-small-model",
            "sequence":sequence,
            "states":report.states,
            "checked_directed_edges":report.checked_edges,
            "partition":str(report.partition),
            "stationary_event_mass":str(report.event_mass),
            "proposal":"positive rational state-dependent geometric heuristic",
            "reversibility_checked":True,
            "elapsed_seconds":round(time.perf_counter()-start,3),
        },sort_keys=True))
    sequence="HHPPHHPP"*4
    start=time.perf_counter()
    trajectory=sample_pivot_chain(
        sequence,Fraction(2),StabilityGuard(structural_policy),
        steps=args.steps,seed=7,
    )
    print(json.dumps({
        "stage":"no-global-enumeration",
        "sequence_length":len(sequence),
        "steps":args.steps,
        "distinct_states_seen":len(set(trajectory)),
        "exact_rational_hastings":True,
        "global_universe_enumerated":False,
        "physical_time_scale":"none",
        "elapsed_seconds":round(time.perf_counter()-start,3),
    },sort_keys=True))
    if args.train:
        from foldbounds.jepa_kinetic_teacher import (
            collect_pivot_examples,fit_pivot_jepa,
        )
        start=time.perf_counter()
        training=build_kernel("HHHHHH",Fraction(2))
        examples=collect_pivot_examples(training,1,max_states=12,seed=3)
        model,loss=fit_pivot_jepa(examples,epochs=4,seed=3)
        heldout="HPHHPPHH"
        policy=StabilityGuard(neural_policy(
            model,heldout,Fraction(2),1,
        ))
        audited=audit_finite_mh(
            heldout,Fraction(2),policy,max_states=1000,
        )
        print(json.dumps({
            "stage":"heldout-learned-policy-audit",
            "training_sequence":"HHHHHH",
            "heldout_sequence":heldout,
            "loss_final":loss[-1],
            "audited_states":audited.states,
            "reversibility_checked":True,
            "frozen_model_required":True,
            "claim_of_better_sampling":False,
            "elapsed_seconds":round(time.perf_counter()-start,3),
        },sort_keys=True))


if __name__=="__main__":
    main()
