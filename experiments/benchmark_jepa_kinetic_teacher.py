"""Prospective, sequence-heldout test of JEPA as a finite pivot-policy student.

All trajectories are exact finite HP Markov computations, NOT MD. The heldout
model receives no class flux labels for ranking. The guard uses an exact oracle
to reject inferior proposals, so its cost is included in scientific claims.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
import time

from foldbounds.jepa_kinetic_teacher import (
    collect_pivot_examples, fit_pivot_jepa, gate_model_proposal,
    rank_classes_from_model,
)
from foldbounds.jepa_reversible_kinetics import (
    build_kernel, optimize_exact_flux,
)

TRAIN = ("HHHHHH", "HHPHHHHH", "HHHHHHHH", "HPHPHHHH")
HOLDOUT = ("HPHHPPHH", "HHPPHHHH", "HPHPHPHH")


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--seed",type=int,default=11)
    parser.add_argument("--epochs",type=int,default=12)
    parser.add_argument("--train-states",type=int,default=24)
    args=parser.parse_args()
    if args.epochs<1 or args.train_states<1:
        parser.error("epochs and train-states must be positive")
    assert not set(TRAIN)&set(HOLDOUT)
    start=time.perf_counter()
    examples=[]
    train_state_count=0
    for sequence in TRAIN:
        kernel=build_kernel(sequence,Fraction(2))
        train_state_count+=len(kernel.states)
        examples.extend(collect_pivot_examples(kernel,1,
                                                max_states=args.train_states,
                                                seed=args.seed))
    model,loss=fit_pivot_jepa(tuple(examples),epochs=args.epochs,seed=args.seed)
    train_time=time.perf_counter()-start
    print(json.dumps({
        "stage":"training", "sequence_train_count":len(TRAIN),
        "train_exact_state_universe":train_state_count,
        "labeled_state_action_examples":len(examples),
        "loss_final":loss[-1],"train_seconds":round(train_time,3),
        "oracle_labels_from_complete_finite_universes":True,
        "real_protein_data":False,
    },sort_keys=True))
    for sequence in HOLDOUT:
        start=time.perf_counter()
        world=build_kernel(sequence,Fraction(2))
        scores=rank_classes_from_model(model,world,1)
        gated=gate_model_proposal(world,1,scores)
        oracle=optimize_exact_flux(world,1)
        assert gated.achieved_flux>=gated.baseline_flux
        assert gated.achieved_flux<=oracle.optimized_flux
        print(json.dumps({
            "stage":"evaluation",
            "heldout_sequence":sequence,
            "states":len(world.states),
            "chosen_move_class":gated.chosen_class,
            "accepted_by_exact_gate":gated.accepted,
            "baseline_flux":str(gated.baseline_flux),
            "student_candidate_flux":str(gated.raw_candidate_flux),
            "achieved_flux":str(gated.achieved_flux),
            "exact_oracle_optimum_flux":str(gated.oracle_flux),
            "achieved_gain":str(gated.achieved_flux/gated.baseline_flux)
            if gated.baseline_flux else None,
            "oracle_gain":str(gated.oracle_flux/gated.baseline_flux)
            if gated.baseline_flux else None,
            "gap_to_oracle":str(gated.oracle_flux-gated.achieved_flux),
            "full_enumeration_used_for_gate":True,
            "elapsed_seconds":round(time.perf_counter()-start,3),
        },sort_keys=True))


if __name__=="__main__":
    main()
