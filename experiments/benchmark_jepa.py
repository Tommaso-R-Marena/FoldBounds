"""Toy-model, whole-sequence holdout JEPA-vs-search benchmark.

Run: python experiments/benchmark_jepa.py [--train] [--rl].
Preparation cost for training data is reported separately, never hidden.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import time

from foldbounds.adaptive import adaptive_refine, verify_refinement_run
from foldbounds.jepa_guidance import (
    collect_checked_decisions, fit_checked_jepa, guided_refine,
    train_online_jepa_rl,
)
from foldbounds.lattice import build_frontier

TRAIN = ("HHHHHH", "HHPHHHHH", "HPHPHPHH")
TEST = ("HHPPHHHH", "HPHHPPHH")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", action="store_true")
    parser.add_argument("--rl", action="store_true")
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--epochs", type=int, default=25)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    if args.steps < 0:
        parser.error("--steps must be nonnegative")
    model = None
    training_calls = 0
    if args.train:
        start = time.perf_counter()
        data = collect_checked_decisions(TRAIN, seed=args.seed, steps=4)
        training_calls += len(TRAIN) + sum(len(d.candidates) for d in data)
        model, losses = fit_checked_jepa(data, epochs=args.epochs, seed=args.seed)
        print(json.dumps({
            "stage": "offline", "decisions": len(data),
            "checking_calls_for_labels": training_calls,
            "last_loss": losses[-1],
            "elapsed_s": round(time.perf_counter() - start, 3),
        }))
    if args.rl:
        start = time.perf_counter()
        model, rewards, calls = train_online_jepa_rl(
            TRAIN, episodes=8, steps=4, seed=args.seed, model=model,
        )
        training_calls += calls
        print(json.dumps({
            "stage": "online-rl", "rewards": rewards,
            "checking_calls": calls,
            "elapsed_s": round(time.perf_counter() - start, 3),
        }))
    for sequence in TEST:
        cert = build_frontier(sequence, Fraction(2), 1, 2)
        for policy, scorer in (
            ("breadth", None),
            ("jepa", model.score if model is not None else None),
        ):
            if policy == "jepa" and model is None:
                continue
            start = time.perf_counter()
            run = guided_refine(cert, steps=args.steps, scorer=scorer)
            print(json.dumps({
                "sequence": sequence, "policy": policy,
                "width": str(run.refinement.history[-1].probability_width),
                "steps": len(run.refinement.history) - 1,
                "checker_calls_including_replay": run.checker_calls,
                "training_calls_excluded": training_calls,
                "elapsed_s": round(time.perf_counter() - start, 3),
            }))
        start = time.perf_counter()
        greedy = adaptive_refine(cert, steps=args.steps, policy="probability")
        verify_refinement_run(greedy)
        calls = 1 + sum(
            sum(node.kind == "bound" for node in frontier.frontier)
            for frontier in greedy.audit_frontiers[:-1]
        ) + len(greedy.audit_frontiers)
        print(json.dumps({
            "sequence": sequence, "policy": "probability-oracle",
            "width": str(greedy.history[-1].probability_width),
            "steps": len(greedy.history) - 1,
            "checker_calls_including_replay": calls,
            "training_calls_excluded": 0,
            "elapsed_s": round(time.perf_counter() - start, 3),
        }))


if __name__ == "__main__":
    main()
