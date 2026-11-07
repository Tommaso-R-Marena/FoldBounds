"""Exact, reproducible HP pivot-Monte-Carlo equilibrium-flux benchmark.

This tests event crossing *per MC proposal*, NOT physical kinetic rates.
All input sequences, temperatures, and thresholds are fixed in source.
Training/test split is prospective for any future JEPA student model.
"""
from __future__ import annotations

from fractions import Fraction
import json
import time

from foldbounds.jepa_reversible_kinetics import (
    build_kernel, event_flux, optimize_exact_flux, verify_kernel,
)

TEST_CASES = (
    ("HHHHHH", 1),
    ("HHHHHH", 2),
    ("HPHHPPHH", 1),
    ("HPHHPPHH", 2),
    ("HHPPHHHH", 1),
    ("HHPPHHHH", 2),
    ("HHHHHHHH", 1),
    ("HHHHHHHH", 2),
    ("HHPHHHHH", 1),
    ("HHPHHHHH", 2),
    ("HPHPHPHH", 1),
    ("HPHPHPHH", 2),
)
# These are fixed before evaluating any future generalizing JEPA policies.
JEPA_HOLDOUT_SEQUENCES = ("HPHHPPHH", "HHPPHHHH", "HPHPHPHH")


def main() -> None:
    for sequence, threshold in TEST_CASES:
        start = time.perf_counter()
        kernel = build_kernel(sequence,Fraction(2))
        verify_kernel(kernel)
        baseline_event,baseline_flux,baseline_exit = event_flux(kernel,threshold)
        result = optimize_exact_flux(kernel,threshold)
        optimized_kernel = build_kernel(
            sequence,Fraction(2),class_probs=result.optimized_weights,
        )
        verify_kernel(optimized_kernel)
        optimized_event,optimized_flux,optimized_exit = event_flux(
            optimized_kernel,threshold,
        )
        assert baseline_event == optimized_event
        assert baseline_flux == result.baseline_flux
        assert optimized_flux == result.optimized_flux
        assert optimized_flux >= baseline_flux
        print(json.dumps({
            "sequence":sequence,
            "threshold":threshold,
            "state_count":len(kernel.states),
            "q":"2",
            "event_stationary_probability":str(baseline_event),
            "uniform_event_flux":str(baseline_flux),
            "optimal_event_flux":str(optimized_flux),
            "uniform_exit_chance":str(baseline_exit),
            "optimal_exit_chance":str(optimized_exit),
            "flux_ratio":str(result.relative_gain) if result.relative_gain else None,
            "best_move_class":result.best_class,
            "positive_class_floor":str(min(result.optimized_weights)),
            "learned_oracle":"EXACT EXHAUSTIVE OPTIMIZER; NOT JEPA",
            "kinetic_time_units":"one algorithmic pivot proposal",
            "wall_time_seconds":round(time.perf_counter()-start,3),
        },sort_keys=True))


if __name__=="__main__":
    main()
