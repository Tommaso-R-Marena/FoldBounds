"""Run fully checked maximum event-flux and reversible witness experiments.

All numbers are exact Fractions within the finite 2D HP *Monte Carlo* model.
This is a retrospective model sanity benchmark, not protein-dynamics evidence.
"""
from fractions import Fraction
import json
import time

from foldbounds.jepa_reversible_kinetics import build_kernel,event_flux
from foldbounds.jepa_flux_maxflow import (
    propose_maxflow_certificate,verify_maxflow_certificate,
    construct_reversible_kernel,event_flux_reversible,
)

CASES = (
    ("HHHHHH",1),
    ("HPHHPPHH",1),
    ("HHPPHHHH",1),
    ("HPHPHPHH",1),
    ("HHPHHHHH",2),
    ("HHHHHHHH",2),
    ("HPHHPPHHH",1),
)


def main():
    for sequence,threshold in CASES:
        started=time.perf_counter()
        base=build_kernel(sequence,Fraction(2))
        baseline=event_flux(base,threshold)[1]
        cert=propose_maxflow_certificate(base,threshold)
        verify_maxflow_certificate(base,cert)
        alpha=Fraction(1,2)
        witness=construct_reversible_kernel(base,cert,blend=alpha)
        blended=event_flux_reversible(witness,threshold)
        optimal=cert.total_mass/base.partition
        assert optimal>=blended>=baseline
        assert blended==(baseline+optimal)/2
        print(json.dumps({
            "sequence":sequence,
            "threshold":threshold,
            "q":"2",
            "states":len(base.states),
            "primal_cross_edges":len(cert.flows),
            "uniform_flux":str(baseline),
            "maxflow_upper_bound":str(optimal),
            "full_flow_witness_flux":str(optimal),
            "half_blend_flux":str(blended),
            "maxflow_gain":str(optimal/baseline) if baseline else None,
            "half_blend_gain":str(blended/baseline) if baseline else None,
            "flow_equals_cut":cert.total_mass==cert.cut_mass,
            "positive_baseline_transition_mass":str(1-alpha),
            "model":"finite canonical HP reversible MC; not molecular folding",
            "time_unit":"one algorithmic pivot attempt",
            "full_enumeration_used":True,
            "elapsed_seconds":round(time.perf_counter()-started,3),
        },sort_keys=True))


if __name__=="__main__":
    main()
