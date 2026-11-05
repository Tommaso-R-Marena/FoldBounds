"""Reproducible compressed thermodynamic response demonstration, 2D HP only."""
from fractions import Fraction
from time import perf_counter

from foldbounds.lattice import build_frontier
from foldbounds.thermal_response import thermal_response_bounds


def experiment(n):
    cert = build_frontier("H"*n, Fraction(1), 1, 2)
    start = perf_counter()
    baseline = thermal_response_bounds(cert)
    enhanced = thermal_response_bounds(cert, use_monotone_spectra=True)
    stronger = thermal_response_bounds(cert, use_monotone_spectra=True,
                                       use_saw_count_cap=True)
    elapsed = perf_counter()-start
    print("n=",n,"q=1 elapsed_s=",round(elapsed,4),
          "compressed_family_paths=",enhanced.monotone_paths)
    print("variance baseline=",baseline.variance_lower,baseline.variance_upper)
    print("variance enhanced=",enhanced.variance_lower,enhanced.variance_upper)
    print("event sensitivity lower=",enhanced.event_susceptibility_lower)
    print("SAW count-cap partition upper=",stronger.partition_upper)
    print("SAW count-cap variance lower=",stronger.variance_lower)
    print("partition cap improvement factor=",enhanced.partition_upper/stronger.partition_upper)
    print("exact primal-dual exchange rounds=",enhanced.variance_dual.iterations)
    assert enhanced.monotone_paths>0
    assert enhanced.variance_lower>0
    assert enhanced.variance_upper<=baseline.variance_upper
    assert stronger.partition_upper<=enhanced.partition_upper
    assert stronger.variance_lower>=enhanced.variance_lower


if __name__=="__main__":
    for length in (10,30,100):
        experiment(length)
