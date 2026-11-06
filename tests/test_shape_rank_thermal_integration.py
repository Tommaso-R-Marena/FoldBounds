"""Shape/energy entropy upper caps combined with existing thermal lower witnesses."""
from fractions import Fraction
import pytest
from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.thermal_response import thermal_response_bounds

def test_shape_rank_refines_thermal_and_contains_exact_oracle():
    for sequence in ("HHHHHHHH", "HPHHPHHH", "HHHHHHHHHH"):
        threshold = 1
        q = Fraction(2)
        cert = build_frontier(sequence, q, threshold, 2)
        opts = dict(use_saw_count_cap=True,
                    use_geometric_contact_cap=True,
                    use_partially_directed_spectra=True)
        old = thermal_response_bounds(cert, **opts)
        new = thermal_response_bounds(cert, shape_memory=4, **opts)
        spectrum = exact_spectrum(sequence)
        z = sum((c*q**k for k,c in spectrum.items()), Fraction(0))
        moment = sum((k*c*q**k for k,c in spectrum.items()), Fraction(0))/z
        event = sum((c*q**k for k,c in spectrum.items() if k>=threshold),Fraction(0))
        assert new.partition_lower <= z <= new.partition_upper
        assert new.partition_upper <= old.partition_upper
        assert new.mean_lower <= moment <= new.mean_upper
        assert new.probability_lower <= event/z <= new.probability_upper
        assert new.probability_lower >= old.probability_lower
        assert new.probability_upper <= old.probability_upper
        assert new.shape_memory == 4

def test_shape_rank_and_horizon_jointly_conservative():
    for sequence in ("HHHHHHHH", "HPHPHHPH"):
        q = Fraction(3, 2)
        cert = build_frontier(sequence, q, 1, 2)
        opts = dict(use_saw_count_cap=True, use_geometric_contact_cap=True)
        prior = thermal_response_bounds(cert, horizon_steps=4, **opts)
        combined = thermal_response_bounds(cert, horizon_steps=4, shape_memory=4, **opts)
        z = sum((c*q**k for k,c in exact_spectrum(sequence).items()),Fraction(0))
        assert combined.partition_upper <= prior.partition_upper
        assert combined.partition_lower <= z <= combined.partition_upper

def test_shape_rank_fail_closed_bad_parameters():
    cert = build_frontier("HHHHHH", Fraction(2), 1, 2)
    for invalid in (True, 1, -1, 8, 2.5):
        with pytest.raises(ValueError):
            thermal_response_bounds(cert, shape_memory=invalid)
