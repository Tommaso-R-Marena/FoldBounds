"""Runs with the full FoldBounds repository: keeps solver unchanged."""
from fractions import Fraction
from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.thermal_response import thermal_response_bounds


def test_opt_in_horizon_tightens_and_stays_sound():
    for n in (8,10,12):
        seq='H'*n
        cert=build_frontier(seq,Fraction(2),1,2)
        baseline=thermal_response_bounds(cert,use_saw_count_cap=True,
                                         use_geometric_contact_cap=True)
        improved=thermal_response_bounds(cert,use_saw_count_cap=True,
                                         use_geometric_contact_cap=True,
                                         horizon_steps=6)
        spectrum=exact_spectrum(seq)
        z=sum((count*Fraction(2)**score for score,count in spectrum.items()),Fraction())
        m=sum((count*Fraction(2)**score*score for score,count in spectrum.items()),Fraction())/z
        assert improved.partition_upper<=baseline.partition_upper
        assert improved.partition_lower<=z<=improved.partition_upper
        assert improved.mean_lower<=m<=improved.mean_upper
        assert improved.variance_lower<=improved.variance_upper
        assert improved.horizon_steps==6


def test_zero_horizon_backwards_compatibility():
    cert=build_frontier('HPHHPHPH',Fraction(3,2),1,2)
    assert thermal_response_bounds(cert)==thermal_response_bounds(cert,horizon_steps=0)
