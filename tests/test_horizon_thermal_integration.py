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


def test_horizon_event_specific_probability_caps():
    for sequence, threshold in (('HHHHHHHH',1), ('HPHHPHHH',1), ('HHHHHHHHHH',2)):
        cert=build_frontier(sequence,Fraction(2),threshold,2)
        no_horizon=thermal_response_bounds(cert,use_saw_count_cap=True,use_geometric_contact_cap=True)
        checked=thermal_response_bounds(cert,use_saw_count_cap=True,use_geometric_contact_cap=True,horizon_steps=6)
        spectrum=exact_spectrum(sequence)
        z=sum((count*Fraction(2)**score for score,count in spectrum.items()),Fraction())
        e=sum((count*Fraction(2)**score for score,count in spectrum.items() if score>=threshold),Fraction())
        p=e/z
        assert checked.probability_lower<=p<=checked.probability_upper
        assert checked.probability_lower>=no_horizon.probability_lower
        assert checked.probability_upper<=no_horizon.probability_upper


def test_horizon_with_directed_spectral_subset_and_overlap_witness():
    sequence = 'HHHHHHHH'
    cert = build_frontier(sequence, Fraction(2), 1, 2)
    common = dict(use_saw_count_cap=True,
                  use_geometric_contact_cap=True,
                  use_partially_directed_spectra=True,
                  witnesses=(tuple((i, 0) for i in range(8)),))
    baseline = thermal_response_bounds(cert, **common)
    refined = thermal_response_bounds(cert, horizon_steps=4, **common)
    exact = exact_spectrum(sequence)
    z = sum((count * Fraction(2)**k for k, count in exact.items()), Fraction(0))
    e = sum((count * Fraction(2)**k for k, count in exact.items() if k >= 1), Fraction(0))
    p = e / z
    assert refined.partition_lower <= z <= refined.partition_upper
    assert refined.partition_upper <= baseline.partition_upper
    assert refined.probability_lower <= p <= refined.probability_upper
    assert refined.probability_lower >= baseline.probability_lower
    assert refined.probability_upper <= baseline.probability_upper
