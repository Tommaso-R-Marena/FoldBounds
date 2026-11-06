"""Independent algebraic oracle for error-correcting NMR geometry fingerprints.

These are SYNTHETIC geometry fixtures. The real-data CI separately
downloads and SHA-256 validates actual 2KDL/2KDM coordinates.
"""
from dataclasses import replace
from fractions import Fraction
from itertools import combinations
from math import isqrt

import pytest

from foldbounds.ga_gb_nmr import GA95, GB95, NMRModel
from foldbounds.ga_gb_error_codes import (
    Reporter, certify_code, classify_readouts, cover_certificate,
    extract_reporters, matching, verify_code, _d2
)

def models(seq, shifted):
    out=[]
    for k in range(20):
        ca={
            i: (round(i*1.1,3),
                round((25.0 + k/100) if (shifted and i>14) else
                      (k/100 if not shifted else 0),3),0.0)
            for i in range(1,57)
        }
        out.append(NMRModel(seq,ca,{}))
    return tuple(out)

@pytest.fixture(scope="module")
def synthetic():
    return models(GA95,False),models(GB95,True)

def test_conservative_integer_band_and_reporter_disjointness(synthetic):
    ga,gb=synthetic
    edges=extract_reporters(ga,gb)
    assert len(edges)>100
    assert all(e.a<=14<e.b and e.a not in (20,30,45)
               and e.b not in (20,30,45) and e.gap_milliang>0 for e in edges)
    for e in edges:
        da=[_d2(m.ca[e.a],m.ca[e.b]) for m in ga]
        db=[_d2(m.ca[e.a],m.ca[e.b]) for m in gb]
        assert e.ga_range_milliang[0]**2<=min(da)
        assert e.ga_range_milliang[1]**2>=max(da)
        assert e.gb_range_milliang[0]**2<=min(db)
        assert e.gb_range_milliang[1]**2>=max(db)
        if e.shorter=="GA95":
            assert e.gb_range_milliang[0]-e.ga_range_milliang[1]==e.gap_milliang
        else:
            assert e.ga_range_milliang[0]-e.gb_range_milliang[1]==e.gap_milliang

@pytest.mark.parametrize("channels",[1,3,5,7,9])
def test_minimax_bottleneck_witness_and_independent_cover(synthetic,channels):
    ga,gb=synthetic
    cert=certify_code(ga,gb,channels=channels)
    assert verify_code(cert,ga,gb)
    assert len(cert.reporters)==channels
    assert len({p for e in cert.reporters for p in e.endpoints})==2*channels
    assert len(cert.stronger_graph_cover_left)+len(cert.stronger_graph_cover_right)<channels
    assert cert.bottleneck_milliang==min(e.gap_milliang for e in cert.reporters)

def test_5_channel_majority_corrects_all_two_bit_errors_for_all_40_models(synthetic):
    ga,gb=synthetic
    cert=certify_code(ga,gb,channels=5)
    for true,models_ in (("GA95",ga),("GB95",gb)):
        for m in models_:
            clean={e.endpoints:Fraction(isqrt(_d2(m.ca[e.a],m.ca[e.b])))
                   for e in cert.reporters}
            assert classify_readouts(cert,clean)==true
            for f in range(3):
                for corrupt in combinations(clean,f):
                    flipped=dict(clean)
                    for key in corrupt:
                        flipped[key]=Fraction(100000 if true=="GA95" else 0)
                    assert classify_readouts(cert,flipped)==true

def test_minimum_cover_certificate_independent_of_matching_algorithm():
    # Exhaustively enumerate matchings for a four-left/four-right toy graph.
    edges=(
        Reporter(1,10,7,"GA95",(10,20),(27,30)),
        Reporter(1,11,3,"GA95",(10,20),(23,30)),
        Reporter(2,10,6,"GA95",(10,20),(26,30)),
        Reporter(2,12,5,"GA95",(10,20),(25,30)),
        Reporter(3,12,4,"GA95",(10,20),(24,30)),
        Reporter(3,13,2,"GA95",(10,20),(22,30)),
        Reporter(4,13,1,"GA95",(10,20),(21,30)),
    )
    for threshold in range(1,9):
        graph=[e for e in edges if e.gap_milliang>=threshold]
        best=0
        for mask in range(1<<len(graph)):
            chosen=[e for i,e in enumerate(graph) if mask>>i&1]
            if len({e.a for e in chosen})==len(chosen) and len({e.b for e in chosen})==len(chosen):
                best=max(best,len(chosen))
        assert len(matching(edges,threshold))==best
        left,right=cover_certificate(edges,threshold-1)
        assert len(left)+len(right)==best
        assert all(e.a in left or e.b in right for e in graph)

def test_corrupt_or_malformed_certificate_is_rejected(synthetic):
    ga,gb=synthetic
    c=certify_code(ga,gb,channels=5)
    for fake in (
        replace(c,bottleneck_milliang=c.bottleneck_milliang+1),
        replace(c,corrupted_channels=3),
        replace(c,coordinate_hashes={"2KDL":"fake","2KDM":"fake"}),
        replace(c,stronger_graph_cover_left=(1,2,3,4,5)),
        replace(c,reporters=c.reporters[:4]),
        replace(c,all_candidate_count=0),
        replace(c,classifier_error_milliang_strict_bound=Fraction(1)),
        replace(c,reporters=(c.reporters[0],)*5),
    ):
        with pytest.raises(ValueError):
            verify_code(fake,ga,gb)

def test_fail_closed_on_bad_distance_models_and_inputs(synthetic):
    ga,gb=synthetic
    with pytest.raises(ValueError):extract_reporters(ga[:-1],gb)
    with pytest.raises(ValueError):extract_reporters(ga,gb,left_end=True)
    with pytest.raises(ValueError):extract_reporters(ga,gb,min_sep=2.0)
    with pytest.raises(ValueError):extract_reporters(ga,gb,distance_min_ang=0)
    with pytest.raises(ValueError):certify_code(ga,gb,channels=4)
    with pytest.raises(ValueError):certify_code(ga,gb,channels=99)
    c=certify_code(ga,gb)
    with pytest.raises(ValueError):classify_readouts(c,{})
    with pytest.raises(ValueError):
        classify_readouts(c,{e.endpoints:1.2 for e in c.reporters})

def test_milliangstrom_grid_rejects_unverified_precision():
    assert _d2((0.001,-1.0,0.0),(1.0,0.0,0.0))==999**2+1000**2
    with pytest.raises(ValueError):
        _d2((0.0001,0,0),(0,0,0))
