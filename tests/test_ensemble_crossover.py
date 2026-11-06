from __future__ import annotations
from fractions import Fraction as F
from itertools import product
from random import Random

from foldbounds.ensemble_crossover import (
    geometry, spectrum, weighted_polys, sign_changes, eval_poly,
    monotone_in_bend_by_coefficients, unique_switch_bracket,
    robust_bend_sign_certificate, probability,
)

SEQ = 'HPHHPHHPHHHH'


def independent_enumeration(seq):
    # Separate recursive construction and contact scoring from production code.
    counts = {}
    walk = [(0, 0), (1, 0)]
    occupied = set(walk)
    directions = ((0,1),(1,0),(0,-1),(-1,0))
    def traverse(bend_count):
        if len(walk)==len(seq):
            hh_contacts = 0
            for j in range(2,len(walk)):
                for i in range(j-1):
                    hh_contacts += (seq[i]=='H' and seq[j]=='H'
                                    and abs(walk[i][0]-walk[j][0])+abs(walk[i][1]-walk[j][1])==1)
            key=(hh_contacts,bend_count)
            counts[key]=counts.get(key,0)+1
            return
        x,y=walk[-1]
        last=(walk[-1][0]-walk[-2][0],walk[-1][1]-walk[-2][1])
        for d in directions:
            nxt=(x+d[0],y+d[1])
            if nxt in occupied:continue
            occupied.add(nxt);walk.append(nxt)
            traverse(bend_count+(last != d))
            walk.pop();occupied.remove(nxt)
    traverse(0)
    return counts


def test_complete_120mer_geometry_and_histogram_oracle():
    h=spectrum(SEQ)
    assert sum(h.values())==30073
    assert dict(h)==independent_enumeration(SEQ)
    assert len(h)==52


def test_exact_probabilities():
    h=spectrum(SEQ)
    for q,b in product([F(3,2),F(2),F(3)], [F(1,4),F(1,2),F(1),F(2),F(4)]):
        e,c=weighted_polys(h,q,2)
        assert probability(h,q,b,2)==eval_poly(e,b)/(eval_poly(e,b)+eval_poly(c,b))
    assert probability(h,F(3,2),F(1),2)==F(23166,50005)


def test_strict_contact_monotonicity_random_ensembles():
    rng=Random(3092)
    for _ in range(40):
        seq=''.join(rng.choice('HP') for i in range(9))
        h=spectrum(seq)
        for b in (F(1,3),F(1),F(3)):
            for t in (1,2,3):
                x=probability(h,F(3,2),b,t)
                y=probability(h,F(2),b,t)
                assert x<=y


def test_exact_positive_derivative_certificate_over_unbounded_quadrant():
    h=spectrum(SEQ)
    count,negative,positive=robust_bend_sign_certificate(h,2,F(3,2),F(1,4))
    assert count==133 and negative==[] and positive==133
    for q in (F(3,2),F(2),F(3),F(7,2),F(5)):
        for b in (F(1,4),F(3,4),F(1),F(3),F(20)):
            ok,_=monotone_in_bend_by_coefficients(h,q,2)
            assert ok
            assert probability(h,q,b,2) <= probability(h,q,b+F(1,100),2)


def test_unique_switch_exact_rational_brackets():
    h=spectrum(SEQ)
    for q,lo_float,hi_float in [(F(3,2),1.173422,1.173423),
                                (F(2),0.647013,0.647014),
                                (F(3),0.263857,0.263858)]:
        e,c=weighted_polys(h,q,2)
        f=tuple(x-y for x,y in zip(e,c))
        assert sign_changes(f)==1
        lo,hi=unique_switch_bracket(h,q,2,F(1,4),F(4),bits=40)
        assert eval_poly(f,lo)<0 and eval_poly(f,hi)>=0
        assert lo_float<float(lo)<float(hi)<hi_float
        assert hi-lo<F(1,10**10)


def test_invalid_data_fail_closed():
    import pytest
    for s in ('A','H','XHH',''):
        with pytest.raises(ValueError):spectrum(s)
    h=spectrum(SEQ)
    with pytest.raises(ValueError):weighted_polys(h,F(0),2)
    with pytest.raises(ValueError):unique_switch_bracket(h,F(1),2,F(1,2),F(1))
    with pytest.raises(ValueError):probability(h,F(2),F(0),2)


def test_strict_threshold_edge_cases():
    for seq in ('HHHHHH','PPPPPP','HPHPHP','HHHHPH'):
        h=spectrum(seq)
        for q,b in product((F(1),F(2)),(F(1,2),F(3))):
            assert probability(h,q,b,-1)==1
            assert probability(h,q,b,100)==0

from dataclasses import replace
from foldbounds.ensemble_crossover import (
    certify_rectangular_population,verify_population_envelope,
    naive_decoupled_population_bounds,
)


def test_attained_sharp_envelope_and_independent_naive_comparison():
    box=certify_rectangular_population(SEQ,2,F(3,2),F(2),F(1,4),F(4))
    assert verify_population_envelope(box)
    assert box.lower_probability==F(84529143,508499807)
    assert box.upper_probability==F(1108499200,1368704059)
    assert box.positive_derivative_terms==133
    assert box.zero_derivative_terms==0
    naive_low,naive_high=naive_decoupled_population_bounds(SEQ,2,F(3,2),F(2),F(1,4),F(4))
    assert naive_low<box.lower_probability
    assert naive_high>box.upper_probability
    assert float(naive_high-naive_low)>0.99999
    assert float(box.upper_probability-box.lower_probability)<0.645


def test_box_enclosure_exact_thousands_of_grid_points():
    box=certify_rectangular_population(SEQ,2,F(3,2),F(2),F(1,4),F(4))
    h=spectrum(SEQ)
    for i in range(41):
        q=F(3,2)+F(i,80)
        for j in range(41):
            b=F(1,4)+F(15*j,160)
            p=probability(h,q,b,2)
            assert box.lower_probability<=p<=box.upper_probability


def test_robust_population_majority_classification():
    h=spectrum(SEQ)
    # For any q in [1.5,2], P increases with q and b in this region.
    # These rational comparisons certify population majorities across
    # continuous 2D rectangles, not merely sampled parameter values.
    assert probability(h,F(2),F(647,1000),2)<F(1,2)
    assert probability(h,F(3,2),F(1174,1000),2)>F(1,2)
    assert certify_rectangular_population(SEQ,2,F(3,2),F(2),F(1,4),F(647,1000)).upper_probability<F(1,2)
    assert certify_rectangular_population(SEQ,2,F(3,2),F(2),F(1174,1000),F(4)).lower_probability>F(1,2)


def test_fail_closed_without_nonnegative_coefficient_certificate():
    import pytest
    with pytest.raises(ValueError):
        certify_rectangular_population('HHHHPHHP',2,F(3,2),F(2),F(1,4),F(4))
    a=certify_rectangular_population(SEQ,2,F(3,2),F(2),F(1,4),F(4))
    assert not verify_population_envelope(replace(a,lower_probability=F(0)))
    assert not verify_population_envelope(replace(a,positive_derivative_terms=132))
    assert not verify_population_envelope(object())