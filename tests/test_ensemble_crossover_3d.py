from fractions import Fraction as F
from foldbounds.ensemble_crossover import (
    probability, unique_switch_bracket,robust_bend_sign_certificate)
from foldbounds.ensemble_crossover_3d import (
    geometry_3d,spectrum_3d,certified_3d_example)


def independent_3d(seq):
    # Independent complete n=8 conformation iterator with a different state
    # representation and direction order than the implementation under test.
    d=((0,0,1),(0,1,0),(1,0,0),(0,0,-1),(0,-1,0),(-1,0,0))
    points=((0,0,0),(1,0,0))
    counts={}
    stack=[(points,0)]
    while stack:
        path,bend=stack.pop()
        if len(path)==len(seq):
            k=sum(seq[i]=='H' and seq[j]=='H' and
                  sum(abs(path[i][z]-path[j][z]) for z in range(3))==1
                  for j in range(2,len(seq)) for i in range(j-1))
            counts[k,bend]=counts.get((k,bend),0)+1
            continue
        last=tuple(path[-1][i]-path[-2][i] for i in range(3))
        for step in d:
            candidate=tuple(path[-1][i]+step[i] for i in range(3))
            if candidate not in path:
                stack.append((path+(candidate,),bend+(step!=last)))
    return counts


def test_independently_enumerated_three_dimensional_ensemble():
    seq='HPHHPHHH'
    h=spectrum_3d(seq)
    assert sum(h.values())==13565
    assert len(h)==22
    assert dict(h)==independent_3d(seq)


def test_3d_exact_continuum_certificate_and_unique_crossovers():
    r=certified_3d_example()
    assert r['conformations']==13565
    assert r['coefficients_total']==55
    assert r['coefficients_positive']==55
    h=spectrum_3d('HPHHPHHH')
    for q in (F(3,2),F(2)):
        lo,hi=r['crossover_brackets'][str(q)]
        assert hi-lo<F(1,10**10)
        assert probability(h,q,lo,1)<F(1,2)
        assert probability(h,q,hi,1)>=F(1,2)
    assert r['interval'][0]<F(1,2)<r['interval'][1]


def test_fail_closed_on_unsupported_3d_sizes_and_alphabet():
    import pytest
    for s in ('H','XYZ','HHHXXHHH',''):
        with pytest.raises(ValueError):spectrum_3d(s)
    with pytest.raises(ValueError):geometry_3d(11)