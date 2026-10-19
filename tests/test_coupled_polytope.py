"""Independent vertex oracle for sharp rational probability intervals."""
from fractions import Fraction
import random
import pytest
from foldbounds.coupled_polytope import tighten_coupled_masses


def independent_vertices(le,ue,lc,uc,lz,uz):
    # Lines: (a,b,c) representing a E + b C = c.
    lines=((1,0,le),(1,0,ue),(0,1,lc),(0,1,uc),
           (1,1,lz),(1,1,uz))
    out=set()
    for i,(a,b,c) in enumerate(lines):
        for d,e,f in lines[i+1:]:
            det=a*e-b*d
            if not det:continue
            E=Fraction(c*e-b*f,det)
            C=Fraction(a*f-c*d,det)
            if le<=E<=ue and lc<=C<=uc and lz<=E+C<=uz:
                out.add((E,C))
    return out


def test_exhaustive_exact_small_box_slab_polytope_corners():
    for le in range(4):
        for ue in range(le,6):
            for lc in range(4):
                for uc in range(lc,6):
                    for lz in range(1,10):
                        for uz in range(lz,11):
                            points=independent_vertices(le,ue,lc,uc,lz,uz)
                            values=(Fraction(le),Fraction(ue),Fraction(lc),
                                    Fraction(uc),Fraction(lz),Fraction(uz))
                            if not points:
                                with pytest.raises(ValueError):
                                    tighten_coupled_masses(*values)
                                continue
                            result=tighten_coupled_masses(*values)
                            r=[E/(E+C) for E,C in points]
                            assert result.probability_lower==min(r)
                            assert result.probability_upper==max(r)
                            for E,C in points:
                                assert result.event_lower<=E<=result.event_upper
                                assert result.complement_lower<=C<=result.complement_upper
                                assert result.partition_lower<=E+C<=result.partition_upper


def test_randomized_fractional_intervals_contain_and_are_optimal():
    rng=random.Random(31177)
    for _ in range(3000):
        E=Fraction(rng.randint(0,30),rng.randint(1,9))
        C=Fraction(rng.randint(0,30),rng.randint(1,9))
        if E+C==0:continue
        le=max(Fraction(0),E-Fraction(rng.randint(0,6),3))
        ue=E+Fraction(rng.randint(0,6),3)
        lc=max(Fraction(0),C-Fraction(rng.randint(0,6),3))
        uc=C+Fraction(rng.randint(0,6),3)
        lz=max(Fraction(1,1000000),E+C-Fraction(rng.randint(0,6),3))
        uz=E+C+Fraction(rng.randint(0,6),3)
        r=tighten_coupled_masses(le,ue,lc,uc,lz,uz)
        vertices=independent_vertices(le,ue,lc,uc,lz,uz)
        assert vertices
        assert r.probability_lower==min(x/(x+y) for x,y in vertices)
        assert r.probability_upper==max(x/(x+y) for x,y in vertices)
        assert r.probability_lower<=E/(E+C)<=r.probability_upper


def test_strictly_tighter_example_than_independent_coupled_ratio():
    # E in [6,100], C in [8,100], Z in [16,25]. Naive P lower=6/106,
    # but total Z <=25 forces C<=19 at E=6, so optimum is 6/25.
    r=tighten_coupled_masses(*map(Fraction,(6,100,8,100,16,25)))
    assert r.probability_lower==Fraction(6,25)
    assert r.probability_upper==Fraction(17,25)
    assert r.event_upper==Fraction(17)
    assert r.complement_upper==Fraction(19)


@pytest.mark.parametrize('input',(True,0,1.5,'1',None))
def test_invalid_types_rejected(input):
    with pytest.raises(ValueError):
        tighten_coupled_masses(input,Fraction(2),Fraction(1),Fraction(3),
                               Fraction(1),Fraction(5))


def test_infeasible_intervals_rejected():
    with pytest.raises(ValueError):
        tighten_coupled_masses(*map(Fraction,(1,2,1,2,10,11)))
    with pytest.raises(ValueError):
        tighten_coupled_masses(*map(Fraction,(5,5,10,10,1,10)))