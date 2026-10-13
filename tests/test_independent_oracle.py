"""Separately implemented lattice walk and contact scorer (no production helpers)."""
from fractions import Fraction
from itertools import product
import pytest

from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead


def independent_spectrum(sequence):
    starts = [((0,0),(1,0))]
    for _ in range(2,len(sequence)):
        next_layer = []
        for walk in starts:
            x,y = walk[-1]
            for x1,y1 in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)):
                if (x1,y1) not in walk:
                    next_layer.append(walk+((x1,y1),))
        starts = next_layer
    histogram = {}
    for walk in starts:
        score = 0
        for i in range(len(sequence)):
            for j in range(i+2,len(sequence)):
                if sequence[i]!='H' or sequence[j]!='H':
                    continue
                x0,y0=walk[i]
                x1,y1=walk[j]
                if abs(x0-x1)+abs(y0-y1)==1:
                    score += 1
        histogram[score] = histogram.get(score,0)+1
    return histogram


@pytest.mark.parametrize('n',range(2,8))
def test_all_short_spectra_match_separate_implementation(n):
    for residues in product('HP',repeat=n):
        sequence=''.join(residues)
        independent=independent_spectrum(sequence)
        assert exact_spectrum(sequence)==independent
        for q,event in ((Fraction(2),1),(Fraction(3,2),2)):
            z=sum((num*q**k for k,num in independent.items()),Fraction(0))
            e=sum((num*q**k for k,num in independent.items() if k>=event),Fraction(0))
            b=verify_lookahead(build_frontier(sequence,q,event,2),lookahead=2)
            assert b.partition_lower <= z <= b.partition_upper
            assert b.probability_lower <= e/z <= b.probability_upper
