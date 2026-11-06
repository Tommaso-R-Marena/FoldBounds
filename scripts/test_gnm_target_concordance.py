"""Property tests of efficient exact interval concordance vs exhaustive comparator."""
import random
from fractions import Fraction
import pytest
from gnm_target_concordance import counts_fast, counts_slow


def test_fast_exactly_matches_brute_force():
    rng=random.Random(20261105)
    for n in range(2,15):
        for _ in range(20):
            lo=[rng.randint(-5,10) for _ in range(n)]
            hi=[v+rng.randint(0,8) for v in lo]
            target=[rng.randint(-2,6) for _ in range(n)]
            if len(set(target))<=1:continue
            fast=counts_fast(lo,hi,target)
            slow=counts_slow(lo,hi,target)
            assert (fast['guaranteed_concordant'],fast['guaranteed_discordant'],
                    fast['uncertain'],fast['target_ties']) == slow
            m=slow[0]+slow[1]+slow[2]
            assert Fraction(fast['lower']) == Fraction(slow[0]-slow[1]-slow[2],m)
            assert Fraction(fast['upper']) == Fraction(slow[0]-slow[1]+slow[2],m)
            assert Fraction(fast['lower'])<=Fraction(fast['upper'])


def test_fail_closed_ties_inversion_type_nan():
    assert counts_fast([1,2,3],[1,2,3],[2,2,3])['target_ties']==1
    with pytest.raises(ValueError):counts_fast([1,3],[2,2],[3,4])
    with pytest.raises(ValueError):counts_fast([1,3],[2,4],[3,float('nan')])
    with pytest.raises(ValueError):counts_fast([1,3],[2,4],[3,float('inf')])
    with pytest.raises(ValueError):counts_fast([1.,3.],[2,4],[3,4])
    with pytest.raises(ValueError):counts_fast([1,3],[2,4],[0,0])


def test_all_points_and_unordered_intervals():
    s=counts_fast([1,2,3,4],[1,2,3,4],[1,2,3,4])
    assert s['guaranteed_concordant']==6
    assert s['guaranteed_discordant']==0 and s['uncertain']==0
    assert s['lower']=='1' and s['upper']=='1'
    s=counts_fast([1,1,1],[2,2,2],[3,2,1])
    assert s['uncertain']==3 and s['lower']=='-1' and s['upper']=='1'