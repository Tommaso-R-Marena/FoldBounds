"""Exact short-walk oracle and adversarial SAW block-certificate checks."""
from fractions import Fraction as Q
from random import Random
import pytest
from foldbounds.lattice import build_frontier
from foldbounds.saw_count_cap import (
    SAWCountCertificate, exact_rooted_block_counts,
    optimal_saw_count_cap, verify_saw_count_cap,
)
from foldbounds.thermal_response import thermal_response_bounds


def independent_walks(n):
    paths = [((0,0),(1,0))]
    for _ in range(2,n):
        following = []
        for path in paths:
            x,y = path[-1]
            for site in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)):
                if site not in path:
                    following.append(path+(site,))
        paths = following
    return paths


def test_exact_block_counts_vs_independent_oracle():
    counts = exact_rooted_block_counts(9)
    for k in range(10):
        assert counts[k] == len(independent_walks(k+2))


def test_products_and_malformed_adversarial_certificates():
    for n in range(101):
        cert=optimal_saw_count_cap(n)
        assert verify_saw_count_cap(cert)
        assert cert.count_upper<=3**n
    for forged in (
        SAWCountCertificate(5,(3,),25,8),
        SAWCountCertificate(4,(4,),27,8),
        SAWCountCertificate(5,(5,),0,8),
    ):
        with pytest.raises(ValueError):
            verify_saw_count_cap(forged)


def test_random_frontiers_still_enclose_independently_enumerated_Z():
    rng=Random(412)
    for _ in range(70):
        n=rng.randint(2,9)
        sequence=''.join(rng.choice('HP') for _ in range(n))
        q=rng.choice((Q(1),Q(3,2),Q(2)))
        z=Q(0)
        for path in independent_walks(n):
            score=sum(sequence[i]==sequence[j]=='H'
                      and abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
                      for i in range(n) for j in range(i+2,n))
            z+=q**score
        cert=build_frontier(sequence,q,1,rng.randrange(2,n+1))
        base=thermal_response_bounds(cert)
        improved=thermal_response_bounds(cert,use_saw_count_cap=True,
                                          use_monotone_spectra=True)
        assert improved.partition_lower<=z<=improved.partition_upper
        assert improved.partition_upper<=base.partition_upper
        assert improved.variance_lower>=base.variance_lower
        assert improved.variance_upper<=base.variance_upper


def test_length_100_strict_reduction_and_nonzero_fluctuations():
    cert=build_frontier('H'*100,Q(1),1,2)
    base=thermal_response_bounds(cert,use_monotone_spectra=True)
    improved=thermal_response_bounds(cert,use_monotone_spectra=True,
                                      use_saw_count_cap=True)
    assert improved.partition_upper<base.partition_upper
    assert improved.variance_lower>base.variance_lower>0
    assert improved.event_susceptibility_lower>base.event_susceptibility_lower>0
    assert improved.variance_upper<=base.variance_upper
