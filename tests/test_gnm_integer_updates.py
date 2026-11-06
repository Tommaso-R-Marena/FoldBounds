from dataclasses import replace
from fractions import Fraction
from itertools import combinations
from random import Random

import pytest
pytest.importorskip('sympy')
from foldbounds.gnm_joint_rank import untrusted_propose_inverse,verify_inverse
from foldbounds.gnm_integer_updates import propose_add_contact,ExactEndpointCache


def test_exact_contact_inverse_rank_one_equal_independent_direct_matrix():
    rng=Random(915)
    for n in (4,5,6,8):
        chain={(i,i+1) for i in range(n-1)}
        optional=[(i,j) for i in range(n) for j in range(i+2,n)]
        for _ in range(10):
            rng.shuffle(optional)
            initial=tuple(sorted(chain|set(optional[:min(2,len(optional))])))
            base=untrusted_propose_inverse(n,initial)
            edge=next(e for e in optional if e not in initial)
            proposed=propose_add_contact(n,initial,base,edge)
            updated=tuple(sorted(set(initial)|{edge}))
            direct=untrusted_propose_inverse(n,updated)
            assert proposed.denominator>0
            assert verify_inverse(n,updated,proposed)
            for i in range(n):
                for j in range(n):
                    assert (proposed.numerator[i][j]*direct.denominator==
                            direct.numerator[i][j]*proposed.denominator)


def test_cache_and_tampering_fail_closed():
    n=7
    chain=tuple((i,i+1) for i in range(n-1))
    c=ExactEndpointCache(n,chain)
    assert c.count==1
    enriched=tuple(sorted(set(chain)|{(0,3),(0,5),(2,6)}))
    a=c.get(enriched)
    assert verify_inverse(n,enriched,a)
    count=c.count
    assert c.get(enriched) is a
    assert c.count==count
    with pytest.raises(ValueError):
        c.get(((0,3),))
    with pytest.raises(ValueError):
        propose_add_contact(n,chain,c.get(chain),(0,1))
    corrupted=replace(c.get(chain),denominator=c.get(chain).denominator+1)
    with pytest.raises(ValueError):
        propose_add_contact(n,chain,corrupted,(0,3))
