"""Independent exhaustive graph-state audit of rank-two GNM certificates."""
from dataclasses import replace
from fractions import Fraction
from itertools import combinations, product
import random

import pytest
pytest.importorskip("sympy")

from foldbounds.gnm_joint_rank import (
    IntegerInverse, classify_contact_boxes, joint_rank_certificate,
    scaled_grounded_laplacian, untrusted_propose_inverse, verify_inverse,
)


def exact_diagonal(n, edges):
    # Separate test oracle: exact Fraction Gauss-Jordan; not the producer,
    # IntegerInverse class or the joint certificate arithmetic.
    a = scaled_grounded_laplacian(n, edges)
    t = [[Fraction(a[i][j]) for j in range(n)] +
         [Fraction(i == j) for j in range(n)] for i in range(n)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if t[i][k]), None)
        assert pivot is not None
        t[k], t[pivot] = t[pivot], t[k]
        z = t[k][k]
        t[k] = [x/z for x in t[k]]
        for i in range(n):
            if i != k:
                coef = t[i][k]
                t[i] = [u-coef*v for u, v in zip(t[i], t[k])]
    # GNM pseudoinverse diagonal = n*B^-1_ii - 1/n.
    return [n * t[i][n+i] - Fraction(1,n) for i in range(n)]


def test_exact_integer_inverse_rejected_after_adversarial_edits():
    n = 6
    e = tuple((i,i+1) for i in range(n-1)) + ((0,3),(1,5))
    proposed = untrusted_propose_inverse(n, e)
    assert verify_inverse(n, e, proposed)
    assert exact_diagonal(n, e)[0] == n * Fraction(proposed.numerator[0][0],proposed.denominator)-Fraction(1,n)
    broken = [list(row) for row in proposed.numerator]
    broken[0][0] += 1
    assert not verify_inverse(n, e, replace(proposed,numerator=tuple(tuple(row) for row in broken)))
    assert not verify_inverse(n, e, replace(proposed,denominator=-proposed.denominator))
    assert not verify_inverse(n, e, replace(proposed,numerator=proposed.numerator[:-1]))
    assert not verify_inverse(n, e, IntegerInverse(((1,),),1))


def test_every_provable_pair_order_holds_for_every_optional_edge_subset():
    randomizer = random.Random(314159)
    successes = 0
    improvements = 0
    for n in (5,6,7,8):
        backbone = {(i,i+1) for i in range(n-1)}
        available = [(i,j) for i in range(n) for j in range(i+2,n)]
        for repetition in range(16):
            randomizer.shuffle(available)
            mandatory = set(backbone)|set(available[:min(3,len(available))])
            optional = tuple(available[min(3,len(available)):min(6,len(available))])
            possible = mandatory | set(optional)
            lo = untrusted_propose_inverse(n,tuple(sorted(mandatory)))
            hi = untrusted_propose_inverse(n,tuple(sorted(possible)))
            cert = joint_rank_certificate(n,mandatory,possible,lo,hi)
            assert cert.simple_certified + cert.gained == cert.joint_certified
            for settings in product((False,True),repeat=len(optional)):
                active = mandatory | {edge for edge, present in zip(optional,settings) if present}
                diag = exact_diagonal(n,tuple(sorted(active)))
                for i,j,sign in cert.certified_orders:
                    assert (diag[i]-diag[j])*sign > 0
            successes+=1
            improvements+=cert.gained
    assert successes == 64
    assert improvements >= 1


def test_distance_boxes_conservative_and_source_only():
    coords=((0,0,0),(3800,0,0),(7600,0,0),(7600,3900,0),(0,0,7900))
    certain,possible=classify_contact_boxes(coords,Fraction(8),Fraction(1,20))
    assert set(certain) <= set(possible)
    assert set((i,i+1) for i in range(len(coords)-1)) <= set(certain)
    assert (0,4) in possible
    assert (0,4) not in certain
    with pytest.raises(ValueError):
        classify_contact_boxes(coords,Fraction(8),0.05)
    with pytest.raises(ValueError):
        classify_contact_boxes(coords,Fraction(8),Fraction(-1,4))


def test_fail_closed_inconsistent_edges_and_proofs():
    n=5
    lo={(i,i+1) for i in range(4)}
    hi=lo|{(0,3)}
    a=untrusted_propose_inverse(n,lo)
    b=untrusted_propose_inverse(n,hi)
    with pytest.raises(ValueError):
        joint_rank_certificate(n,hi,lo,b,a)
    with pytest.raises(ValueError):
        joint_rank_certificate(n,lo,hi,a,a)
    with pytest.raises(ValueError):
        joint_rank_certificate(n,lo,hi,a,b,exact_verification=False)
    with pytest.raises(ValueError):
        scaled_grounded_laplacian(n,lo|{(0,0)})
    with pytest.raises(ValueError):
        scaled_grounded_laplacian(n,{(0,4)})
