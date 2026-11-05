"""Hermetic tests for prefix-cover and source-independent bound composition."""
import random
import pytest
from fractions import Fraction
from gnm_prefix_tree_independent_replay import verify_prefix_cover


def test_full_depth_coverage():
    for n in range(1,8):
        x=[(n,m) for m in range(1<<n)]
        assert verify_prefix_cover(x,n)==list(range(1<<n))


def test_variable_depth_coverage():
    # One leaf with first bit zero; two leaves with first bit one.
    # 001, 011, 101, 111 route by their low-order prefix bits.
    x=[(1,0),(2,1),(2,3)]
    assert verify_prefix_cover(x,6)==[0 if m%2==0 else (1 if m%4==1 else 2)
                                      for m in range(64)]


def test_duplicate_overlap_and_missing_fail_closed():
    invalid=[[(0,0),(1,0)],[(1,0)],[(1,0),(1,0)],[(1,0),(2,1)],
             [(2,0),(2,1),(2,2)],[(1,2)],[(7,1)],[(2,-1)],[(1,True)]]
    for x in invalid:
        with pytest.raises(ValueError):verify_prefix_cover(x,6)


def test_random_prefix_trees_cover_every_world_once():
    rng=random.Random(13579)
    for _ in range(100):
        def tree(depth,mask):
            if depth==6 or (depth>0 and rng.random()<.45):
                return [(depth,mask)]
            return tree(depth+1,mask)+tree(depth+1,mask|(1<<depth))
        x=tree(0,0)
        rng.shuffle(x)
        owners=verify_prefix_cover(x,6)
        for world,index in enumerate(owners):
            d,m=x[index]
            assert world&((1<<d)-1)==m


def test_leaf_certification_composes_for_all_boolean_networks():
    rng=random.Random(918)
    for _ in range(60):
        tree=[]
        def make(d,m):
            if d==5 or rng.random()<.25:
                tree.append((d,m))
            else:
                make(d+1,m); make(d+1,m|(1<<d))
        make(0,0)
        owners=verify_prefix_cover(tree,5)
        # Toy monotone contact response is bounded by its case extreme edges:
        # a graph with optional edges S has q=1/(1+sum_{e in S} weights[e]).
        weights=[Fraction(rng.randint(1,7),10) for _ in range(5)]
        for world,index in enumerate(owners):
            depth,mask=tree[index]
            forced_on=sum(weights[e] for e in range(depth) if mask>>e&1)
            remaining=sum(weights[e] for e in range(depth,5))
            actual=sum(weights[e] for e in range(5) if world>>e&1)
            assert Fraction(1,1)/(1+forced_on+remaining) <= Fraction(1,1)/(1+actual)
            assert Fraction(1,1)/(1+actual) <= Fraction(1,1)/(1+forced_on)