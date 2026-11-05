"""Independent all-walk oracle and tampering suite for ranked shape envelopes."""
from fractions import Fraction
from random import Random
from dataclasses import replace
import pytest
from foldbounds.shape_rank_upper import (build_shape_rank_certificate,
    verify_shape_rank_certificate, ShapeClass)
from foldbounds.horizon_certificate import ROOT, _children
from foldbounds.saw_count_cap import optimal_saw_count_cap


def walks(n):
    w=[ROOT]
    for _ in range(2,n):w=[v for p in w for v in _children(p)]
    return w

def score(seq,w):
    return sum(seq[i]==seq[j]=='H' and j-i>=2 and abs(w[i][0]-w[j][0])+abs(w[i][1]-w[j][1])==1
               for i in range(len(w)) for j in range(i+1,len(w)))

@pytest.mark.parametrize('n',range(2,11))
def test_independent_exact_population_enclosed(n):
    allwalks=walks(n)
    for seq in ('H'*n,'P'*n,('HPPH'*n)[:n],('HP'*n)[:n]):
        prefix=allwalks[len(allwalks)//2][:min(n,5)]
        for p in (ROOT,prefix):
            cases=[w for w in allwalks if w[:len(p)]==p]
            for memory in (2,4):
                cert=build_shape_rank_certificate(seq,p,suffix_memory=memory)
                for q in (Fraction(1),Fraction(3,2),Fraction(2)):
                    res=verify_shape_rank_certificate(cert,q,1)
                    z=sum((q**score(seq,w) for w in cases),Fraction())
                    event=sum((q**score(seq,w) for w in cases if score(seq,w)>=1),Fraction())
                    assert z<=res.partition_upper
                    assert event<=res.event_upper
                    assert z-event<=res.complement_upper


def test_random_and_tamper():
    rng=Random(2206)
    for _ in range(40):
        n=rng.randrange(5,11);seq=''.join(rng.choice('HP') for _ in range(n))
        p=rng.choice(walks(n))[:rng.randrange(2,n+1)]
        cert=build_shape_rank_certificate(seq,p,suffix_memory=4)
        b=verify_shape_rank_certificate(cert,Fraction(2))
        z=sum((Fraction(2)**score(seq,w) for w in walks(n) if w[:len(p)]==p),Fraction())
        assert z<=b.partition_upper
        if cert.histogram:
            changed=replace(cert,histogram=(replace(cert.histogram[0],relaxed_count=cert.histogram[0].relaxed_count+1),)+cert.histogram[1:])
            with pytest.raises(ValueError):verify_shape_rank_certificate(changed,Fraction(2))
        with pytest.raises(ValueError):verify_shape_rank_certificate(replace(cert,true_count_upper=0),Fraction(2))


def test_ranked_shape_upper_never_weaker_than_uniform_count_box():
    from foldbounds.horizon_certificate import _score_interval
    from foldbounds.saw_count_cap import optimal_saw_count_cap
    for n in (8,12,16,20):
        sequence='H'*n
        cap=_score_interval(sequence,ROOT)[1]
        count=optimal_saw_count_cap(n-2).count_upper
        cert=build_shape_rank_certificate(sequence,suffix_memory=4)
        for q in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(3)):
            result=verify_shape_rank_certificate(cert,q)
            assert result.partition_upper <= count*q**cap
            assert result.event_upper <= result.partition_upper
            assert result.complement_upper <= result.partition_upper


def test_parity_sparse_patterns_and_higher_thresholds():
    from foldbounds.horizon_certificate import _score_interval
    for seq in ('HPHPHPHPH','HPPHHPPHH','HPPPPPPPH'):
        w=walks(len(seq))
        for t in (0,1,2,4):
            c=build_shape_rank_certificate(seq,suffix_memory=4)
            b=verify_shape_rank_certificate(c,Fraction(3,2),t)
            z=sum((Fraction(3,2)**score(seq,v) for v in w),Fraction())
            e=sum((Fraction(3,2)**score(seq,v) for v in w if score(seq,v)>=t),Fraction())
            assert b.partition_upper>=z
            assert b.event_upper>=e
            assert b.complement_upper>=z-e


def test_separate_exhaustive_relaxed_trajectory_histograms():
    """Independently enumerate finite-memory *invalid* walks, not only true SAWs."""
    from collections import Counter
    from foldbounds.shape_rank_upper import ShapeClass
    moves=((1,0),(0,1),(-1,0),(0,-1))
    for n in range(3,10):
        true_walks=walks(n)
        for seq in ('H'*n,('HPPH'*n)[:n]):
            potential=sum(seq[i]==seq[j]=='H' and (j-i)%2==1
                          for i in range(n) for j in range(i+3,n))
            for prefix in (ROOT,true_walks[len(true_walks)//3][:min(5,n)]):
                used_prefix=set(prefix)
                fixed=0
                potential_suffix=0
                for i in range(n):
                    if seq[i]!='H':continue
                    for j in range(i+3,n,2):
                        if seq[j]!='H':continue
                        if j>=len(prefix):potential_suffix+=1
                        elif (abs(prefix[i][0]-prefix[j][0])+
                              abs(prefix[i][1]-prefix[j][1])==1):fixed+=1
                for memory in (2,4):
                    hist=Counter()
                    def visit(path):
                        if len(path)==n:
                            xs=[p[0] for p in path];ys=[p[1] for p in path]
                            w=max(xs)-min(xs)+1;h=max(ys)-min(ys)+1
                            cap=min(fixed+potential_suffix,max(0,n+1-w-h))
                            if w*h>=n and cap>=fixed:
                                hist[cap]+=1
                            return
                        x,y=path[-1]
                        for dx,dy in moves:
                            dest=(x+dx,y+dy)
                            if dest in used_prefix or dest in path[-memory:]:
                                continue
                            visit(path+(dest,))
                    visit(prefix)
                    expected=tuple(ShapeClass(k,hist[k]) for k in sorted(hist,reverse=True))
                    got=build_shape_rank_certificate(seq,prefix,suffix_memory=memory)
                    assert got.histogram==expected
