"""Independent exhaustive-oracle tests for weighted monotone conformation families."""
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.lattice import CertificateError, ROOT, build_frontier
from foldbounds.lookahead import verify_lookahead
from foldbounds.monotone import monotone_family_mass
from foldbounds.adaptive import adaptive_refine, verify_refinement_run

DIRECTIONS=((1,0),(-1,0),(0,1),(0,-1))
QUADRANTS=(
    {(1,0),(0,1)}, {(1,0),(0,-1)},
    {(-1,0),(0,1)}, {(-1,0),(0,-1)},
)


def independent_walks(n):
    paths=[ROOT]
    for _ in range(2,n):
        nxt=[]
        for path in paths:
            x,y=path[-1]
            for dx,dy in DIRECTIONS:
                new=(x+dx,y+dy)
                if new not in path:
                    nxt.append(path+(new,))
        paths=nxt
    return paths


def independent_contacts(seq,path):
    return sum(seq[i]=='H' and seq[j]=='H'
               and abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
               for i in range(len(seq)) for j in range(i+2,len(seq)))


def monotone_suffix(path,m):
    moves={tuple(path[j][k]-path[j-1][k] for k in range(2))
           for j in range(m,len(path))}
    return any(moves<=cone for cone in QUADRANTS)


def independent_masses(seq,prefix,q,threshold,complete):
    matching=(walk for walk in complete if walk[:len(prefix)]==prefix
              and monotone_suffix(walk,len(prefix)))
    event=other=Fraction(0)
    count=0
    for walk in matching:
        k=independent_contacts(seq,walk)
        w=q**k
        count+=1
        if k>=threshold: event+=w
        else: other+=w
    return event,other,count


@pytest.mark.parametrize('n',range(2,9))
def test_all_restricted_completions_against_independent_exhaustive_oracle(n):
    walks=independent_walks(n)
    sequences=['H'*n,'P'*n,('HP'*n)[:n],('PH'*n)[:n],('HHPHHPPH'*n)[:n]]
    if n<=5:
        sequences=[''.join(p) for p in product('HP',repeat=n)]
    for length in sorted({2,max(2,n//2),max(2,n-1),n}):
        prefixes=sorted({walk[:length] for walk in walks})
        if len(prefixes)>12:
            prefixes=prefixes[:6]+prefixes[len(prefixes)//2:len(prefixes)//2+3]+prefixes[-3:]
        for prefix in prefixes:
            for seq in sequences:
                for q in (Fraction(1),Fraction(3,2),Fraction(2)):
                    for t in (0,1,2,4):
                        actual=monotone_family_mass(seq,prefix,q,t)
                        e,b,c=independent_masses(seq,prefix,q,t,walks)
                        assert (actual.event,actual.complement,actual.distinct_paths)==(e,b,c)


@pytest.mark.parametrize('seq', ['HHHHHHHH','HPHHPHHP','HHPHHHPHHH'])
@pytest.mark.parametrize('q',[Fraction(1),Fraction(3,2),Fraction(2)])
@pytest.mark.parametrize('threshold',[0,1,2,3])
def test_monotone_option_cannot_weaken_existing_bounds(seq,q,threshold):
    cert=build_frontier(seq,q,threshold,2)
    old=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                         witness_nodes=0,use_monotone=False)
    new=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                         witness_nodes=0,use_monotone=True)
    assert old.partition_lower<=new.partition_lower
    assert new.partition_upper==old.partition_upper
    assert old.event_lower<=new.event_lower<=new.event_upper<=old.event_upper
    assert old.complement_lower<=new.complement_lower<=new.complement_upper<=old.complement_upper
    assert old.probability_lower<=new.probability_lower<=new.probability_upper<=old.probability_upper


def test_nontrivial_lower_bound_improvement():
    cert=build_frontier('H'*12,Fraction(2),1,2)
    old=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                         witness_nodes=0)
    new=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                         witness_nodes=0,use_monotone=True)
    assert new.partition_lower>old.partition_lower
    assert new.event_lower>old.event_lower
    assert new.partition_upper==old.partition_upper


def test_blocked_prefix_gives_zero():
    prefix=((0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2),(0,1),(1,1))
    m=monotone_family_mass('H'*10,prefix,Fraction(2),1)
    assert m.distinct_paths==0
    assert m.event==m.complement==0


def test_adaptive_replay_records_mode():
    cert=build_frontier('H'*10,Fraction(2),1,2)
    run=adaptive_refine(cert,steps=3,lookahead=0,use_matching=True,
                        use_transfer=False,witness_nodes=0,
                        use_monotone=True,policy='mass')
    assert run.use_monotone is True
    assert verify_refinement_run(run)==run.bounds
    with pytest.raises(ValueError):
        verify_refinement_run(replace(run,use_monotone=False))


@pytest.mark.parametrize('bad',[True,-1,1.5,'yes'])
def test_invalid_threshold_rejected(bad):
    if type(bad) is int:
        # Negative thresholds are legitimate: event is the entire space.
        assert monotone_family_mass('HHHH',ROOT,Fraction(2),bad).complement==0
    else:
        with pytest.raises(ValueError):
            monotone_family_mass('HHHH',ROOT,Fraction(2),bad)


@pytest.mark.parametrize('bad',[Fraction(0),Fraction(1,2),2.0])
def test_invalid_q_rejected(bad):
    with pytest.raises(ValueError):
        monotone_family_mass('HHHH',ROOT,bad,1)


def test_invalid_geometry_rejected():
    with pytest.raises(CertificateError):
        monotone_family_mass('HHHH',ROOT+((0,0),),Fraction(2),1)


@pytest.mark.parametrize('bad',[0,1,'yes'])
def test_boolean_mode_is_strict(bad):
    cert=build_frontier('HHHH',Fraction(2),1,2)
    with pytest.raises(ValueError):
        verify_lookahead(cert,use_monotone=bad)
