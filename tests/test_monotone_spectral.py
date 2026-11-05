"""Independent exact SAW oracle for disjoint compressed monotone families."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
import pytest

from foldbounds.lattice import build_frontier
from foldbounds.monotone_spectral import (
    monotone_spectral_histogram, is_monotone_cone_suffix,
)
from foldbounds.thermal_response import thermal_response_bounds


def walks(n):
    states = [((0,0),(1,0))]
    for _ in range(2,n):
        nxt = []
        for path in states:
            x,y = path[-1]
            for site in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)):
                if site not in path:
                    nxt.append(path+(site,))
        states = nxt
    return tuple(states)


def independent_score(seq, path):
    return sum(seq[i]==seq[j]=='H'
               and abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
               for i in range(len(seq)) for j in range(i+2,len(seq)))


def independent_cone(prefix, path):
    changes = set()
    for (x,y),(xx,yy) in zip(path[len(prefix)-1:],path[len(prefix):]):
        changes.add((xx-x,yy-y))
    return not ({(1,0),(-1,0)}<=changes or {(0,1),(0,-1)}<=changes)


def exact_observables(seq, q, threshold, all_paths):
    z=m1=m2=em=e=Q(0)
    for path in all_paths:
        score = independent_score(seq,path)
        mass = q**score
        z+=mass;m1+=mass*score;m2+=mass*score*score
        if score>=threshold:
            e+=mass
            em+=mass*score
    mu=m1/z
    return z,mu,m2/z-mu*mu,(em*z-e*m1)/z**2


def test_all_cone_histograms_vs_independently_filtered_full_paths():
    for n in range(3,9):
        all_paths=walks(n)
        prefixes={path[:d] for path in all_paths
                  for d in (2,min(4,n),min(6,n))}
        for seq in ('H'*n,'P'*n,('HP'*n)[:n]):
            for prefix in prefixes:
                expected=Counter(independent_score(seq,path) for path in all_paths
                                 if path[:len(prefix)]==prefix
                                 and (len(prefix)==n or independent_cone(prefix,path)))
                assert dict(monotone_spectral_histogram(seq,prefix))==dict(expected)
                if len(prefix)<n:
                    for path in all_paths:
                        if path[:len(prefix)]==prefix:
                            assert is_monotone_cone_suffix(prefix,path)==independent_cone(prefix,path)


def test_exhaustive_hp_thermodynamics_from_partial_cones():
    for n in range(2,8):
        all_paths=walks(n)
        for seq_tuple in product('HP',repeat=n):
            seq=''.join(seq_tuple)
            for q in (Q(1),Q(3,2)):
                z,mean,var,susc=exact_observables(seq,q,1,all_paths)
                for depth in (2,n):
                    result=thermal_response_bounds(build_frontier(seq,q,1,depth),
                                                   use_monotone_spectra=True)
                    assert result.partition_lower<=z<=result.partition_upper
                    assert result.mean_lower<=mean<=result.mean_upper
                    assert result.variance_lower<=var<=result.variance_upper
                    assert result.event_susceptibility_lower<=susc<=result.event_susceptibility_upper


def test_witness_cone_overlap_subtracted_exactly():
    seq='HHHHHHH'
    all_paths=walks(len(seq))
    cert=build_frontier(seq,Q(1),1,2)
    witness=all_paths[:60]
    base=thermal_response_bounds(cert,witnesses=witness)
    improved=thermal_response_bounds(cert,witnesses=witness,
                                     use_monotone_spectra=True)
    assert improved.partition_lower>=base.partition_lower
    assert improved.partition_upper==base.partition_upper
    assert improved.variance_lower>=base.variance_lower
    assert improved.variance_upper<=base.variance_upper
    assert improved.monotone_paths>0
    assert improved.witness_paths==len(witness)
    z,_,variance,_=exact_observables(seq,Q(1),1,all_paths)
    assert improved.partition_lower<=z<=improved.partition_upper
    assert improved.variance_lower<=variance<=improved.variance_upper


def test_length_30_positive_response_without_full_state_enumeration():
    cert=build_frontier('H'*30,Q(1),1,2)
    old=thermal_response_bounds(cert)
    improved=thermal_response_bounds(cert,use_monotone_spectra=True)
    assert improved.monotone_paths>2**28
    assert improved.variance_lower>0
    assert improved.event_susceptibility_lower>0
    assert improved.variance_upper<=old.variance_upper
    assert improved.mean_lower>=old.mean_lower
    assert improved.mean_upper<=old.mean_upper


def test_rejects_duplicate_or_invalid_witnesses():
    cert=build_frontier('HHHHHH',Q(2),1,2)
    path=walks(6)[0]
    for witnesses in ((path,path),(path[:-1],),(path[:-1]+((99,99),),)):
        with pytest.raises(ValueError):
            thermal_response_bounds(cert,witnesses=witnesses,
                                    use_monotone_spectra=True)
