"""Full thermodynamic regression for disjoint column-progressive witnesses.

Exact reference enumerator does not use production contact-range or scorer.
"""
from fractions import Fraction as Q
from itertools import product
from random import Random

import pytest

from foldbounds.lattice import build_frontier
from foldbounds.partial_directed import partially_directed_spectrum
from foldbounds.thermal_response import thermal_response_bounds


def walks(n):
    all_walks=[((0,0),(1,0))]
    for _ in range(2,n):
        next_walks=[]
        for walk in all_walks:
            x,y=walk[-1]
            for site in ((x,y+1),(x,y-1),(x+1,y),(x-1,y)):
                if site not in walk:
                    next_walks.append(walk+(site,))
        all_walks=next_walks
    return all_walks


def score(seq,path):
    return sum(seq[i]=='H' and seq[j]=='H' and
               abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
               for i in range(len(seq)) for j in range(i+2,len(seq)))


def observables(seq,q,threshold):
    z=Q(0);z1=Q(0);z2=Q(0);a=Q(0);a1=Q(0)
    for path in walks(len(seq)):
        k=score(seq,path)
        w=q**k
        z+=w;z1+=w*k;z2+=w*k*k
        if k>=threshold:
            a+=w;a1+=w*k
    mean=z1/z
    return z,mean,z2/z,z2/z-mean**2,(a1*z-a*z1)/z**2,a/z


def check_bounds(ans,exact):
    z,mean,m2,var,susc,prob=exact
    assert ans.partition_lower<=z<=ans.partition_upper
    assert ans.mean_lower<=mean<=ans.mean_upper
    assert ans.second_moment_lower<=m2<=ans.second_moment_upper
    assert ans.variance_lower<=var<=ans.variance_upper
    assert ans.event_susceptibility_lower<=susc<=ans.event_susceptibility_upper
    assert ans.probability_lower<=prob<=ans.probability_upper


def test_exhaustive_short_sequences_validated_against_full_model():
    for n in range(2,7):
        for seq_tuple in product('HP',repeat=n):
            seq=''.join(seq_tuple)
            for q in (Q(1),Q(3,2)):
                exact=observables(seq,q,1)
                for depth in (2,n):
                    cert=build_frontier(seq,q,1,depth)
                    ans=thermal_response_bounds(cert,
                        use_partially_directed_spectra=True,
                        use_geometric_contact_cap=True,
                        use_saw_count_cap=True)
                    check_bounds(ans,exact)


def test_randomized_mixed_frontiers_thresholds_and_many_witnesses():
    rng=Random(92182)
    for _ in range(70):
        n=rng.randrange(3,10)
        seq=''.join(rng.choice('HP') for _ in range(n))
        q=rng.choice((Q(1), Q(3,2), Q(2), Q(5,3)))
        threshold=rng.randrange(0,4)
        depth=rng.randrange(2,n+1)
        cert=build_frontier(seq,q,threshold,depth)
        allwalks=walks(n)
        # All terminal witnesses whose prefixes fall in an unknown subtree.
        bounded={node.path for node in cert.frontier if node.kind=='bound' and len(node.path)<n}
        eligible=[w for w in allwalks if any(w[:len(pre)]==pre for pre in bounded)]
        rng.shuffle(eligible)
        witnesses=tuple(eligible[:rng.randrange(0,min(12,len(eligible))+1)])
        ans=thermal_response_bounds(cert,witnesses=witnesses,
            use_partially_directed_spectra=True,
            use_saw_count_cap=True,use_geometric_contact_cap=True)
        check_bounds(ans,observables(seq,q,threshold))
        assert ans.witness_paths==len(witnesses)


def test_compact_geometry_and_directed_family_provide_strict_improvements():
    for n in (12,20,30,40):
        cert=build_frontier('H'*n,Q(1),1,2)
        cone=thermal_response_bounds(cert,use_monotone_spectra=True,
                                     use_saw_count_cap=True)
        directed=thermal_response_bounds(cert,use_partially_directed_spectra=True,
                                         use_saw_count_cap=True,
                                         use_geometric_contact_cap=True)
        assert directed.partially_directed_paths>cone.monotone_paths
        assert directed.partition_lower>cone.partition_lower
        assert directed.partition_upper<=cone.partition_upper
        assert directed.variance_lower>cone.variance_lower
        assert directed.variance_upper < cone.variance_upper
        assert directed.mean_upper < cone.mean_upper
        assert directed.support_upper <= n
        assert directed.partially_directed_subtrees>=1
        assert directed.geometric_contact_cap_enabled


def test_invalid_overlap_modes_and_incorrect_options_fail_closed():
    cert=build_frontier('HHHHH',Q(2),1,2)
    with pytest.raises(ValueError,match='overlap'):
        thermal_response_bounds(cert,use_monotone_spectra=True,
                                use_partially_directed_spectra=True)
    for kwargs in ({'use_partially_directed_spectra':1},
                   {'use_geometric_contact_cap':1}):
        with pytest.raises(ValueError):
            thermal_response_bounds(cert,**kwargs)


def test_partition_into_disjoint_frontier_prefixes_no_duplicate_masses():
    rng=Random(808)
    for n in range(3,9):
        seq=''.join(rng.choice('HP') for _ in range(n))
        full=dict(partially_directed_spectrum(seq))
        for depth in range(2,n+1):
            cert=build_frontier(seq,Q(1),1,depth)
            aggregate={}
            for node in cert.frontier:
                if node.kind=='dead':
                    continue
                hist=dict(partially_directed_spectrum(seq,node.path))
                for k,count in hist.items():
                    aggregate[k]=aggregate.get(k,0)+count
            assert aggregate==full
