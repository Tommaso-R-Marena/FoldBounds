"""Independent adversarial validation of fixed-composition 3D HP switches."""
from dataclasses import replace
from fractions import Fraction as F
import pytest
from foldbounds.composition_swap import (DEFAULT,certify_swap,verify_swap,
    certify_landscape,verify_landscape,certify_error_budget,verify_error_budget,
    sequence)

D=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))

def independently_enumerate(parent,mutant):
    """Second full cubic DFS; no production geometry/fingerprint/energy code."""
    n=len(parent); t=2 if n==8 else 3
    path=[(0,0,0),(1,0,0)];visited=set(path);events=[0,0];totals=[0,0];count=0
    seqs=(parent,mutant)
    def walk():
        nonlocal count
        if len(path)==n:
            count+=1;pairs=[]
            for j in range(2,n):
                for i in range(j-1):
                    if sum(abs(path[i][k]-path[j][k]) for k in range(3))==1:
                        pairs.append((i,j))
            is_event=len(pairs)>=t
            for a,seq in enumerate(seqs):
                k=sum(seq[i]=='H' and seq[j]=='H' for i,j in pairs)
                w=8**k
                totals[a]+=w
                if is_event:events[a]+=w
            return
        x,y,z=path[-1]
        for dx,dy,dz in D:
            q=(x+dx,y+dy,z+dz)
            if q in visited:continue
            visited.add(q);path.append(q)
            walk()
            path.pop();visited.remove(q)
    walk()
    return count,(F(events[0],totals[0]),F(events[1],totals[1]))

@pytest.mark.parametrize('a,b,n,exact,err',[
    ('PPPHPHHH','HPPHPHPH',13565,(F(2516,7919),F(16152,22423)),F(19420473,11376575)),
    ('PPPHPHHHPH','HPPHPHPHPH',308981,(F(95124,332659),F(1264000,1997213)),F(2983192704,2031712291))])
def test_full_separate_oracle_and_error_tolerance(a,b,n,exact,err):
    cert=certify_swap(a,b)
    assert verify_swap(cert)
    assert cert.sites==(1,7) and cert.conformations==n
    assert a.count('H')==b.count('H')
    assert [sum(seq[i]=='H' for i in range(par,len(seq),2)) for par in (0,1) for seq in (a,b)][0::2]==[sum(seq[i]=='H' for i in range(par,len(seq),2)) for par in (0,1) for seq in (a,b)][1::2]
    assert cert.parent_robust_upper<F(1,2)<cert.mutant_robust_lower
    count,prob=independently_enumerate(a,b)
    assert count==n and prob==exact
    root=certify_error_budget(a,b)
    assert verify_error_budget(root) and root.max_symmetric_multiplier_squared==err
    assert err>F(36,25)
    assert not verify_error_budget(replace(root,max_symmetric_multiplier_squared=F(2)))

def test_full_mutation_hypercubes_no_duplicates_or_missing_edges():
    for n,edges,switches,orbits in [(8,768,12,6),(10,5120,4,2)]:
        result=certify_landscape(n)
        assert (result.sequences,result.allowed_edges,
                len(result.certified_switches),result.reversal_orbits)==(2**n,edges,switches,orbits)
        assert verify_landscape(result)
        assert all(verify_swap(certify_swap(sequence(a,n),sequence(b,n)))
                   for a,b,_,_ in result.certified_switches)

def test_forged_claims_and_wrong_composition_fail_closed():
    c=certify_swap('PPPHPHHHPH','HPPHPHPHPH')
    assert not verify_swap(replace(c,separation=F(100)))
    assert not verify_swap(replace(c,mutant_robust_lower=F(1)))
    assert not verify_swap(replace(c,sites=(2,7)))
    assert not verify_landscape(replace(certify_landscape(8),allowed_edges=10))
    for a,b in [('PPPHPHHH','HPPHPHHH'),
                ('PPPHPHHH','PPPHPHHH'),
                ('PHPHPHPH','HPPHPHPH'),
                ('PPPHPHHHPH','HPPHPHPH')]:
        with pytest.raises(ValueError):certify_swap(a,b)

def test_n10_error_bound_cannot_be_overextended():
    a,b='PPPHPHHHPH','HPPHPHPHPH'
    with pytest.raises(ValueError):
        certify_swap(a,b,(F(8),F(12),F(3,4),F(5,4),F(4,5),F(5,4)))

@pytest.mark.parametrize('bad',[
    (F(12),F(8),F(3,4),F(5,4),F(5,6),F(6,5)),
    (F(8),F(12),F(0),F(5,4),F(5,6),F(6,5)),
    (F(8),F(12),F(3,4),F(5,4),F(0),F(6,5)),
    (F(8),F(12),F(3,4),F(5,4),F(7,5),F(6,5))])
def test_bad_boxes_rejected(bad):
    with pytest.raises(ValueError):
        certify_swap('PPPHPHHH','HPPHPHPH',bad)
