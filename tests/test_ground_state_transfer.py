"""Independent exact small-grid oracles and fail-closed ground-state regression."""
from dataclasses import replace
from fractions import Fraction
import os
import pytest

from foldbounds.ground_state_transfer import (
    count_from_start, build_ground_state_certificate,
    verify_ground_state_certificate, certify_ground_population,
    directed_ground_union,
)
from foldbounds.partial_directed import is_partially_directed_path

MOVES=((1,0),(0,1),(-1,0),(0,-1))


def _direct_rectangle_paths(w,h):
    """Direct geometry-first exhaustive DFS; no DP recurrence or helpers."""
    cells={(x,y) for x in range(w) for y in range(h)}
    results=[]
    for sx,sy in cells:
        if (sx+1,sy) not in cells:
            continue
        path=[(sx,sy),(sx+1,sy)]
        used=set(path)
        def dfs():
            if len(path)==w*h:
                results.append(tuple((x-sx,y-sy) for x,y in path))
                return
            x,y=path[-1]
            for dx,dy in MOVES:
                nxt=(x+dx,y+dy)
                if nxt in cells and nxt not in used:
                    used.add(nxt);path.append(nxt)
                    dfs()
                    path.pop();used.remove(nxt)
        dfs()
    return results


def _independent_contact_score(path):
    return sum(
        abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
        for i in range(len(path)) for j in range(i+2,len(path))
    )


def _independent_all_saws(n):
    path=[(0,0),(1,0)]
    used=set(path)
    result=[]
    def dfs():
        if len(path)==n:
            result.append(tuple(path))
            return
        x,y=path[-1]
        for dx,dy in MOVES:
            nxt=(x+dx,y+dy)
            if nxt not in used:
                used.add(nxt);path.append(nxt)
                dfs()
                path.pop();used.remove(nxt)
    dfs()
    return result


@pytest.mark.parametrize("w,h",[(2,2),(2,3),(3,3),(4,3),(4,4)])
def test_direct_dfs_vs_transfer_for_all_starts(w,h):
    paths=_direct_rectangle_paths(w,h)
    by_start={}
    for path in paths:
        origin=(-min(x for x,y in path),-min(y for x,y in path))
        by_start[origin]=by_start.get(origin,0)+1
    for y in range(h):
        for x in range(w-1):
            assert count_from_start(w,h,x,y)==by_start.get((x,y),0)
    if w==h:
        cert=build_ground_state_certificate(w)
        assert cert.degeneracy==len(paths)
        assert verify_ground_state_certificate(cert)


def test_square_counts_vs_independent_published_oeis_a096969():
    # OEIS A096969 oriented Hamiltonian totals divided by four directions.
    for m,expected in ((2,2),(3,10),(4,138),(5,2162),
                       (6,114674),(7,6767640)):
        cert=build_ground_state_certificate(m)
        assert cert.max_contacts==(m-1)**2
        assert cert.degeneracy==expected
        assert verify_ground_state_certificate(cert)


def test_exact_global_coefficient_small():
    for side in (2,3):
        cert=build_ground_state_certificate(side)
        all_paths=_independent_all_saws(side*side)
        scores=[_independent_contact_score(p) for p in all_paths]
        assert max(scores)==cert.max_contacts
        assert scores.count(cert.max_contacts)==cert.degeneracy


def test_disjoint_directed_and_ground_union_small():
    for side in (2,3):
        cert=build_ground_state_certificate(side)
        all_paths=_independent_all_saws(side*side)
        ground={p for p in all_paths if _independent_contact_score(p)==cert.max_contacts}
        directed={p for p in all_paths if is_partially_directed_path(p)}
        for q in (Fraction(1),Fraction(3,2),Fraction(10)):
            combo=directed_ground_union(cert,q)
            assert combo.directed_overlap==len(ground & directed)
            expected=sum((q**_independent_contact_score(p) for p in ground|directed),Fraction(0))
            assert combo.union_lower==expected


def test_gibbs_population_enclosures_against_all_short_walks():
    for side in (2,3):
        cert=build_ground_state_certificate(side)
        paths=_independent_all_saws(side*side)
        hist={}
        for path in paths:
            k=_independent_contact_score(path)
            hist[k]=hist.get(k,0)+1
        for q in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(5)):
            z=sum((c*q**k for k,c in hist.items()),Fraction(0))
            for t in range(1,cert.max_contacts+1):
                event=sum((c*q**k for k,c in hist.items() if k>=t),Fraction(0))
                p=certify_ground_population(cert,q,t)
                assert p.probability_lower<=event/z
                robust=certify_ground_population(cert,q,t,
                    correction_lower=Fraction(1,2),correction_upper=Fraction(2))
                adversarial=Fraction(1,2)*event/(Fraction(1,2)*event+2*(z-event))
                assert robust.probability_lower<=adversarial


def test_all_real_q_tail_gate_monotonic_for_sampled_q():
    cert=build_ground_state_certificate(4)
    for t in (1,3,6,9):
        values=[certify_ground_population(cert,q,t).probability_lower
                for q in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(5),Fraction(10))]
        assert values==sorted(values)


def test_certificate_tampering_and_bad_arguments_rejected():
    c=build_ground_state_certificate(4)
    assert not verify_ground_state_certificate(replace(c,degeneracy=c.degeneracy+1))
    assert not verify_ground_state_certificate(replace(c,max_contacts=8))
    assert not verify_ground_state_certificate(replace(c,per_start=c.per_start[1:]))
    assert not verify_ground_state_certificate(replace(c,side=True))
    assert not verify_ground_state_certificate(replace(c,per_start=((True,0,1),)))
    with pytest.raises(ValueError):build_ground_state_certificate(True)
    with pytest.raises(ValueError):build_ground_state_certificate(11)
    with pytest.raises(ValueError):certify_ground_population(c,Fraction(1,2),1)
    with pytest.raises(ValueError):certify_ground_population(c,Fraction(2),0)
    with pytest.raises(ValueError):certify_ground_population(c,Fraction(2),10)
    with pytest.raises(ValueError):certify_ground_population(c,Fraction(2),2,
        correction_lower=Fraction(3),correction_upper=Fraction(1))
    with pytest.raises(ValueError):directed_ground_union(c,Fraction(1,2))


@pytest.mark.skipif(os.environ.get("FOLDBOUNDS_LONG_TESTS")!="1",
    reason="expensive 8x8 full transfer replay; opt in explicitly")
def test_64_residue_complete_coefficient_and_robust_majority():
    c=build_ground_state_certificate(8)
    assert c.degeneracy==1511656642
    assert verify_ground_state_certificate(c)
    z=directed_ground_union(c,Fraction(10))
    assert z.directed_overlap==2
    assert z.union_lower/z.directed_lower>1000000
    ordinary=certify_ground_population(c,Fraction(10),30)
    assert ordinary.walk_count_cap==5231553077432424492299000247
    assert ordinary.probability_lower>Fraction(96,100)
    robust=certify_ground_population(c,Fraction(10),30,
        correction_lower=Fraction(1,2),correction_upper=Fraction(2))
    assert robust.probability_lower>Fraction(87,100)
