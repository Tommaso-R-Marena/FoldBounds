"""Independent synthetic oracles for mutation-invariant geometrical code."""
from dataclasses import replace
from fractions import Fraction
from itertools import combinations, product

import pytest

from foldbounds.ga_gb_nmr import GA95, GB95, NMRModel
from foldbounds.ga_gb_error_codes import Reporter, matching, extract_reporters
from foldbounds.ga_gb_mutation_safe import (
    FORBIDDEN_MUTATION_SITES, MutationSafeCode,
    make_safe_code, check_safe_code, reporter_midpoint,
    channel_vote, decode_with_erasures, certified_budget,
    coordinate_distance_readings, finite_margin_to_truth,
)


def model_bank(seq,gb=False):
    answer=[]
    for k in range(20):
        coords={}
        for i in range(1,57):
            # Nondegenerate exact milli-Angstrom coordinates.
            x=round(1.1*i,3)
            y=round(k/100+((23. if gb else 0.) if i>14 else 0.),3)
            coords[i]=(x,y,0.)
        answer.append(NMRModel(seq,coords,{}))
    return tuple(answer)


@pytest.fixture(scope="module")
def banks():
    return model_bank(GA95),model_bank(GB95,True)


def test_conserved_policy_excludes_mutation_25_and_original_sites(banks):
    a,b=banks
    assert FORBIDDEN_MUTATION_SITES=={20,25,30,45}
    edges=extract_reporters(a,b)
    assert any(e.a==25 or e.b==25 for e in edges)
    for k in (1,3,5,7):
        code=make_safe_code(a,b,k)
        assert check_safe_code(code,a,b)
        assert all(e.a not in FORBIDDEN_MUTATION_SITES
                   and e.b not in FORBIDDEN_MUTATION_SITES
                   for e in code.reference.reporters)
        assert code.eligible_edges<len(edges)
        assert len({p for e in code.reference.reporters for p in e.endpoints})==2*k


def test_verifiable_maximum_matching_threshold_on_predeclared_graph(banks):
    a,b=banks
    for k in (1,3,5,7):
        code=make_safe_code(a,b,k)
        allowed=[e for e in extract_reporters(a,b)
                 if not ({e.a,e.b} & FORBIDDEN_MUTATION_SITES)]
        best=code.reference.bottleneck_milliang
        assert len(matching(allowed,best))>=k
        assert len(matching(allowed,best+1))<k
        left=set(code.reference.stronger_graph_cover_left)
        right=set(code.reference.stronger_graph_cover_right)
        assert len(left)+len(right)<k
        assert all(e.a in left or e.b in right
                   for e in allowed if e.gap_milliang>best)


def test_reject_forged_mutation_sites_and_strong_margin_certificate(banks):
    a,b=banks
    c=make_safe_code(a,b,3)
    unsafe=replace(c,forbidden_sites=frozenset({20,30,45}))
    with pytest.raises(ValueError):
        check_safe_code(unsafe,a,b)
    forged=replace(c,reference=replace(c.reference,bottleneck_milliang=10**10))
    with pytest.raises(ValueError):
        check_safe_code(forged,a,b)
    forged=replace(c,reference=replace(c.reference,stronger_graph_cover_left=()))
    # Some feasible matching optima have all stronger edges on right;
    # force a universally invalid candidate count instead.
    forged=replace(c,eligible_edges=0)
    with pytest.raises(ValueError):
        check_safe_code(forged,a,b)
    with pytest.raises(ValueError):
        make_safe_code(a,b,2)
    with pytest.raises(ValueError):
        make_safe_code(a,b,1001)


@pytest.mark.parametrize("k", (1,3,5,7,9))
def test_error_erasure_budget_is_sharp_for_all_k(k):
    for e in range(k+1):
        for t in range(k+1):
            should=(e+t<=k and 2*t+e<k)
            assert certified_budget(k,e,t)==should
            correct_votes=k-e-t
            if should:
                assert correct_votes>t
            elif e+t<=k:
                assert correct_votes<=t


def test_all_budget_allowable_erasure_error_patterns_independently():
    for k in (1,3,5,7,9):
        reports=tuple(Reporter(i,100+i,10_000,"GA95",
                               (10_000,12_000),(22_000,25_000))
                      for i in range(1,k+1))
        from foldbounds.ga_gb_error_codes import CodeCertificate
        code=MutationSafeCode(
            CodeCertificate(reports,k,(k-1)//2,10_000,(),(),k,
                            {"2KDL":"pinned","2KDM":"pinned"},14,20,
                            Fraction(5000)),
            FORBIDDEN_MUTATION_SITES,k
        )
        full={r.endpoints:Fraction(11_000) for r in reports}
        for e in range(k+1):
            for absent in combinations(reports,e):
                remaining=[r for r in reports if r not in absent]
                for t in range(len(remaining)+1):
                    if not certified_budget(k,e,t):
                        continue
                    for wrong in combinations(remaining,t):
                        obs=dict(full)
                        for r in absent:
                            obs.pop(r.endpoints)
                        for r in wrong:
                            obs[r.endpoints]=Fraction(100_000)
                        assert decode_with_erasures(code,obs)=="GA95"
        for r in reports:
            missing=frozenset({r.a})
            result=decode_with_erasures(code,full,missing_sites=missing)
            if k>1:
                assert result=="GA95"
            else:
                assert result is None


def test_unknown_channels_inexact_readouts_and_ties_abstain():
    r=Reporter(1,9,10_000,"GA95",(10_000,12_000),(22_000,25_000))
    from foldbounds.ga_gb_error_codes import CodeCertificate
    cert=CodeCertificate((r,),1,0,10_000,(),(),1,{},14,20,Fraction(5000))
    code=MutationSafeCode(cert,FORBIDDEN_MUTATION_SITES,1)
    assert decode_with_erasures(code,{}) is None
    assert decode_with_erasures(code,{r.endpoints:None}) is None
    assert decode_with_erasures(code,{r.endpoints:Fraction(17_000)}) is None
    with pytest.raises(ValueError):
        decode_with_erasures(code,{(2,11):Fraction(0)})
    with pytest.raises(ValueError):
        decode_with_erasures(code,{r.endpoints:18_000.0})
    with pytest.raises(ValueError):
        decode_with_erasures(code,{r.endpoints:Fraction(-1)})
    with pytest.raises(ValueError):
        certified_budget(True,0,0)
    assert reporter_midpoint(r)==Fraction(17_000)
    assert channel_vote(r,Fraction(10_000))=="GA95"
    assert channel_vote(r,Fraction(25_000))=="GB95"


def test_integer_bounds_on_conformer_distances(banks):
    a,b=banks
    c=make_safe_code(a,b,3)
    for expected,models in (("GA95",a),("GB95",b)):
        for m in models:
            readings=coordinate_distance_readings(c,m)
            assert decode_with_erasures(c,readings)==expected
            assert finite_margin_to_truth(c,m,expected)>0
    with pytest.raises(ValueError):
        finite_margin_to_truth(c,a[0],"OtherFold")
