"""Independent integer oracles for exact geometric order-statistic robustness."""
from dataclasses import replace
from fractions import Fraction
from itertools import combinations, product

import pytest

from foldbounds.ga_gb_nmr import GA95, NMRModel
from foldbounds.ga_gb_error_codes import CodeCertificate, Reporter
from foldbounds.ga_gb_mutation_safe import (
    MutationSafeCode,FORBIDDEN_MUTATION_SITES,decode_with_erasures
)
from foldbounds.mutation_safe_radius import (
    risk_channels_for_budget, sorted_signed_margins,
    certify_perturbation_radius,check_perturbation_certificate,
    lossless_radius_A,
)


def fake_code(margins_mA):
    k=len(margins_mA)
    edges=tuple(Reporter(i,15+i,18_000,"GA95",(10_000,11_000),
                         (29_000,30_000)) for i in range(1,k+1))
    cc=CodeCertificate(edges,k,(k-1)//2,18_000,(),(),k,{},14,20,
                       Fraction(9000))
    safe=MutationSafeCode(cc,FORBIDDEN_MUTATION_SITES,k)
    ca={i:(float(i),500.,0.) for i in range(1,57)}
    for i,margin in enumerate(margins_mA,1):
        a,b=i,15+i
        ca[a]=(0.,float(i*4),0.)
        ca[b]=((20_000-margin)/1000.,float(i*4),0.)
    return safe,NMRModel(GA95,ca,{})


def test_exact_quantile_radii_and_sharp_tradeoff():
    code,model=fake_code((1000,2000,3000,9000,9000))
    m=sorted_signed_margins(code,model,"GA95")
    assert m==tuple(map(Fraction,(1000,2000,3000,9000,9000)))
    data=(("GA95",model),)
    pristine=certify_perturbation_radius(code,data)
    assert pristine.geometry_risk_channels==2
    assert pristine.strict_pair_distance_radius_milliangstrom==3000
    assert pristine.strict_per_CA_site_radius_milliangstrom==1500
    assert pristine.holds_for_all_submitted_conformers
    assert check_perturbation_certificate(code,data,pristine)
    missing_and_wrong=certify_perturbation_radius(
        code,data,erased_channels=2,additional_wrong_channels=1)
    assert missing_and_wrong.geometry_risk_channels==0
    assert missing_and_wrong.strict_pair_distance_radius_milliangstrom==1000
    assert missing_and_wrong.strict_per_CA_site_radius_milliangstrom==500
    assert check_perturbation_certificate(code,data,missing_and_wrong)
    assert lossless_radius_A(Fraction(3000))=="3/1"
    with pytest.raises(ValueError):
        check_perturbation_certificate(
            code,data,replace(pristine,strict_pair_distance_radius_milliangstrom=Fraction(5000)))


def brute_force_majority_guarantee(margins,delta,e,t):
    """Independent nonoptimized exhaustive attacker over wrong/erased channels."""
    k=len(margins)
    for erased in combinations(range(k),e):
        survivors=[i for i in range(k) if i not in erased]
        for new_wrong in combinations(survivors,t):
            for noisy_mask in product((False,True),repeat=len(survivors)):
                wrong=set(new_wrong)
                for i,flip in zip(survivors,noisy_mask):
                    if flip and margins[i]<=delta:
                        wrong.add(i)
                votes_correct=len(survivors)-len(wrong)
                if votes_correct<=len(wrong):
                    return False
    return True


def test_order_statistic_is_a_sharp_adversarial_threshold():
    margins=(1000,2000,3000,9000,9000)
    for e,t in ((0,0),(0,1),(1,0),(2,0),(2,1)):
        b=risk_channels_for_budget(5,e,t)
        bound=sorted(margins)[b]
        assert brute_force_majority_guarantee(margins,Fraction(bound)-Fraction(1,2),e,t)
        assert not brute_force_majority_guarantee(margins,Fraction(bound),e,t)


def test_negative_preexisting_reporter_can_still_allow_conservative_majority():
    code,model=fake_code((-1000,2000,3000,9000,9000))
    out=certify_perturbation_radius(code,(("GA95",model),))
    assert out.geometry_risk_channels==2
    assert out.strict_pair_distance_radius_milliangstrom==3000
    assert decode_with_erasures(code,{
        e.endpoints:Fraction(
            20_000-int(x)
        ) for e,x in zip(code.reference.reporters,(-1000,2000,3000,9000,9000))
    })=="GA95"
    strict=certify_perturbation_radius(
        code,(("GA95",model),),erased_channels=2,additional_wrong_channels=1)
    assert strict.strict_pair_distance_radius_milliangstrom==0
    assert not strict.holds_for_all_submitted_conformers


@pytest.mark.parametrize("bad", ((5,2,1),(5,0,3),(3,3,0),(3,-1,0),
                                  (3,0,-1),(0,0,0),(True,0,0)))
def test_invalid_guarantees_rejected(bad):
    with pytest.raises(ValueError):
        risk_channels_for_budget(*bad)


def test_global_minimum_is_over_every_declared_conformer():
    code,a=fake_code((1000,2000,3000,9000,9000))
    _,b=fake_code((1000,2000,2500,9000,9000))
    cert=certify_perturbation_radius(code,(("GA95",a),("GA95",b)))
    assert cert.conformer_count==2
    assert cert.strict_pair_distance_radius_milliangstrom==2500
    assert check_perturbation_certificate(code,(("GA95",a),("GA95",b)),cert)
    with pytest.raises(ValueError):
        check_perturbation_certificate(code,(("GA95",a),),cert)
