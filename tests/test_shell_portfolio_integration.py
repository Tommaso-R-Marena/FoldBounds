"""Regression tests for independently checked variable-depth portfolio integration.

Requires FoldBounds main's enhanced lookahead/adaptive modules; these tests
are intentionally run by GitHub Python 3.11–3.13 CI on the stacked branch.
"""
from dataclasses import replace
from fractions import Fraction
import pytest

from foldbounds.lattice import ROOT, build_frontier, CertificateError, exact_spectrum
from foldbounds.shell_portfolio import (optimize_shell_portfolio,
                                        verify_shell_portfolio,
                                        PortfolioCertificate)
from foldbounds.lookahead import verify_lookahead
from foldbounds.adaptive import adaptive_refine,verify_refinement_run


@pytest.mark.parametrize('seq,q,threshold',[('HPHHPHHPHHHH',Fraction(3,2),2),
                                           ('H'*12,Fraction(2),2),
                                           ('HPHHPHHP',Fraction(2),1),
                                           ('P'*9,Fraction(3,2),2)])
def test_both_independently_certified_event_classes_tighten_coupled_probability(seq,q,threshold):
    cert=build_frontier(seq,q,threshold,2)
    kwargs=dict(lookahead=0,use_matching=True,use_transfer=True,
                use_turn_family=True,use_monotone=True)
    base=verify_lookahead(cert,**kwargs)
    plans=tuple(optimize_shell_portfolio(seq,ROOT,q,threshold,leaf_budget=25,
                   max_core_length=min(9,len(seq)),objective=kind)
                for kind in ('event','complement'))
    certs=tuple(plan.certificate for plan in plans)
    improved=verify_lookahead(cert,shell_portfolios=certs,**kwargs)
    assert improved.event_lower>=base.event_lower
    assert improved.complement_lower>=base.complement_lower
    assert improved.partition_lower>=base.partition_lower
    assert improved.event_upper==base.event_upper
    assert improved.complement_upper==base.complement_upper
    assert improved.partition_upper==base.partition_upper
    assert improved.probability_lower>=base.probability_lower
    assert improved.probability_upper<=base.probability_upper
    histogram=exact_spectrum(seq)
    event=sum((num*q**k for k,num in histogram.items() if k>=threshold),Fraction(0))
    other=sum((num*q**k for k,num in histogram.items() if k<threshold),Fraction(0))
    assert improved.event_lower<=event<=improved.event_upper
    assert improved.complement_lower<=other<=improved.complement_upper
    assert improved.partition_lower<=event+other<=improved.partition_upper
    assert improved.probability_lower<=event/(event+other)<=improved.probability_upper


def test_portfolio_parameter_substitution_rejected():
    seq='HPHHPHHP'; q=Fraction(2);t=1
    p=optimize_shell_portfolio(seq,ROOT,q,t,leaf_budget=10,max_core_length=7)
    cert=build_frontier(seq,q,t,2)
    for bad in [replace(p.certificate,q=Fraction(3,2)),
                replace(p.certificate,sequence='H'*len(seq)),
                replace(p.certificate,event_at_least=2),
                replace(p.certificate,prefix=ROOT+((2,0),))]:
        with pytest.raises(CertificateError):
            verify_lookahead(cert,lookahead=0,shell_portfolios=(bad,))
    with pytest.raises(ValueError):
        verify_lookahead(cert,lookahead=0,shell_portfolios=[p.certificate])


def test_adaptive_portfolio_recorded_replayed_and_detects_tampering():
    seq='HPHHPHHPHHHH';q=Fraction(3,2);t=2
    p=optimize_shell_portfolio(seq,ROOT,q,t,leaf_budget=25,max_core_length=9)
    initial=build_frontier(seq,q,t,2)
    run=adaptive_refine(initial,steps=2,lookahead=0,use_matching=True,
                        use_transfer=True,use_monotone=True,
                        shell_portfolios=(p.certificate,),policy='mass')
    assert run.shell_portfolios==(p.certificate,)
    assert verify_refinement_run(run)==run.bounds
    with pytest.raises(ValueError):
        verify_refinement_run(replace(run,shell_portfolios=()))
    with pytest.raises(CertificateError):
        verify_refinement_run(replace(run,shell_portfolios=(
            replace(p.certificate,sequence='H'*len(seq)),)))