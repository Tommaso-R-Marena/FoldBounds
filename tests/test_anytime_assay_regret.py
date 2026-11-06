"""Diverse synthetic/hard adversarial tests; no network or copied S669 rows."""
from dataclasses import replace
from decimal import Decimal as D
from fractions import Fraction as F
from itertools import product
import random

import pytest

from foldbounds.assay_regret_certificates import Mutation, AssayCertificateError, rank_select
from foldbounds.anytime_assay_regret import (
    DEFAULT_BUDGETS, AnytimeCertificate, check_grid, evidence_for_panel,
    assay_group_evidence, group_score, propose_certificate,
    signed_conformal_cutoff, robust_noise_radius, bound_from_assayed,
    retrospectively_audit,
)
from foldbounds.anytime_assay_regret_verifier import verify_anytime_certificate


def example_protein(i, n=8, conditions=1):
    # Unique wildtypes of fixed length across all i, no accidentally shared
    # motifs. Distinct predictions/outcomes exercise bad and good rankings.
    alphabet='ACDEFGHIKLMNPQRSTVWY'
    seq=''.join(random.Random(70000+i).choices(alphabet,k=32))
    out=[]
    for condition in range(conditions):
        for j in range(n):
            wt=seq
            mut=wt[:j%len(wt)]+('H' if wt[j%len(wt)]!='H' else 'A')+wt[j%len(wt)+1:]
            y=D((j*17+i*11+condition*5)%29-15)/D(7)
            p=D((j*13+i*7+condition*3)%31-16)/D(8)
            out.append(Mutation(f'P{i}_c{condition}_m{j}',wt,mut,D('7'),
                                D(20+condition*10),y,p))
    return out


def sample(proteins=23):
    return [x for i in range(proteins) for x in example_protein(i,n=3+i%9,conditions=1+i%2)]


def future_panel(panel):
    return tuple(replace(x,ddg=None) for x in panel)


def selected_observations(panel,budget,epsilon=F(0),sign=1):
    ans={}
    for x in rank_select(panel,budget):
        # Fake a measurement reported within a bounded deterministic epsilon
        ans[x.name]=x.ddg+D(epsilon.numerator*sign)/D(epsilon.denominator)
    return ans


def test_exact_construction_and_independent_audit():
    rows=sample()
    cal=sorted(assay_group_evidence(rows,DEFAULT_BUDGETS))[:16]
    c=propose_certificate(rows,cal)
    result=verify_anytime_certificate(c,rows)
    assert result['verified'] and not result['abstains']
    assert c.cutoff==sorted(c.exact_group_scores)[15]
    assert len(c.budgets)==4


@pytest.mark.parametrize('alpha',[F(1,10), F(1,5), F(1,2)])
def test_exact_signed_finite_sample_rank_and_abstention(alpha):
    a=[F(-10),F(-5),F(0),F(2),F(3),F(5),F(9),F(10),F(11)]
    k=((len(a)+1)*(alpha.denominator-alpha.numerator)+alpha.denominator-1)//alpha.denominator
    assert signed_conformal_cutoff(a,alpha)==(a[k-1] if k<=len(a) else None)
    if alpha==F(1,10):
        assert signed_conformal_cutoff(a[:4],alpha) is None
    assert signed_conformal_cutoff([F(-3)]*19,F(1,10))==F(-3)


@pytest.mark.parametrize('budget',[F(1,20),F(1,10),F(1,5),F(1,2),F(1)])
def test_model_selection_is_label_blind_and_declared_regret_is_exact(budget):
    p=tuple(example_protein(5,9))
    evidence=evidence_for_panel(p,budget)
    k=max(1,(len(p)*budget.numerator+budget.denominator-1)//budget.denominator)
    assert evidence.assayed==k and evidence.total==len(p)
    assert evidence.selected==tuple(x.name for x in rank_select(p,budget))
    all_best=max(x.ddg for x in p)
    selected_best=max(x.ddg for x in rank_select(p,budget))
    assert evidence.actual_regret==F(all_best)-F(selected_best)
    changed=tuple(replace(x,ddg=D(-10000+i)) for i,x in enumerate(p))
    assert evidence_for_panel(changed,budget).selected==evidence.selected
    if budget==1:
        assert evidence.observed_margin==evidence.actual_regret==0


def test_complete_coverage_of_budgets_in_finite_rank_proof():
    # With n=20 distinct *complete* groups and alpha=1/10,
    # leaving any one out makes at least 18/20 exchangeably possible
    # held-out groups covered simultaneously at ALL four budgets.
    rows=sample(20)
    evidence=assay_group_evidence(rows,DEFAULT_BUDGETS)
    keys=tuple(sorted(evidence))
    passes=0
    for heldout in keys:
        c=propose_certificate(rows,(g for g in keys if g!=heldout))
        passes+=int(retrospectively_audit(c,evidence[heldout])['simultaneous_coverage'])
    assert passes>=18


def test_selected_future_measurements_and_anytime_stopping():
    rows=sample(24)
    full=assay_group_evidence(rows,DEFAULT_BUDGETS)
    holdout=sorted(full)[-1]
    c=propose_certificate(rows,(g for g in sorted(full) if g!=holdout))
    p=tuple(x for x in rows if x.wt_seq==holdout)
    # Multi-condition future experiment: evaluate separate panels, never mix.
    by_cond={}
    for x in p:
        by_cond.setdefault((x.ph,x.temp),[]).append(x)
    for panel in by_cond.values():
        truth=tuple(panel)
        unobserved=future_panel(truth)
        for b in DEFAULT_BUDGETS:
            obs=selected_observations(truth,b)
            bound=bound_from_assayed(c,unobserved,b,obs)
            assert bound is None or bound>=0
            assert len(obs)==max(1,(len(panel)*b.numerator+b.denominator-1)//b.denominator)
            # A *retrospective* failing group may violate the nominal
            # certificate, which is not a conditional guarantee.
            if retrospectively_audit(c,full[holdout])['simultaneous_coverage']:
                assert evidence_for_panel(truth,b).actual_regret<=bound
    if retrospectively_audit(c,full[holdout])['simultaneous_coverage']:
        # Any selection of a budget based on observed assay results inherits
        # the same simultaneous event; choosing the stop does not alter q.
        for panel in by_cond.values():
            available=tuple(panel)
            stop=next((b for b in DEFAULT_BUDGETS
                       if bound_from_assayed(c,future_panel(available),b,
                           selected_observations(available,b))<=F(5)),DEFAULT_BUDGETS[-1])
            assert evidence_for_panel(available,stop).actual_regret<=bound_from_assayed(
                c,future_panel(available),stop,selected_observations(available,stop))


def test_independent_verifier_rejects_tampered_rank_and_scores():
    rows=sample(24)
    group_ids=sorted(assay_group_evidence(rows,DEFAULT_BUDGETS))[:18]
    cert=propose_certificate(rows,group_ids)
    assert verify_anytime_certificate(cert,rows)['verified']
    corrupt=[
        replace(cert,cutoff=F(-999)),
        replace(cert,exact_group_scores=cert.exact_group_scores[:-1]+(F(999999),)),
        replace(cert,calibration_group_ids=cert.calibration_group_ids[1:]),
        replace(cert,calibration_group_ids=cert.calibration_group_ids+(cert.calibration_group_ids[0],)),
        replace(cert,alpha=F(0)),
        replace(cert,weights=(F(0),)+cert.weights[1:]),
        replace(cert,budgets=tuple(reversed(cert.budgets))),
        replace(cert,lambda_margin=F(-1)),
    ]
    for c in corrupt:
        with pytest.raises(AssayCertificateError):
            verify_anytime_certificate(c,rows)


def test_deployment_does_not_look_at_unmeasured_truth():
    rows=sample(24)
    group_ids=sorted(assay_group_evidence(rows,DEFAULT_BUDGETS))[:19]
    cert=propose_certificate(rows,group_ids)
    panel=tuple(example_protein(30,10))
    b=F(1,5)
    measured=selected_observations(panel,b)
    trial=future_panel(panel)
    valid=bound_from_assayed(cert,trial,b,measured)
    assert valid is not None and valid>=0
    with pytest.raises(AssayCertificateError):
        bound_from_assayed(cert,panel,b,measured)
    with pytest.raises(AssayCertificateError):
        bound_from_assayed(cert,trial,b,{})
    with pytest.raises(AssayCertificateError):
        bound_from_assayed(cert,trial,b,{**measured,'forged_unselected':D(123)})
    with pytest.raises(AssayCertificateError):
        bound_from_assayed(cert,trial,F(3,7),measured)
    with pytest.raises(AssayCertificateError):
        bound_from_assayed(cert,trial,b,{k:D('NaN') for k in measured})


def test_assay_conditions_physical_guard():
    p=tuple(example_protein(2,5,conditions=2))
    with pytest.raises(AssayCertificateError,match='condition'):
        evidence_for_panel(p,F(1,5))
    with pytest.raises(AssayCertificateError):
        evidence_for_panel(p[:0],F(1,5))
    bad=tuple(example_protein(3,4))
    with pytest.raises(AssayCertificateError):
        evidence_for_panel(bad+(bad[0],),F(1,5))


def test_noise_perturbation_constant_exact_and_sharp_bound():
    rows=sample(24)
    cal=sorted(assay_group_evidence(rows,DEFAULT_BUDGETS))[:20]
    c=propose_certificate(rows,cal,lambda_margin=F(1))
    eps=F(1,4)
    assert robust_noise_radius(c,eps)==F(3,4)
    # With unit weights the total conservatism is (2+2lambda)*epsilon.
    b=F(1,5)
    panel=tuple(example_protein(90,8))
    observed=selected_observations(panel,b,epsilon=eps,sign=-1)
    plain=bound_from_assayed(c,future_panel(panel),b,
                            selected_observations(panel,b))
    noisy=bound_from_assayed(c,future_panel(panel),b,observed,
                            assay_error_epsilon=eps)
    assert noisy is not None and noisy>=0
    assert noisy>=plain or plain is None  # selected underreported => larger margin
    # Under an arbitrary plus/minus eps measurement perturbation, the true
    # regret bound is sound whenever the original simultaneous score event holds.
    e=evidence_for_panel(panel,b)
    true_margin=e.observed_margin
    obs_margin=evidence_for_panel(future_panel(panel),b,measured=observed,
                                  strict_deployment=True).observed_margin
    assert abs(obs_margin-true_margin)<=eps


@pytest.mark.parametrize('lambda_value',[F(0),F(1,4),F(1),F(3)])
@pytest.mark.parametrize('weight_spec',[(F(1),)*4,(F(2),F(3,2),F(1),F(1,2))])
def test_exact_producer_verifier_alternate_weights_and_margins(lambda_value,weight_spec):
    rows=sample(17)
    keys=sorted(assay_group_evidence(rows,DEFAULT_BUDGETS))[:12]
    cert=propose_certificate(rows,keys,lambda_margin=lambda_value,weights=weight_spec)
    assert verify_anytime_certificate(cert,rows)['verified']
    assert robust_noise_radius(cert,F(1,2))==max((F(2)+lambda_value)*F(1,2)/w for w in weight_spec)


def test_motif_connected_components_independently_audited():
    rows=sample(25)
    # Some generated wildtypes have long shared motifs; independent auditor
    # must reproduce transitive grouping without calling proposal grouping.
    ev=assay_group_evidence(rows,DEFAULT_BUDGETS,motif_length=12)
    c=propose_certificate(rows,tuple(ev)[:-1],motif_length=12)
    assert verify_anytime_certificate(c,rows)['verified']


@pytest.mark.parametrize('b,w,lam',[
    ((F(1,5),F(1,5)),(F(1),F(1)),F(1)),
    ((F(1,5),),(F(0),),F(1)),
    ((F(0),),(F(1),),F(1)),
    ((F(1,5),),(F(1),),F(-1)),
    ((F(1,5),),(F(1,2),F(2)),F(1)),
])
def test_bad_budget_grid_fail_closed(b,w,lam):
    with pytest.raises(AssayCertificateError):
        check_grid(b,w,lam)


def test_calibration_too_small_must_abstain():
    rows=sample(20)
    groups=sorted(assay_group_evidence(rows,DEFAULT_BUDGETS))[:3]
    c=propose_certificate(rows,groups,alpha=F(1,100))
    assert c.cutoff is None
    assert retrospectively_audit(c,assay_group_evidence(rows,DEFAULT_BUDGETS)[groups[0]])['calibration_abstained']
    p=future_panel(tuple(example_protein(88,8)))
    assert bound_from_assayed(c,p,F(1,5),{x.name:D(1) for x in rank_select(p,F(1,5))}) is None


def test_no_spurious_protein_guarantee_for_nonexchangeable_groups():
    # A distribution shift can have arbitrarily larger regret than ANY
    # interval computed on the calibration source.
    rows=sample(24)
    c=propose_certificate(rows,sorted(assay_group_evidence(rows,DEFAULT_BUDGETS))[:20])
    assert verify_anytime_certificate(c,rows)['verified']
    # The checker cannot and must not claim future exchangeability.
    assert 'exchangeab' in verify_anytime_certificate(c,rows)['checker']

@pytest.mark.parametrize('epsilon',[F(0),F(1,10),F(1,4),F(1,2)])
@pytest.mark.parametrize('lam',[F(0),F(1,4),F(1)])
def test_adversarial_measurement_perturbation_cannot_undercut_robust_bound(epsilon,lam):
    # Perturb ALL calibrating measurements adversarially. A single known
    # uniform worst-case error epsilon propagates to the conformal rank and
    # to future measured anchors. This is a deterministic sensitivity proof.
    truth=sample(22)
    labels=lambda item: 1 if sum(map(ord,item.name))%2 else -1
    measured=[replace(x,ddg=x.ddg+D(labels(x)*epsilon.numerator)/D(epsilon.denominator))
              for x in truth]
    keys=sorted(assay_group_evidence(truth,DEFAULT_BUDGETS))[:19]
    true_cert=propose_certificate(truth,keys,lambda_margin=lam)
    measured_cert=propose_certificate(measured,keys,lambda_margin=lam)
    assert true_cert.cutoff is not None and measured_cert.cutoff is not None
    eta=robust_noise_radius(measured_cert,epsilon)
    assert abs(measured_cert.cutoff-true_cert.cutoff)<=eta
    heldout=sorted(set(x.wt_seq for x in truth)-set(keys))[0]
    panels={}
    for x in truth:
        if x.wt_seq==heldout:
            panels.setdefault((x.ph,x.temp),[]).append(x)
    # When the true score has the nominal simultaneous coverage event,
    # every exposed prospective horizon must obey the robust bound even
    # after adversarial ±epsilon perturbations in all selected assays.
    true_group=assay_group_evidence(truth,DEFAULT_BUDGETS)[heldout]
    if retrospectively_audit(true_cert,true_group)['simultaneous_coverage']:
        for panel in panels.values():
            for b in DEFAULT_BUDGETS:
                selection=rank_select(panel,b)
                obs={x.name: x.ddg+D(labels(x)*epsilon.numerator)/D(epsilon.denominator)
                     for x in selection}
                upper=bound_from_assayed(measured_cert,future_panel(panel),b,obs,
                                        assay_error_epsilon=epsilon)
                assert evidence_for_panel(tuple(panel),b).actual_regret<=upper


def test_verified_deployment_wrapper_detects_tampering_and_exact_sequence_leakage():
    from foldbounds.anytime_assay_regret_verifier import verified_deployment_bound
    rows=sample(25)
    group_ids=sorted(assay_group_evidence(rows,DEFAULT_BUDGETS))[:19]
    c=propose_certificate(rows,group_ids)
    panel=tuple(example_protein(99,10))
    future=future_panel(panel)
    measurements=selected_observations(panel,F(1,5))
    assert verified_deployment_bound(c,rows,future,F(1,5),measurements) is not None
    with pytest.raises(AssayCertificateError):
        verified_deployment_bound(replace(c,cutoff=F(-999)),rows,future,F(1,5),measurements)
    existing=tuple(x for x in rows if x.wt_seq==group_ids[0])
    with pytest.raises(AssayCertificateError,match='same wildtype'):
        verified_deployment_bound(c,rows,future_panel(existing),F(1,5),
                                  selected_observations(existing,F(1,5)))


def test_verified_deployment_wrapper_rejects_motif_but_not_unrelated_fold():
    from foldbounds.anytime_assay_regret_verifier import verified_deployment_bound
    rows=sample(30)
    ev=assay_group_evidence(rows,DEFAULT_BUDGETS,motif_length=12)
    ids=tuple(sorted(ev)[:20])
    cert=propose_certificate(rows,ids,motif_length=12)
    assert verify_anytime_certificate(cert,rows)['verified']
    test_p=tuple(example_protein(299,8))
    assert verified_deployment_bound(cert,rows,future_panel(test_p),F(1,5),
                                      selected_observations(test_p,F(1,5))) is not None
    known_seq=ids[0]
    related=tuple(replace(x,wt_seq=x.wt_seq[:12]+known_seq[:12],
                          name='related_'+x.name) for x in test_p)
    with pytest.raises(AssayCertificateError,match='motif'):
        verified_deployment_bound(cert,rows,future_panel(related),F(1,5),
                                  selected_observations(related,F(1,5)))

@pytest.mark.parametrize('n',[30,55])
def test_end_to_end_benchmark_handles_bonferroni_impossibility_without_forgery(n):
    import runpy
    from pathlib import Path
    study=runpy.run_path(str(Path(__file__).resolve().parents[1]/'experiments/benchmark_anytime_assay_regret.py'))['study']
    result=study(sample(n),seed=20261105,motif_length=None,check_independent=True)
    assert result['protein_or_motif_groups']==n
    assert all(0<=r['joint_groups_covered']<=n for r in result['budget_results'].values())
    assert sum(result['adaptive_stops']['2']['counts'].values())==n
    if n==30:
        assert all(r['mean_bonf'] is None for r in result['budget_results'].values())
    if n==55:
        assert all(r['mean_bonf'] is not None for r in result['budget_results'].values())


def test_common_mode_offset_invariant_mean_anchor_is_a_no_go_for_topk():
    """Independent mathematical no-go: centered gap ALWAYS exactly zero.

    Let unselected_max_prediction U, mean selected prediction s_bar,
    mean selected measurement y_bar, best selected measurement ymax.
    Since U<=s_bar and y_bar<=ymax, U+(y_bar-s_bar)-ymax<=0.
    This precludes claims that such mean offset correction helps top-k.
    """
    for i in range(40):
        p=tuple(example_protein(100+i,3+i%10))
        for b in DEFAULT_BUDGETS:
            sel=rank_select(p,b)
            remaining=[x for x in p if x.name not in {r.name for r in sel}]
            if not remaining:
                continue
            u=max(F(x.pred) for x in remaining)
            sm=sum((F(x.pred) for x in sel),F(0))/len(sel)
            ym=sum((F(x.ddg) for x in sel),F(0))/len(sel)
            ymax=max(F(x.ddg) for x in sel)
            assert u<=sm and ym<=ymax
            assert max(F(0),u+ym-sm-ymax)==0
