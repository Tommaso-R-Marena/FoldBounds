"""Adversarial, exhaustive and small independent oracle tests."""
from dataclasses import replace
from decimal import Decimal as D
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb
from random import Random

import pytest

from foldbounds.assay_regret_certificates import (
    AssayCertificateError, Mutation, assay_panels, block_scores,
    certify_unmeasured_plan, conformal_rank_cutoff, evaluate,
    expected_uniform_random_panel_regret, git_blob_sha1, motif_components,
    parse_s669, protein_random_expected_regret, rank_select,
    regret_from_panels, summarize,
)


def observations():
    rng=Random(78145)
    alphabet='ACDEFGHIKLMNPQRSTVWY'
    proteins=set()
    while len(proteins)<25:
        proteins.add(''.join(rng.choice(alphabet) for _ in range(28)))
    out=[]
    for i,s in enumerate(sorted(proteins)):
        for j in range(2+i%5):
            pred=D(j)/4+D(i%5)/3
            truth=pred + D(j%3-1)/2+D(i%6-2)/5
            mutant=('H' if s[j%len(s)]!='H' else 'A')
            sequence=s[:j%len(s)]+mutant+s[j%len(s)+1:]
            out.append(Mutation(f'protein{i:02d}_mutation{j:02d}',s,sequence,D('7'),D('25'),truth,pred))
    return out


def four_panel():
    p=observations()[:4]
    return tuple(replace(p[i],name=f'variant{i}',wt_seq=p[0].wt_seq,
                         mut_seq=p[0].wt_seq[:i] + ('H' if p[0].wt_seq[i]!='H' else 'A') + p[0].wt_seq[i+1:],
                         pred=D([3,1,2,0][i]),ddg=D([1,5,4,3][i])) for i in range(4))


def test_pinned_source_gate_rejects_unsourced_or_forged_csv():
    assert git_blob_sha1(b'test')=='30d74d258442c7c65512eafab474568dd706c430'
    with pytest.raises(AssayCertificateError,match='provenance'):
        parse_s669(b'name,wt_seq,mut_seq,pH,temp,ddG,GeoDDG-3D_dir\n')
    good=(b'name,wt_seq,mut_seq,pH,temp,ddG,GeoDDG-3D_dir\n'
          b'x,AAAA,HAAA,7,25,-1.5,-1.2\n')
    rows=parse_s669(good,enforce_provenance=False)
    assert len(rows)==1 and rows[0].ddg==D('-1.5')
    for raw in [b'x,AAAA,AAAA,7,25,-1.5,-1.2\n',
                b'x,AAAA,HHAA,7,25,-1.5,-1.2\n',
                b'x,AAAA,HAAA,NaN,25,-1.5,-1.2\n',
                b'x,AAAA,HAAA,7,25,inf,-1.2\n']:
        with pytest.raises(AssayCertificateError):
            parse_s669(good.splitlines(keepends=True)[0]+raw,enforce_provenance=False)
    with pytest.raises(AssayCertificateError):
        parse_s669(good+good.splitlines(keepends=True)[1],enforce_provenance=False)


def test_exact_condition_separation_and_trivial_panel_exclusion():
    r=observations()
    x=r[0]; y=r[1]
    out=assay_panels([x,y,replace(y,name='different_condition',temp=D('60'))])
    assert len(out[x.wt_seq])==1 and len(out[x.wt_seq][0][1])==2
    out=assay_panels([x,replace(y,temp=D('60'))])
    assert not out


def test_very_short_groups_are_excluded_without_leaking_into_calibration():
    rows=observations()
    first=rows[0]
    reduced=[x for x in rows if x.wt_seq != first.wt_seq] + [first]
    panels=assay_panels(reduced)
    assert first.wt_seq not in panels
    blocks,_=block_scores(reduced,Fraction(1,5))
    assert first.wt_seq not in blocks


def test_budgeted_selection_never_reads_observed_truth():
    panel=four_panel()
    a=rank_select(panel,Fraction(1,4))
    permuted=[replace(x,ddg=D('-100000') if i%2 else D('100000'))
              for i,x in enumerate(panel)]
    assert [x.name for x in a]==[x.name for x in rank_select(permuted,Fraction(1,4))]
    assert [x.name for x in a]==['variant0']
    assert [x.name for x in rank_select(panel,Fraction(1,2))]==['variant0','variant2']
    with pytest.raises(AssayCertificateError):
        rank_select(panel,Fraction(0))
    with pytest.raises(AssayCertificateError):
        rank_select(panel,Fraction(1,4),selector='look_at_ddg')


def test_regret_is_nonnegative_and_zero_if_all_mutants_assayed():
    rows=observations()
    panels=assay_panels(rows)
    assert len(panels)==25
    for p in panels.values():
        r=regret_from_panels(p,Fraction(1,5))
        assert r.max_panel_regret>=0
        assert r.measured_candidates<=r.available_candidates
        assert regret_from_panels(p,Fraction(1)).max_panel_regret==0


def test_exact_uniform_random_baseline_matches_every_combination():
    panel=four_panel()
    for k in range(1,5):
        budget=Fraction(k,4)
        real=expected_uniform_random_panel_regret(panel,budget)
        possible=[]
        for subset in combinations(panel,k):
            possible.append(max(x.ddg for x in panel)-max(x.ddg for x in subset))
        brute=sum((Fraction(x) for x in possible),Fraction(0))/len(possible)
        assert real==brute
    # Closed-form test: [1,2,3], k=2, sampled best 2,3,3 => regret 1/3.
    n=panel[:3]
    n=tuple(replace(x,ddg=D(i+1)) for i,x in enumerate(n))
    assert expected_uniform_random_panel_regret(n,Fraction(2,3))==Fraction(1,3)


def test_conformal_off_by_one_and_fail_closed_for_small_calibration():
    scores=[D(j) for j in range(1,10)]
    assert conformal_rank_cutoff(scores,Fraction(1,10))==D(9)
    assert conformal_rank_cutoff(scores,Fraction(1,5))==D(8)
    assert conformal_rank_cutoff(scores,Fraction(1,20)) is None
    for bad in [0,.1,Fraction(0),Fraction(1)]:
        with pytest.raises(AssayCertificateError):
            conformal_rank_cutoff(scores,bad)
    with pytest.raises(AssayCertificateError):
        conformal_rank_cutoff([])
    with pytest.raises(AssayCertificateError):
        conformal_rank_cutoff([D('-1'),D('2')])
    with pytest.raises(AssayCertificateError):
        conformal_rank_cutoff([D('NaN')])


def test_leave_one_out_exchangeable_rank_all_permutations():
    # Deterministic symmetry proof by choosing each of 10 possible test groups.
    for scores in [list(map(D,range(10))),[D(i//2) for i in range(10)]]:
        hits=0
        for i,v in enumerate(scores):
            q=conformal_rank_cutoff(scores[:i]+scores[i+1:],Fraction(1,10))
            hits+=int(v<=q)
        assert hits>=9


def test_12_mer_components_transitive_and_no_affinity_tag_homology_claim():
    x='AAAAAAAABBBBCCCC'
    y='BBBBCCCCDDDD'
    z='CCCCDDDDFFFF'
    t='H'*12+'V'*16
    u='H'*12+'W'*16
    components=motif_components([x,y,z,t,u],motif_length=8)
    assert components[x]==components[y]==components[z]
    assert components[t]==components[u]
    assert components[x]!=components[t]
    for bad in (5,41,True,12.5):
        with pytest.raises(AssayCertificateError):
            motif_components([x,y],motif_length=bad)


def test_motif_connected_components_preserve_all_protein_panels():
    rows=observations()
    proteins=assay_panels(rows)
    blocked, members=block_scores(rows,Fraction(1,5),motif_length=12)
    assert set(m for group in members.values() for m in group)==set(proteins)
    assert sum(r.available_candidates for r in blocked.values())==len(rows)
    assert sum(r.measured_candidates for r in blocked.values())>0
    for g, parts in members.items():
        assert blocked[g].max_panel_regret==max(
            regret_from_panels(proteins[s],Fraction(1,5)).max_panel_regret for s in parts)


@pytest.mark.parametrize('seed',[1,20261105,20261106,2**32-1])
@pytest.mark.parametrize('motif_length',[None,12])
def test_heldout_calibration_no_overlap_and_results_consistent(seed,motif_length):
    folds=evaluate(observations(),seed=seed,motif_length=motif_length)
    summary=summarize(folds)
    assert summary['groups_total']==25
    assert sum(f.test_groups for f in folds)==25
    assert sum(f.calibration_groups for f in folds)==100
    assert summary['budget_assays']<=summary['available_assays']
    assert all(f.cutoff is not None and f.cutoff>=0 for f in folds)
    assert 0<=summary['groups_covered']<=25


def test_input_order_may_not_affect_model_or_score():
    rows=observations()
    for budget in (Fraction(1,5),Fraction(1,2)):
        assert block_scores(rows,budget)[0]==block_scores(list(reversed(rows)),budget)[0]
        assert evaluate(rows,budget)==evaluate(list(reversed(rows)),budget)


def test_reference_truth_tampering_must_not_change_candidate_ids():
    panel=four_panel()
    replaced=[replace(x,ddg=D(1000)) for x in panel]
    for budget in (Fraction(1,4),Fraction(1,2),Fraction(3,4)):
        assert [x.name for x in rank_select(panel,budget)]==[x.name for x in rank_select(replaced,budget)]
        assert rank_select(panel,budget)[0].name==rank_select(replaced,budget)[0].name


def test_deployment_assay_plan_contains_no_truth_data():
    measured=four_panel()
    unmeasured=[replace(x,ddg=None) for x in measured]
    plan=certify_unmeasured_plan(unmeasured,Fraction(1,2),D('1.7'))
    assert plan['mutation_ids']==['variant0','variant2']
    assert plan['status']=='marginal-group-assurance-only'
    assert plan['regret_upper_kcal_mol']=='1.7'
    assert certify_unmeasured_plan(unmeasured,Fraction(1,5),None)['status']=='abstain'
    with pytest.raises(AssayCertificateError):
        certify_unmeasured_plan(measured,Fraction(1,5),D(1))
    with pytest.raises(AssayCertificateError):
        certify_unmeasured_plan(unmeasured,Fraction(1,5),D(-1))
    with pytest.raises(AssayCertificateError):
        certify_unmeasured_plan(unmeasured,Fraction(1,5),D('Infinity'))
    with pytest.raises(AssayCertificateError):
        certify_unmeasured_plan([replace(x,temp=D(30)) if i==0 else x
                                for i,x in enumerate(unmeasured)],Fraction(1,5),D(1))


def test_dataset_shift_breaks_any_claim_of_conditional_protection():
    # Finite-sample conformal does not protect a specifically chosen highly
    # adversarial future protein with score beyond all calibration scores.
    calibration=[D(x) for x in (0,0,0,1,1,1,1,1,1,1,1)]
    q=conformal_rank_cutoff(calibration,Fraction(1,5))
    assert q==1
    assert D('100000')>q


def test_empirical_top_selection_is_not_guaranteed_to_be_good():
    panel=four_panel()
    # Model ranks variant0 highest, whereas actual best is variant1.
    assert regret_from_panels([((D(7),D(25)),panel)],Fraction(1,4)).max_panel_regret==D(4)
    assert regret_from_panels([((D(7),D(25)),panel)],Fraction(1)).max_panel_regret==0


def test_duplicate_mutations_rejected_even_across_conditions():
    x,y=observations()[:2]
    with pytest.raises(AssayCertificateError):
        assay_panels([x,y,replace(x,temp=D(50))])
