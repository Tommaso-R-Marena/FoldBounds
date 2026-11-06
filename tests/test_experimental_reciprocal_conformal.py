"""Unit/property/adversarial checks; no external network or S669 CSV required."""
from dataclasses import replace
from decimal import Decimal as D
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.experimental_reciprocal_conformal import (
    MODELS, AssuranceError, Measurement, antisymmetric_projection,
    reciprocity_pythagoras, git_blob_sha1, load_s669, protein_folds,
    choose_predictor, conformal_cutoff, evaluate, summarize,
)


def observation(i, group, truth, fwd, rev, *, other=None):
    w=group
    m=('H' if w[0]!='H' else 'A')+w[1:]
    pred={n+'_dir':D(str(fwd if other is None else other)) for n in MODELS}
    pred.update({n+'_inv':D(str(rev if other is None else -D(str(other)))) for n in MODELS})
    pred['GeoDDG-3D_dir']=D(str(fwd))
    pred['GeoDDG-3D_inv']=D(str(rev))
    return Measurement(str(i),w,m,D('7'),D('25'),D(str(truth)),pred)


def sample():
    base='ACDEFGHIKLMNPQRSTVWY'
    return [observation(f'{i}_{j}',base[i%20:]+'A'*i,
                        D(j-2)/2+D(i%3)/3,D(j-2)/2,D(2-j)/2)
            for i in range(25) for j in range(1+(i%5))]


@pytest.mark.parametrize("forward,reverse,true",product(
    [D('-3.125'),D(0),D('1.2')],
    [D('-2.4'),D(0),D('3.2')],
    [D('-10.1'),D(0),D('0.9')]))
def test_orthogonal_thermodynamic_projection(forward,reverse,true):
    projected,negative=antisymmetric_projection(forward,reverse)
    assert projected==-negative
    assert reciprocity_pythagoras(forward,reverse,true)
    assert 2*(projected-true)**2 <= (forward-true)**2+(reverse+true)**2


def test_forward_mae_can_worsen_despite_paired_l2_guarantee():
    projected,_=antisymmetric_projection(D(1),D(-10))
    assert abs(projected-D(1)) > abs(D(1)-D(1))


def test_git_blob_provenance_reference():
    assert git_blob_sha1(b'test')=='30d74d258442c7c65512eafab474568dd706c430'


def test_synthetic_source_can_be_loaded_but_never_pretend_it_is_s669(tmp_path):
    fields=['name','wt_seq','mut_seq','pH','temp','ddG']
    fields += [n+'_'+mode for n in MODELS for mode in ('dir','inv')]
    header=','.join(fields)+'\n'
    values=['variant1','AAAA','HAAA','7','25','-1']+['-0.7','0.7']*len(MODELS)
    path=tmp_path/'sample.csv'
    path.write_text(header+','.join(values)+'\n',encoding='utf-8')
    assert len(load_s669(path,verify_provenance=False))==1
    with pytest.raises(AssuranceError,match='pinned'):
        load_s669(path)
    for bad in ('AAAA','HAAH'):
        values[2]=bad
        path.write_text(header+','.join(values)+'\n')
        with pytest.raises(AssuranceError,match='single-site'):
            load_s669(path,verify_provenance=False)
    values[2]='HAAA';values[5]='NaN'
    path.write_text(header+','.join(values)+'\n')
    with pytest.raises(AssuranceError):
        load_s669(path,verify_provenance=False)


@pytest.mark.parametrize('seed',[1,20261105,20261106,2**32-1])
def test_train_calibration_test_protein_sets_disjoint(seed):
    mapping,groups=protein_folds(sample(),seed)
    assert len(groups)==25
    for fold in range(5):
        test={g for g in groups if mapping[g]==fold}
        calibration={g for g in groups if mapping[g]==(fold+1)%5}
        train=set(groups)-test-calibration
        assert not test&calibration and not test&train and not calibration&train
        assert len(test|calibration|train)==25
        assert len(test)==len(calibration)==5


@pytest.mark.parametrize('bad_seed',[0,-1,2**32,2.3,True])
def test_seed_rejections(bad_seed):
    with pytest.raises(AssuranceError):
        protein_folds(sample(),bad_seed)


def test_exact_group_conformal_quantile_and_infinite_abstention():
    records=[observation(i,f'P{i}AAA',D(i),D(0),D(0)) for i in range(1,11)]
    groups={r.wt:[r] for r in records}
    keys=list(groups)
    assert conformal_cutoff(keys[:9],groups,'GeoDDG-3D_dir',Fraction(1,10))==9
    assert conformal_cutoff(keys[:9],groups,'GeoDDG-3D_dir',Fraction(1,5))==8
    assert conformal_cutoff(keys[:9],groups,'GeoDDG-3D_dir',Fraction(1,20)) is None
    for alpha in (0.1,Fraction(0),Fraction(1),Fraction(-1,10)):
        with pytest.raises(AssuranceError):
            conformal_cutoff(keys,groups,'GeoDDG-3D_dir',alpha)
    with pytest.raises(AssuranceError):
        conformal_cutoff(keys[:2]*2,groups,'GeoDDG-3D_dir')


def test_exchangeability_rank_argument_by_enumeration():
    data=[observation(i,f'P{i}AAA',i,0,0) for i in range(1,11)]
    groups={r.wt:[r] for r in data}
    successes=0
    for candidate in data:
        keys=[r.wt for r in data if r.wt!=candidate.wt]
        bound=conformal_cutoff(keys,groups,'GeoDDG-3D_dir',Fraction(1,10))
        successes+=int(abs(candidate.ddg-candidate.predict('GeoDDG-3D_dir'))<=bound)
    assert successes==9


@pytest.mark.parametrize('strategy',['fixed','selected_phys','selected_any'])
def test_all_folds_and_exact_group_coverage(strategy):
    original=sample()
    f=evaluate(original,alpha=Fraction(1,10),strategy=strategy)
    assert len(f)==5
    s=summarize(f)
    assert s['proteins_total']==25 and s['mutations_total']==len(original)
    assert sum(v['cal_proteins'] for v in f)==25
    assert all(v['test_proteins']==5 for v in f)
    if strategy=='fixed':
        assert all(v['predictor']=='GeoDDG-3D_dir' for v in f)
    if strategy=='selected_phys':
        assert all(v['predictor'].endswith('_sym') for v in f)


def test_input_permutation_does_not_change_report():
    for strategy in ('fixed','selected_phys','selected_any'):
        a=summarize(evaluate(sample(),strategy=strategy))
        b=summarize(evaluate(list(reversed(sample())),strategy=strategy))
        assert a==b


def test_test_target_tampering_cannot_change_training_selection_or_calibration():
    observations=sample()
    mapping,groups=protein_folds(observations)
    train=[g for g in groups if mapping[g] in (2,3,4)]
    cal=[g for g in groups if mapping[g]==1]
    heldout=next(g for g in groups if mapping[g]==0)
    choice=choose_predictor(train,groups)
    bound=conformal_cutoff(cal,groups,choice)
    revised=[replace(x,ddg=D(1000000)) if x.wt==heldout else x for x in observations]
    _,changed=protein_folds(revised)
    assert choose_predictor(train,changed)==choice
    assert conformal_cutoff(cal,changed,choice)==bound


def test_naive_mutation_calibration_does_not_claim_group_coverage():
    for alpha in (Fraction(1,10),Fraction(1,5)):
        group=evaluate(sample(),alpha=alpha)
        row=evaluate(sample(),alpha=alpha,calibration_unit='mutation')
        for g,r in zip(group,row):
            if g['cutoff'] is not None and r['cutoff'] is not None:
                assert D(r['cutoff'])<=D(g['cutoff'])
            assert r['proteins_covered']<=g['proteins_covered']


def test_untrusted_modes_fail_closed():
    with pytest.raises(AssuranceError):
        evaluate(sample(),strategy='test_optimized')
    with pytest.raises(AssuranceError):
        evaluate(sample(),calibration_unit='wishful')
    with pytest.raises(AssuranceError):
        conformal_cutoff([],{},'GeoDDG-3D_dir')
    with pytest.raises(AssuranceError):
        choose_predictor([],{})
