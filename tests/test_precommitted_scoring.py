import importlib.util, pathlib, tempfile, pytest, csv

SCRIPT=pathlib.Path(__file__).resolve().parents[1]/'experiments'/'score_precommitted_fold_predictions.py'
spec=importlib.util.spec_from_file_location('score_precommitted_fold_predictions',SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def test_precise_multiclass_metrics_and_deterministic_order():
    p={'a':{'alpha':.8,'beta':.15,'unfolded':.05},'b':{'alpha':.1,'beta':.8,'unfolded':.1}}
    y={'a':{'alpha':1.0,'beta':0.,'unfolded':0.},'b':{'alpha':0.,'beta':.9,'unfolded':.1}}
    s=mod.score(p,y)
    assert s['variants']==2 and s['dominant_fold_accuracy']==1
    assert s['mean_cross_entropy_nats']>0
    assert .02<s['mean_multiclass_brier']<.06

def test_malformed_prediction_must_fail():
    with tempfile.TemporaryDirectory() as tmp:
        p=pathlib.Path(tmp)/'bad.csv'
        p.write_text('id,alpha,beta,unfolded,model_name,model_git_sha,trained_until\nX,0.9,0.9,0.1,m,hash,2026-01-01\n')
        with pytest.raises(ValueError):mod.load(p,True)
        p.write_text('id,alpha,beta,unfolded,model_name,model_git_sha,trained_until\nX,0.8,0.1,0.1,m,hash,2026-01-01\nX,0.8,0.1,0.1,m,hash,2026-01-01\n')
        with pytest.raises(ValueError):mod.load(p,True)
        p.write_text('id,alpha,beta,unfolded,model_name,model_git_sha,trained_until\nX,0.8,0.1,0.1,m,,2026-01-01\n')
        with pytest.raises(ValueError):mod.load(p,True)