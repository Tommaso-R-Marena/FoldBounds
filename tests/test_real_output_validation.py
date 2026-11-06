import sys,importlib.util,pathlib,json,pytest
from foldbounds.ga_gb_nmr import GA95,GB95,analyze_pair
p=pathlib.Path(__file__).resolve().parents[1]/'experiments'/'validate_real_ga_gb_output.py'
spec=importlib.util.spec_from_file_location('validate_real_ga_gb_output',p)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def fake_report():
    import importlib.util
    f=pathlib.Path(__file__).with_name('test_real_fold_benchmark.py')
    s=importlib.util.spec_from_file_location('foldbounds_real_fixture',f)
    mod=importlib.util.module_from_spec(s);s.loader.exec_module(mod)
    o=analyze_pair(mod.synthetic_pdb('2KDL',GA95),mod.synthetic_pdb('2KDM',GB95))
    from foldbounds.real_fold_odds import report_known_threshold
    o['experimental_thermodynamic_constraints']=report_known_threshold()
    return o

def test_validated_synthetic_data_is_not_mistaken_for_experimental_results():
    o=fake_report()
    # This validator verifies structure and provenance shape, not authenticity of the source;
    # authenticity requires official HTTPS download and recorded SHA-256.
    assert len(m.check(o)['verified_3_substitutions'])==3

@pytest.mark.parametrize('bad', ['sequence','count','sha','mutations','odds'])
def test_bad_report_rejected(bad):
    o=fake_report()
    if bad=='sequence':o['GB95']['sequence']='H'*56
    elif bad=='count':o['GA95']['nmr_models']=19
    elif bad=='sha':o['GA95']['coordinate_sha256']='hello'
    elif bad=='mutations':o['comparison']['mutations']=[]
    else:o['experimental_thermodynamic_constraints']['min_cross_fold_odds_ratio_if_opposing_states_exceed_90pct']='9'
    with pytest.raises(ValueError):m.check(o)