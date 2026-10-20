from math import sqrt
from foldbounds.ga_gb_nmr import NMRModel,GA95,GB95
from foldbounds.real_aromatic_geometry import classify_ring_pair,ring_normal,jaccard,network_discrimination,aromatic_summary

def model_for_test(seq=GB95,split=5.0):
    # Unit geometry controls; not physical NMR coordinates.
    side={}
    for i,x in [(29,0.0),(30,split)]:
        side[i]={'CG':(x,0.0,0.0),'CD1':(x+1,0.0,0.0),'CD2':(x,1,0.0),
                 'CE1':(x+1,1,0.0),'CE2':(x+.5,1.5,0.0),'CZ':(x+.2,1.0,0.0)}
    ca={i:(1.5*i,0,0) for i in range(1,57)}
    return NMRModel(seq,ca,side)

def test_ring_normal_and_proximity_only_when_physically_possible():
    m=model_for_test()
    n=ring_normal(m,29)
    assert n is not None and abs(n[2])>0.999
    d=classify_ring_pair(m,29,30)
    assert d is not None and 4<d['centroid_distance_A']<6
    assert d['orientation_class']=='close_other_orientation' # coplanar, not stacked
    assert ring_normal(model_for_test(seq=GA95),30) is None

def test_contact_network_jaccard_math_and_degenerate_case():
    assert jaccard({(1,5)}, {(1,5),(2,6)})==.5
    assert jaccard(set(),set())==1

def test_robust_model_separation_is_not_same_as_fold_population():
    a=model_for_test();b=model_for_test()
    d=network_discrimination((a,a),(b,b))
    assert d['cross_fold_pairs']==4 and not d['strict_model_level_separation']

def test_aromatic_absent_vs_present():
    assert aromatic_summary((model_for_test(seq=GA95),))[0]['valid_models']==0
    assert aromatic_summary((model_for_test(seq=GB95),))[0]['valid_models']==1