import pytest
from fractions import Fraction
from foldbounds.ga_gb_nmr import GA95,GB95,NMRModel
from foldbounds.real_fold_reporters import reporter_candidates,fret_ca_proxy_envelopes,invert_two_state_emission

def model(seq,dist_11_21):
    ca={i:(i,0,0) for i in range(1,57)}
    ca[11]=(0,0,0);ca[21]=(dist_11_21,0,0)
    return NMRModel(seq,ca,{})

def test_synthetic_pair_ranking_is_strict_and_no_mutation_sites():
    a=(model(GA95,14),model(GA95,16))
    b=(model(GB95,32),model(GB95,34))
    pairs=reporter_candidates(a,b,min_distance=12,max_distance=40)
    p=next(x for x in pairs if (x.i,x.j)==(11,21))
    assert p.worst_case_separation_A==16
    assert p.max_symmetric_CA_to_label_error_A==8
    assert p.state_with_shorter_distance=='GA95'
    assert all(x.i not in (20,30,45) and x.j not in (20,30,45) for x in pairs)

def test_fret_ca_is_proxy_not_calibration():
    r=reporter_candidates((model(GA95,14),model(GA95,16)),
                          (model(GB95,32),model(GB95,34)),12,40)
    x=next(x for x in r if (x.i,x.j)==(11,21))
    a=fret_ca_proxy_envelopes(x,25)['CA_distance_only_not_actual_dye_efficiency']
    assert a['GA95'][0]>a['GB95'][1]
    with pytest.raises(ValueError):fret_ca_proxy_envelopes(x,25,localization_A=-1)

@pytest.mark.parametrize('expected',[(Fraction(2,5),Fraction(3,5)),(Fraction(1,10),Fraction(9,10))])
def test_exact_two_state_population_inverse_bounds_contain_truth(expected):
    p=Fraction(1,2);alo,ahi=expected
    # A is high-emission state (80-90%), B low (10-20%).
    ylo=p*Fraction(4,5)+(1-p)*Fraction(1,10)
    yhi=p*Fraction(9,10)+(1-p)*Fraction(1,5)
    lo,hi=invert_two_state_emission(ylo,yhi,Fraction(4,5),Fraction(9,10),Fraction(1,10),Fraction(1,5))
    assert lo<=p<=hi

def test_three_state_or_ambiguous_emission_is_rejected():
    with pytest.raises(ValueError):invert_two_state_emission(Fraction(1,3),Fraction(2,3),Fraction(1,2),Fraction(4,5),Fraction(2,5),Fraction(3,5))
    with pytest.raises(ValueError):invert_two_state_emission(Fraction(1,100),Fraction(2,100),Fraction(4,5),Fraction(9,10),Fraction(1,10),Fraction(1,5))

from foldbounds.real_fold_reporters import invert_three_state_emission

def test_exact_three_state_lp_equals_two_state_when_no_unfolded():
    args=(Fraction(7,20),Fraction(11,20),Fraction(4,5),Fraction(9,10),Fraction(1,10),Fraction(1,5))
    x=invert_three_state_emission(*args,Fraction(0),Fraction(1),Fraction(0))
    a,b=invert_two_state_emission(*args)
    assert (x['fold_A_lower'],x['fold_A_upper'])==(a,b)
    assert x['unfolded_lower']==x['unfolded_upper']==0

def test_unfolded_state_uncertainty_cannot_sharpen_populations():
    args=(Fraction(2,5),Fraction(3,5),Fraction(4,5),Fraction(9,10),Fraction(1,10),Fraction(1,5),Fraction(0),Fraction(1))
    small=invert_three_state_emission(*args,Fraction(1,10))
    large=invert_three_state_emission(*args,Fraction(1))
    assert large['fold_A_lower']<=small['fold_A_lower']
    assert large['fold_A_upper']>=small['fold_A_upper']

def test_randomized_three_state_forward_inverse_containment():
    import random
    rng=random.Random(60241)
    for _ in range(400):
        pA=Fraction(rng.randrange(11),10)
        pU=Fraction(rng.randrange(11-int(pA*10)),10)
        pB=1-pA-pU
        # Emission order A high, B low, unfolded unknown.
        eA=Fraction(rng.randrange(8,10),10)
        eB=Fraction(rng.randrange(1,3),10)
        eU=Fraction(rng.randrange(11),10)
        y=pA*eA+pB*eB+pU*eU
        s=invert_three_state_emission(y,y,Fraction(4,5),Fraction(9,10),
                                      Fraction(1,10),Fraction(1,5),Fraction(0),Fraction(1),pU)
        assert s['fold_A_lower']<=pA<=s['fold_A_upper']
        assert s['unfolded_lower']<=pU<=s['unfolded_upper']

def test_incompatible_three_state_readout_fails_closed():
    with pytest.raises(ValueError):
        invert_three_state_emission(Fraction(0),Fraction(0),Fraction(4,5),Fraction(9,10),
                                    Fraction(1,10),Fraction(1,5),Fraction(1,2),Fraction(3,5),Fraction(1,10))