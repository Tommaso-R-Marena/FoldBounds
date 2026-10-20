import pytest,math
from fractions import Fraction
from foldbounds.ga_gb_nmr import GA95,GB95,MUTATIONS,NMRModel,analyze_pair,parse_pdb,contact_pairs,local_ca_geometry
from foldbounds.real_fold_odds import minimum_cross_fold_odds_ratio,robust_population,min_gap_shift_kcal, max_admissible_independent_state_energy_error

def synthetic_pdb(entry,seq,n=20):
    from foldbounds.ga_gb_nmr import RES3
    one_to_three={v:k for k,v in RES3.items()}
    header=f'HEADER    SYNTHETIC MODEL                          20-OCT-26   {entry}\nEXPDTA    SOLUTION NMR\n'
    lines=[header]
    for model in range(1,n+1):
        lines.append(f'MODEL     {model:4d}\n')
        for i,aa in enumerate(seq,1):
            atom=f'ATOM   {i:4d}   CA {one_to_three[aa]} A{i:4d}    {i*1.5:8.3f}{model*.05:8.3f}{(i%3)*.3:8.3f}{1.:6.2f}{20.:6.2f}           C\n'
            lines.append(atom)
        lines.append('ENDMDL\n')
    return ''.join(lines).encode()

def test_three_exact_canonical_real_mutations():
    assert len(GA95)==len(GB95)==56
    assert [(i,a,b) for i,(a,b) in enumerate(zip(GA95,GB95),1) if a!=b]==[(20,'L','A'),(30,'I','F'),(45,'L','Y')]

def test_synthetic_parser_multi_nmr_models():
    s=synthetic_pdb('2KDL',GA95)
    m=parse_pdb(s,'2KDL',GA95)
    assert len(m)==20 and all(len(t.ca)==56 for t in m)
    assert analyze_pair(s,synthetic_pdb('2KDM',GB95))['comparison']['mutations']==[
        {'position':20,'GA95':'L','GB95':'A'},{'position':30,'GA95':'I','GB95':'F'},{'position':45,'GA95':'L','GB95':'Y'}]

@pytest.mark.parametrize('what',['badid','wrongseq','nummodel','dropped_CA'])
def test_parser_fails_closed(what):
    p=synthetic_pdb('2KDL',GA95)
    if what=='badid':p=p.replace(b'2KDL',b'2KXX')
    if what=='wrongseq':p=synthetic_pdb('2KDL',GB95)
    if what=='nummodel':p=synthetic_pdb('2KDL',GA95,19)
    if what=='dropped_CA':p=p.replace(b'   CA ',b'   XX ',1)
    with pytest.raises(ValueError):parse_pdb(p,'2KDL',GA95)

@pytest.mark.parametrize('p,target,expected',[(Fraction(9,10),Fraction(9,10),0),
    (Fraction(19,20),Fraction(9,10),1), (Fraction(99,100),Fraction(9,10),1)])
def test_energy_error_tolerances(p,target,expected):
    d=max_admissible_independent_state_energy_error(p,target)
    assert (d>0)==bool(expected)
    if expected:assert robust_population(p,Fraction(1,1000000),Fraction(1000000))[0] < target

def test_actual_90pct_opposing_folds_imply_at_least_81_odds_ratio():
    assert minimum_cross_fold_odds_ratio(Fraction(9,10))==81
    assert 2.5<min_gap_shift_kcal(Fraction(9,10))<2.7

def test_synthetic_contact_occupancy_is_not_boltzmann_population():
    m=parse_pdb(synthetic_pdb('2KDL',GA95),'2KDL',GA95)
    assert all(isinstance(contact_pairs(x),set) for x in m)
    assert all(0<=local_ca_geometry(x)['helix_like_i_to_i3_fraction']<=1 for x in m)