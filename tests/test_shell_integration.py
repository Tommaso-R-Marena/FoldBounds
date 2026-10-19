"""Integration and audit-replay checks for independently counted shell masses."""
from dataclasses import replace
from fractions import Fraction

import pytest

from foldbounds.lattice import ROOT, build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead
from foldbounds.adaptive import adaptive_refine, verify_refinement_run


@pytest.mark.parametrize('seq,q,t,core', [
    ('H'*8,Fraction(2),1,6),
    ('HPHHPHHP',Fraction(3,2),1,6),
    ('H'*12,Fraction(2),2,10),
    ('HPHHPHHPHHHH',Fraction(3,2),2,10),
])
def test_shell_strengthening_is_sound_against_exact_independent_spectrum(seq,q,t,core):
    cert=build_frontier(seq,q,t,2)
    old=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                         use_monotone=True,use_turn_family=True,transfer_memory=4)
    new=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                         use_monotone=True,use_turn_family=True,transfer_memory=4,
                         shell_core=core)
    dos=exact_spectrum(seq)
    exact_e=sum((q**k*v for k,v in dos.items() if k>=t),Fraction(0))
    exact_n=sum((q**k*v for k,v in dos.items() if k<t),Fraction(0))
    assert old.event_lower<=new.event_lower<=exact_e<=new.event_upper==old.event_upper
    assert old.complement_lower<=new.complement_lower<=exact_n<=new.complement_upper==old.complement_upper
    assert old.partition_lower<=new.partition_lower<=exact_e+exact_n<=new.partition_upper
    assert new.partition_upper==old.partition_upper
    actual=exact_e/(exact_e+exact_n)
    assert old.probability_lower<=new.probability_lower<=actual<=new.probability_upper<=old.probability_upper


def test_mixed_chain_event_lower_mass_has_strict_gain():
    s='HPHHPHHPHHHH'
    cert=build_frontier(s,Fraction(3,2),2,2)
    old=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                         use_monotone=True,use_turn_family=True,transfer_memory=4)
    new=verify_lookahead(cert,lookahead=0,use_matching=True,use_transfer=True,
                         use_monotone=True,use_turn_family=True,transfer_memory=4,
                         shell_core=10)
    assert old.event_lower>=Fraction(459,8)
    assert new.event_lower>=6867
    assert new.event_lower>old.event_lower
    assert new.probability_lower>old.probability_lower


def test_adaptive_shell_mode_is_independently_replayed_and_tampering_rejected():
    cert=build_frontier('H'*10,Fraction(2),1,2)
    run=adaptive_refine(cert,steps=2,lookahead=0,use_matching=True,
                        use_transfer=False,use_monotone=True,
                        use_turn_family=True,shell_core=8,policy='mass')
    assert run.shell_core==8
    assert verify_refinement_run(run)==run.bounds
    with pytest.raises((ValueError,AssertionError)):
        verify_refinement_run(replace(run,shell_core=None))


@pytest.mark.parametrize('bad',[True,False,1,0,-1,2.5,'8',None])
def test_shell_core_rejected_when_invalid(bad):
    if bad is None:
        return
    cert=build_frontier('H'*8,Fraction(2),1,2)
    with pytest.raises(ValueError):
        verify_lookahead(cert,shell_core=bad)