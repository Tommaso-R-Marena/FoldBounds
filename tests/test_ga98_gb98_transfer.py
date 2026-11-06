"""Adversarial unit controls for mutation-heldout 3D NMR fold transfer."""
import json
from dataclasses import replace
from fractions import Fraction

import pytest

from foldbounds.ga_gb_nmr import NMRModel
from foldbounds.ga_gb_error_codes import CodeCertificate, Reporter
from experiments.ga98_gb98_transfer import (
    GA98, GB98, _geometric_class, heldout_transfer, verify_frozen_reference,
    trusted_test_structures,
)


def artificial_code():
    # Three coordinate-disjoint, direction-consistent geometric reporter
    # pairs. This is SYNTHETIC and never substituted for real PDB data.
    sites=((5,25),(10,21),(11,22))
    pairs=tuple(Reporter(a,b,10000,"GA95",(12000,14000),(24000,26000))
                for a,b in sites)
    return CodeCertificate(pairs,3,1,10000,(),(),len(pairs),{
        "2KDL":"fb75cc55eb8b89e0cec5e2be69901b348e60916375b9c410a0573edfca1db591",
        "2KDM":"1c49e5faa25b443285432d935c621151294ed3b636b0106e6b1a7f7ccf2d00ec"
    },14,20,Fraction(5000))


def make_structure(sequence,shorter=True,flip=()):
    ca={i:(float(i),0.,0.) for i in range(1,57)}
    for a,b in ((5,25),(10,21),(11,22)):
        # Use distinct residue y strips so unrelated labels never collide.
        ca[a]=(0.,float(a*3),0.)
        distance=26. if (a,b) in flip else (13. if shorter else 25.)
        ca[b]=(distance,float(a*3),0.)
    return NMRModel(sequence,ca,{})


def test_test_sequences_have_one_exact_helix_beta_mutation():
    assert len(GA98)==len(GB98)==56
    diff=[(i+1,a,b) for i,(a,b) in enumerate(zip(GA98,GB98)) if a!=b]
    assert diff==[(45,"L","Y")]
    assert GA98!=GB98


def test_independent_geometry_voting_on_reference_and_unseen():
    cert=artificial_code()
    for seq,short,truth in ((GA98,True,"GA95"),(GB98,False,"GB95")):
        model=make_structure(seq,short)
        got,readouts,margins,votes=_geometric_class(cert,model)
        assert got==truth
        assert votes[truth]==3
        assert len(readouts)==3
        assert all(x>=0 for x in margins.values())


def test_nontrivial_transfer_error_is_not_suppressed():
    cert=artificial_code()
    a_good=make_structure(GA98,True)
    b_good=make_structure(GB98,False)
    results=heldout_transfer({3:cert},(a_good,)*20,(b_good,)*20)
    correct=results["3"]
    assert correct["transfer_clean_correct"]==40
    assert correct["transfer_corrupted_trials_correct"]==160
    assert correct["minimum_correct_votes_per_model"]==3
    assert correct["minimum_signed_reporter_decision_margin_A"]>0

    # One GA-family reporter becomes GB-like in unseen geometry.
    flip_one=make_structure(GA98,True,flip=((5,25),))
    results2=heldout_transfer({3:cert},(flip_one,)*20,(b_good,)*20)
    got=results2["3"]
    assert got["transfer_clean_correct"]==40
    assert got["transfer_corrupted_trials_correct"]<160
    assert got["minimum_correct_votes_per_model"]==2
    assert got["minimum_signed_reporter_decision_margin_A"]<0


def test_full_unseen_fold_failure_is_reported():
    cert=artificial_code()
    ga_wrong=make_structure(GA98,True,flip=((5,25),(10,21)))
    b=make_structure(GB98,False)
    res=heldout_transfer({3:cert},(ga_wrong,)*20,(b,)*20)["3"]
    assert res["transfer_clean_correct"]==20
    assert res["by_entry"]["2LHC"]["clean_correct"]==0
    assert res["by_entry"]["2LHD"]["clean_correct"]==20


def test_ambiguous_sub_milliangstrom_threshold_fails_closed():
    cert=artificial_code()
    # First pair midpoint is 19 Å. Exact equality cannot be assigned
    # a binary class by a positive strict margin certificate.
    path=make_structure(GA98,True)
    ca=dict(path.ca)
    ca[25]=(19.,15.,0.)
    with pytest.raises(ValueError,match="ambiguous|precision"):
        _geometric_class(cert,replace(path,ca=ca))


def test_untrusted_training_or_target_bytes_rejected(tmp_path):
    ref=tmp_path/"fake.json"
    ref.write_text(json.dumps({"PDB_ensemble_sources":{"2KDL":"fake"},
                               "codes":{str(k):{} for k in (1,3,5,7,9)}}))
    with pytest.raises(ValueError):
        verify_frozen_reference(b"fake",b"fake",ref)
    with pytest.raises(ValueError):
        trusted_test_structures(b"missing",b"missing")
