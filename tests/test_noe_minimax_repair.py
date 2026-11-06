"""Adversarial, independent integer replay of real NOE minimax proposals."""
from dataclasses import replace
from pathlib import Path
from math import isqrt
import hashlib,json
import pytest

from foldbounds.noe_restraint_audit import AtomRef,NOERestraint,load_pinned_reference
from foldbounds.noe_minimax_repair import (
    RepairCertificate,certificate_from_json,verify_repair,SCALE)
from foldbounds.noe_holdout import frozen_fold,five_fold_audit

ROOT=Path(__file__).resolve().parents[1]
WITNESSES=ROOT/'data/noe_minimax_repair_witnesses.json'
REFERENCE=ROOT/'data/pinned_unambiguous_noe_2png_6nbn.json'


def make_single(hi=8,lo=0):
    return NOERestraint(1,AtomRef(1,'H'),(0,0,0),AtomRef(2,'H'),(10,0,0),lo,hi)


def test_exact_synthetic_upper_violation_minimax_bracket():
    # Source separation 10 mÅ, allowed max 8 => optimum radius 1 mÅ.
    r=make_single()
    cert=RepairCertificate('toy',1,1,999,1000,
        ((r.left,(1000,0,0)),(r.right,(-1000,0,0))))
    verified=verify_repair([r],cert)
    assert verified is not None
    assert verified.actual_max_movement_squared_ticks==1000**2
    assert verified.optimality_gap_angstrom==1e-6
    assert verify_repair([r],replace(cert,strict_lower_radius_ticks=1000)) is None
    assert verify_repair([r],replace(cert,upper_radius_ticks=999)) is None
    assert verify_repair([r],replace(cert,moves=cert.moves[:1])) is None
    assert verify_repair([r],replace(cert,moves=cert.moves+cert.moves[:1])) is None
    assert verify_repair([r],replace(cert,lower_restraint_id=2)) is None


def test_exact_synthetic_below_lower_and_unrelated_constraint():
    left=AtomRef(1,'H');right=AtomRef(2,'H')
    r=NOERestraint(1,left,(0,0,0),right,(2,0,0),4,10)
    c=RepairCertificate('toy',1,1,999,1000,
        ((left,(-1000,0,0)),(right,(1000,0,0))))
    assert verify_repair((r,),c) is not None
    extra=NOERestraint(2,left,(0,0,0),AtomRef(3,'H'),(0,0,0),0,0)
    assert verify_repair((r,extra),replace(c,source_count=2)) is None
    assert verify_repair((r,extra),replace(c,source_count=2,
                                moves=c.moves+((AtomRef(999,'H'),(1,0,0)),))) is None


def test_strict_input_parsing_and_errors():
    j=json.loads(WITNESSES.read_text())
    assert j['schema']=='FoldBounds.NOEMinimaxRepair.v1'
    for id,obj in j['samples'].items():
        c=certificate_from_json(id,obj)
        assert c.source_count in (304,553)
        assert c.upper_radius_ticks-c.strict_lower_radius_ticks<=2
        with pytest.raises(ValueError):
            certificate_from_json(id,{**obj,'moves':[['bad/no',1,0,0]]})
        with pytest.raises(ValueError):
            certificate_from_json(id,{**obj,'moves':[["1:H",1,2]]})
        assert verify_repair([],c) is None


def test_split_is_deterministic_and_independent_of_coordinates():
    ids=[frozen_fold('2png',i) for i in range(1,100)]
    assert set(ids)==set(range(5))
    assert ids==[frozen_fold('2png',i) for i in range(1,100)]
    with pytest.raises(ValueError):frozen_fold('',1)
    with pytest.raises(ValueError):frozen_fold('x',True)


def test_negative_holdout_toy():
    aa=AtomRef(1,'H');bb=AtomRef(2,'H');cc=AtomRef(3,'H')
    test=[NOERestraint(i+1,a,(0,0,0),b,(0,0,0),1,8000)
          for i,(a,b) in enumerate([(aa,bb),(bb,cc),(aa,cc),(aa,bb),(aa,cc)])]
    # Duplicate pairs have separate valid experimental IDs but metric path
    # verification will refuse direct observed repetition for queried pairs.
    assert len(test)==5
    assert isinstance(frozen_fold('toy',1),int)


@pytest.mark.parametrize('sample,expected',[
    ('2png',(304,170,27283,27284,10)),
    ('6nbn',(553,381,368720,368722,9)),
])
def test_real_exact_integer_global_minimax_certificate(sample,expected):
    if not REFERENCE.exists():pytest.skip('requires source-pinned PR #80 fixture')
    fixture_bytes=REFERENCE.read_bytes()
    digest=hashlib.sha1((f'blob {len(fixture_bytes)}\0').encode()+fixture_bytes).hexdigest()
    j=json.loads(WITNESSES.read_text())
    assert digest==j['source_github_file_blob_sha']
    ref=load_pinned_reference(REFERENCE)
    cert=certificate_from_json(sample,j['samples'][sample])
    out=verify_repair(ref[sample],cert)
    assert out is not None
    assert (out.restraints,out.atoms,out.strict_lower_radius_ticks,
            out.upper_radius_ticks,out.moved_atoms)==expected
    assert out.optimality_gap_angstrom<=2e-6+1e-12
    # Deliberate tampering must be rejected.
    assert verify_repair(ref[sample],replace(cert,upper_radius_ticks=cert.strict_lower_radius_ticks)) is None
    assert verify_repair(ref[sample],replace(cert,moves=cert.moves[:-1])) is None
    assert verify_repair(ref[sample],replace(cert,source_count=1)) is None
    x=list(cert.moves);a,delta=x[0];x[0]=(a,(delta[0]+50000,delta[1],delta[2]))
    assert verify_repair(ref[sample],replace(cert,moves=tuple(x))) is None


@pytest.mark.parametrize('sample,expected',[
 ('2png',(304,257,93,3,0)),('6nbn',(553,348,87,2,0))
])
def test_frozen_five_fold_real_noe_upper_only_holdout(sample,expected):
    if not REFERENCE.exists():pytest.skip('requires source-pinned PR #80 fixture')
    r=load_pinned_reference(REFERENCE)[sample]
    folds=five_fold_audit(sample,r)
    actual=(sum(x.heldout for x in folds),sum(x.connected for x in folds),
            sum(x.proximity_le_8a for x in folds),
            sum(x.inferred_upper_stronger_than_heldout for x in folds),
            sum(x.entire_interval_guaranteed_by_upper_only for x in folds))
    assert actual==expected


def test_independent_untrusted_proposal_must_pass_exact_checker():
    import sys
    sys.path.insert(0,str(ROOT/'experiments'))
    from propose_noe_minimax_repair import propose_repair
    r=make_single()
    proposed=propose_repair('toy',(r,))
    assert verify_repair((r,),proposed) is not None
    with pytest.raises(ValueError):propose_repair('toy',())


@pytest.mark.parametrize('sample',('2png','6nbn'))
def test_real_generator_reproduces_verified_minimax_bracket(sample):
    if not REFERENCE.exists():pytest.skip('requires source-pinned PR #80 fixture')
    import sys
    sys.path.insert(0,str(ROOT/'experiments'))
    from propose_noe_minimax_repair import propose_repair
    ref=load_pinned_reference(REFERENCE)
    proposed=propose_repair(sample,ref[sample])
    assert verify_repair(ref[sample],proposed) is not None
    # The accepted checked certificate is held to an independently useful
    # finite objective bound; the producer is not trusted for optimality.
    assert proposed.upper_radius_ticks-proposed.strict_lower_radius_ticks <= 80
