"""Sparse NOE repair: rigorous minimum count and minimal-count radius bracket."""
from pathlib import Path
from dataclasses import replace
import hashlib,json
import pytest
from foldbounds.noe_restraint_audit import AtomRef,NOERestraint,load_pinned_reference
from foldbounds.noe_sparse_repair import (SparseRepairCertificate,
    verify_sparse_repair,sparse_from_json)
from foldbounds.noe_minimax_repair import RepairCertificate

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'data/pinned_unambiguous_noe_2png_6nbn.json'
WIT=ROOT/'data/noe_minimal_support_repair_witnesses.json'


def test_one_edge_requires_exactly_one_atom_and_full_gap():
    a=AtomRef(1,'H');b=AtomRef(2,'H')
    data=(NOERestraint(1,a,(0,0,0),b,(10,0,0),0,8),)
    c=RepairCertificate('toy',1,1,999,2000,((a,(2000,0,0)),))
    s=SparseRepairCertificate(c,(1,),1,1999)
    out=verify_sparse_repair(data,s)
    assert out is not None
    assert out.minimum_number_of_atoms_moved==1
    assert out.sparse_upper_radius_ticks==2000
    assert verify_sparse_repair(data,replace(s,strict_sparse_lower_radius_ticks=2000)) is None
    assert verify_sparse_repair(data,replace(s,disjoint_violated_ids=(1,1))) is None
    assert verify_sparse_repair(data,replace(s,lower_witness_id=2)) is None
    assert verify_sparse_repair(data,replace(s,repair=replace(c,moves=()))) is None


def test_matching_validator_rejects_overlapping_nontrivial_witness():
    aa=AtomRef(1,'H');bb=AtomRef(2,'H');cc=AtomRef(3,'H')
    records=(NOERestraint(1,aa,(0,0,0),bb,(10,0,0),0,8),
             NOERestraint(2,bb,(10,0,0),cc,(20,0,0),0,8))
    # Two matching IDs that SHARE bb cannot certify 2 distinct changed atoms.
    c=RepairCertificate('toy',2,1,999,2000,
          ((aa,(2000,0,0)),(cc,(-2000,0,0))))
    s=SparseRepairCertificate(c,(1,2),1,1999)
    assert verify_sparse_repair(records,s) is None


@pytest.mark.parametrize('sample,expected',[
 ('2png',(304,170,5,54567,54569)),
 ('6nbn',(553,381,4,737441,737443)),
])
def test_source_pinned_minimum_support_real_data(sample,expected):
    if not REF.exists():pytest.skip('requires inherited source-pinned PR #80 fixture')
    original=REF.read_bytes()
    gitsha=hashlib.sha1(f'blob {len(original)}\0'.encode()+original).hexdigest()
    j=json.loads(WIT.read_text())
    assert gitsha==j['source_github_file_blob_sha']
    records=load_pinned_reference(REF)[sample]
    cert=sparse_from_json(sample,j['samples'][sample])
    verified=verify_sparse_repair(records,cert)
    assert verified is not None
    assert (verified.checked_repair.restraints,verified.checked_repair.atoms,
            verified.minimum_number_of_atoms_moved,
            verified.strict_sparse_lower_radius_ticks,
            verified.sparse_upper_radius_ticks)==expected
    assert verified.sparse_upper_radius_ticks-verified.strict_sparse_lower_radius_ticks<=2
    assert verify_sparse_repair(records,replace(cert,disjoint_violated_ids=cert.disjoint_violated_ids[:-1])) is None
    assert verify_sparse_repair(records,replace(cert,strict_sparse_lower_radius_ticks=cert.repair.upper_radius_ticks)) is None
    assert verify_sparse_repair(records,replace(cert,repair=replace(cert.repair,moves=cert.repair.moves[:-1]))) is None


def test_source_schema_witness_is_independent_of_optimizer():
    j=json.loads(WIT.read_text())
    assert j['schema']=='FoldBounds.NOESparseRepair.v1'
    for sample,v in j['samples'].items():
        cert=sparse_from_json(sample,v)
        assert cert.lower_witness_id in cert.disjoint_violated_ids
        assert len(cert.disjoint_violated_ids)==len(cert.repair.moves)
        assert verify_sparse_repair([],cert) is None
        with pytest.raises(ValueError):sparse_from_json(sample,{**v,'matching_ids':None})
