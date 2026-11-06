"""Independent tests for source-only triangle-metric inference."""
from pathlib import Path
from dataclasses import replace
import pytest

from foldbounds.noe_restraint_audit import AtomRef,NOERestraint,load_pinned_reference
from foldbounds.noe_metric_closure import (
    discover_proximities,verify_proximity_witness,ProximityWitness)

FIXTURE=Path(__file__).resolve().parents[1]/'data/pinned_unambiguous_noe_2png_6nbn.json'


def make(k,a,b,upper):
    return NOERestraint(k,AtomRef(a,'H'),(a,0,0),
                        AtomRef(b,'H'),(b,0,0),0,upper)


def test_proof_carrying_unmeasured_contact():
    data=[make(1,1,2,3000),make(2,2,3,4000)]
    proposals=discover_proximities(data,cutoff_mA=8000)
    assert len(proposals)==1
    w=proposals[0]
    assert (w.restraint_ids,w.upper_bound_mA)==((1,2),7000)
    assert verify_proximity_witness(data,w)
    assert not verify_proximity_witness(data,replace(w,upper_bound_mA=6999))
    assert not verify_proximity_witness(data,replace(w,restraint_ids=(2,1)))
    assert not verify_proximity_witness(data,w,cutoff_mA=6000)


def test_no_direct_restatement():
    data=[make(1,1,2,3000),make(2,2,3,4000),make(3,1,3,7800)]
    assert discover_proximities(data,cutoff_mA=8000)==()
    w=ProximityWitness((AtomRef(1,'H'),AtomRef(3,'H')),(3,),7800)
    assert not verify_proximity_witness(data,w)
    assert verify_proximity_witness(data,w,require_unmeasured=False)


def test_full_nmr_real_pair_closure():
    if not FIXTURE.exists():
        pytest.skip("source-pinned fixture on PR branch")
    data=load_pinned_reference(FIXTURE)
    expectations={'2png':(110,1),'6nbn':(359,39)}
    for name,(n,long_range) in expectations.items():
        witnesses=discover_proximities(data[name],cutoff_mA=8000)
        assert len(witnesses)==n
        assert sum(w.sequence_separation>=5 for w in witnesses)==long_range
        assert len(set(frozenset(w.pair) for w in witnesses))==len(witnesses)
        assert all(verify_proximity_witness(data[name],w) for w in witnesses)
