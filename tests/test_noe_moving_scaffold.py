"""Independent assumptions and fail-closed checks for moving-scaffold NOE no-go."""
from __future__ import annotations
import json
from dataclasses import replace
from pathlib import Path
from fractions import Fraction
from math import dist
from random import Random
import pytest

from foldbounds.noe_scaffold_dual import verify_no_go
from foldbounds.noe_scaffold_branch_bound import default_val69_branch_bound
from foldbounds.noe_moving_scaffold import (
    HeterogeneousScaffoldBudget, UniformMovingScaffoldProof,
    verify_uniform_moving_scaffold,verify_heterogeneous_moving_scaffold,
)

RAW=json.loads((Path(__file__).resolve().parents[1]/'data'/'val69_local_noe_dual.json').read_text())
CERT=RAW


def test_uniform_continuum_and_boundary():
    bb=default_val69_branch_bound()
    for anchor, local in [(0,199),(10,189),(50,149),(100,99),(150,49),(199,0)]:
        r=verify_uniform_moving_scaffold(RAW,UniformMovingScaffoldProof(anchor,local,bb))
        assert r is not None and r.effective_fixed_frame_delta_mA==199
        assert r.branch_result.objective_upper_mA2 < r.branch_result.necessary_lower_mA2
    assert verify_uniform_moving_scaffold(RAW,UniformMovingScaffoldProof(50,150,bb)) is None


def test_uniform_rejects_bad_budgets_and_proofs():
    bb=default_val69_branch_bound()
    for bad in [-1,True,1.5,'10',None]:
        assert verify_uniform_moving_scaffold(RAW,UniformMovingScaffoldProof(bad,0,bb)) is None
    assert verify_uniform_moving_scaffold(RAW,UniformMovingScaffoldProof(0,0,replace(bb,h=replace(bb.h,cap_mA2=10**20)))) is None
    forged=json.loads(json.dumps(RAW));forged['processed_restraint_upper_mA']=6000
    assert verify_uniform_moving_scaffold(forged,UniformMovingScaffoldProof(0,0,bb)) is None


def test_heterogeneous_missing_nonproof_anchors_is_safe():
    # At δ=0.05 Å, each of five *used* anchors may shift by 0.04 Å.
    # The other two source-model heavy atoms can be entirely unconstrained.
    used=['N','CA','CB','CG2','O']
    b=HeterogeneousScaffoldBudget(
        {k:40 for k in used},
        {label:{k:50 for k in used} for label in ('H','HB')},
    )
    result=verify_heterogeneous_moving_scaffold(CERT,b)
    assert result is not None
    assert result.constrained_anchor_names==tuple(sorted(used))
    assert result.possible_support_mA2 < result.necessary_support_mA2
    assert 'C' not in result.constrained_anchor_names and 'CG1' not in result.constrained_anchor_names


def test_heterogeneous_fails_closed_and_mixed_budgets():
    used=['N','CA','CB','CG2','O']
    a={k:40 for k in used}
    b={s:{k:50 for k in used} for s in ('H','HB')}
    assert verify_heterogeneous_moving_scaffold(CERT,HeterogeneousScaffoldBudget(a,b)) is not None
    assert verify_heterogeneous_moving_scaffold(CERT,HeterogeneousScaffoldBudget({k:200 for k in used},b)) is None
    assert verify_heterogeneous_moving_scaffold(CERT,HeterogeneousScaffoldBudget({k:40 for k in used if k!='N'},b)) is None
    assert verify_heterogeneous_moving_scaffold(CERT,HeterogeneousScaffoldBudget({k:(True if k=='N' else 40) for k in used},b)) is None
    broken={s:dict(t) for s,t in b.items()};broken['HB'].pop('O')
    assert verify_heterogeneous_moving_scaffold(CERT,HeterogeneousScaffoldBudget(a,broken)) is None
    # A nonzero weight mutation must still reproduce the exact dual identity.
    false_dual=json.loads(json.dumps(CERT));false_dual['proof']['H'][0][3]+=1
    assert verify_heterogeneous_moving_scaffold(false_dual,HeterogeneousScaffoldBudget(a,b)) is None


def test_exact_triangle_transfer_independent_random_oracle():
    rng=Random(2718)
    # Independently sample completely arbitrary real anchor motions and H motions.
    # The derived old-reference shell inequality must hold in each instance.
    for _ in range(500):
        anchor=[rng.uniform(-10,10) for _ in range(3)]
        old=[rng.uniform(-10,10) for _ in range(3)]
        new=[rng.uniform(-10,10) for _ in range(3)]
        shifted=[anchor[i]+rng.uniform(-.5,.5) for i in range(3)]
        eta=dist(anchor,shifted)
        delta=abs(dist(new,shifted)-dist(old,anchor))
        effective=abs(dist(new,anchor)-dist(old,anchor))
        assert effective <= eta+delta+1e-10


def test_heterogeneous_is_not_physical_claim():
    # Deliberately omit local-distance assumptions for one actually used anchor:
    # neither method can infer a conclusion from missing experimental inputs.
    b=HeterogeneousScaffoldBudget({}, {'H':{},'HB':{}})
    assert verify_heterogeneous_moving_scaffold(CERT,b) is None
