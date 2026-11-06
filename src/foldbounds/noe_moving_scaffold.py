"""Checked moving-anchor generalization of the source-pinned Val69 NOE no-go.

Each anchor may move by eta, and each new hydrogen-to-moving-anchor
separation may deviate from its source value by delta. Triangle inequality
makes the new hydrogen-to-ORIGINAL-anchor deviation at most eta+delta.
All inputs are stipulated millangstrom budgets, not experimental error bars.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from typing import Mapping

from .noe_scaffold_dual import (
    _ceil, _dot, _sub, verify_no_go, with_delta, _validate,
)
from .noe_scaffold_branch_bound import (
    BranchBoundProof, BranchBoundResult, verify_branch_bound_no_go,
)


def nat(x):
    return type(x) is int and x >= 0


@dataclass(frozen=True)
class UniformMovingScaffoldProof:
    max_anchor_motion_mA: int
    max_hydrogen_anchor_distance_change_mA: int
    fixed_frame_proof: BranchBoundProof


@dataclass(frozen=True)
class VerifiedUniformMovingScaffoldNoGo:
    anchor_motion_mA: int
    local_distance_change_mA: int
    effective_fixed_frame_delta_mA: int
    branch_result: BranchBoundResult


def verify_uniform_moving_scaffold(raw: dict, proof: UniformMovingScaffoldProof):
    if type(proof) is not UniformMovingScaffoldProof:
        return None
    eta=proof.max_anchor_motion_mA
    delta=proof.max_hydrogen_anchor_distance_change_mA
    if not nat(eta) or not nat(delta) or type(proof.fixed_frame_proof) is not BranchBoundProof:
        return None
    if eta+delta>proof.fixed_frame_proof.delta_mA:
        return None
    checked=verify_branch_bound_no_go(raw,proof.fixed_frame_proof)
    if checked is None:
        return None
    return VerifiedUniformMovingScaffoldNoGo(eta,delta,eta+delta,checked)


@dataclass(frozen=True)
class HeterogeneousScaffoldBudget:
    # Missing anchor budgets mean UNCONSTRAINED, never zero.
    anchor_motion_mA: Mapping[str,int]
    hydrogen_anchor_change_mA: Mapping[str,Mapping[str,int]]


@dataclass(frozen=True)
class VerifiedHeterogeneousScaffoldNoGo:
    necessary_support_mA2: int
    possible_support_mA2: Fraction
    h_support_mA2: Fraction
    hb_support_mA2: Fraction
    constrained_anchor_names: tuple[str,...]


def _dual_support(endpoint,direction,anchors,rows,local,movement):
    if type(rows) is not list or not 1<=len(rows)<=50:
        return None
    normal=[Fraction(0)]*3
    cost=Fraction(0)
    used=set()
    for row in rows:
        if (type(row) is not list or len(row)!=5 or type(row[0]) is not int
            or row[0] not in (-1,1) or type(row[1]) is not str
            or type(row[2]) is not str or row[1]==row[2]
            or not nat(row[3]) or type(row[4]) is not int
            or row[4]<=0):
            return None
        pa=anchors.get(row[1])
        pb=anchors.get(row[2])
        if pa is None or pb is None:
            return None
        e=[]
        for name,point in ((row[1],pa),(row[2],pb)):
            if name not in local or name not in movement:
                return None
            if not nat(local[name]) or not nat(movement[name]):
                return None
            used.add(name)
            d=local[name]+movement[name]
            p=_sub(endpoint,point)
            e.append(2*_ceil(_dot(p,p))*d+d*d)
        weight=Fraction(row[3],row[4])
        n=_sub(pa,pb)
        for j in range(3):
            normal[j]+=weight*row[0]*n[j]
        cost+=weight*Fraction(e[0]+e[1],2)
    if tuple(normal)!=tuple(direction):
        return None
    return cost,used


def verify_heterogeneous_moving_scaffold(cert: dict,
                                         budget: HeterogeneousScaffoldBudget):
    try:
        if (not _validate(cert) or
            type(budget) is not HeterogeneousScaffoldBudget or
            not isinstance(budget.anchor_motion_mA,Mapping) or
            not isinstance(budget.hydrogen_anchor_change_mA,Mapping) or
            set(budget.hydrogen_anchor_change_mA)!={'H','HB'} or
            any(not isinstance(x,Mapping) for x in budget.hydrogen_anchor_change_mA.values())):
            return None
        # Replays all input-vector identities and the base source conflict.
        if verify_no_go(with_delta(cert,0)) is None:
            return None
        v=_sub(cert['target_hydrogens_val69']['HB'],cert['target_hydrogens_val69']['H'])
        source=_dot(v,v)-cert['processed_restraint_upper_mA']*_ceil(_dot(v,v))
        if source<=0:
            return None
        a=_dual_support(cert['target_hydrogens_val69']['H'],v,cert['fixed_heavy_atoms_val69'],cert['proof']['H'],
                        budget.hydrogen_anchor_change_mA['H'],budget.anchor_motion_mA)
        b=_dual_support(cert['target_hydrogens_val69']['HB'],tuple(-x for x in v),
                        cert['fixed_heavy_atoms_val69'],cert['proof']['HB'],
                        budget.hydrogen_anchor_change_mA['HB'],budget.anchor_motion_mA)
        if a is None or b is None or a[0]+b[0]>=source:
            return None
        return VerifiedHeterogeneousScaffoldNoGo(
            source,a[0]+b[0],a[0],b[0],tuple(sorted(a[1]|b[1])))
    except (TypeError,ValueError,KeyError,AttributeError,OverflowError,ZeroDivisionError):
        return None
