"""Locally checked, state-dependent JEPA-guided Metropolis-Hastings.

This does NOT enumerate the conformational universe for a single step.
Any FROZEN deterministic policy returning normalized positive rational
pivot-class weights is admissible. The locally checked Hastings ratio
preserves the target Boltzmann measure algebraically. Policy changes during
a trajectory invalidate the time-homogeneous stationary-kernel argument.
This is Monte Carlo, not molecular dynamics or a physical folding rate.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite, lcm
import random
from typing import Callable

from .lattice import (
    Path, _check_q, _check_sequence, _contacts, _possible_hh_pairs,
    _validate_path,
)
from .jepa_reversible_kinetics import all_states, move_classes, pivot, _angles

Policy = Callable[[Path], tuple[Fraction,...]]


def quantize_ranks(scores: tuple[float,...]) -> tuple[Fraction,...]:
    """Exact positive rank weights without trusting float exponentials."""
    if not scores or any(
        type(score) not in (int,float) or not isfinite(score)
        for score in scores
    ):
        raise ValueError("invalid learned score vector")
    ordering=sorted(range(len(scores)),key=lambda i:(-scores[i],i))
    weights=[0]*len(scores)
    for rank,index in enumerate(ordering):
        weights[index]=len(scores)-rank
    z=sum(weights)
    return tuple(Fraction(weight,z) for weight in weights)


def validate_class_probs(path:Path,probs:tuple[Fraction,...])->None:
    if (len(probs)!=len(move_classes(len(path)))
        or any(type(value) is not Fraction or value<=0 for value in probs)
        or sum(probs,Fraction())!=1):
        raise ValueError("frozen class weights must be exact positive Fractions summing to 1")


class StabilityGuard:
    """Detect observed policy drift for any previously queried conformation.

    Caching alone is not a universal proof of purity on unseen states.
    The caller must also freeze the complete neural model and settings.
    """
    def __init__(self, producer:Policy):
        self.producer=producer
        self.seen:dict[Path,tuple[Fraction,...]]={}

    def __call__(self,path:Path)->tuple[Fraction,...]:
        proposed=tuple(self.producer(path))
        validate_class_probs(path,proposed)
        old=self.seen.get(path)
        if old is not None and old!=proposed:
            raise ValueError("policy changed on a previously queried conformation")
        self.seen[path]=proposed
        return proposed


def proposal_row(path:Path,policy:Policy)->dict[Path,Fraction]:
    """Exact class-weighted pivot proposal; collisions become self-loops."""
    classes=move_classes(len(path))
    probs=tuple(policy(path))
    validate_class_probs(path,probs)
    result:dict[Path,Fraction]={}
    for (k,family),class_prob in zip(classes,probs):
        angles=_angles(family)
        for angle in angles:
            target=pivot(path,k,angle)
            if target is None:
                target=path
            result[target]=result.get(target,Fraction())+class_prob/len(angles)
    if sum(result.values(),Fraction())!=1:
        raise AssertionError("candidate proposal row is not stochastic")
    return result


@dataclass(frozen=True)
class EdgeReceipt:
    from_path:Path
    to_path:Path
    proposal_forward:Fraction
    proposal_reverse:Fraction
    boltzmann_from:Fraction
    boltzmann_to:Fraction
    transition_forward:Fraction
    transition_reverse:Fraction
    common_flow_mass:Fraction


def propose_edge_receipt(
    sequence:str,q:Fraction,x:Path,y:Path,policy:Policy,
)->EdgeReceipt:
    """Construct an edge's Hastings proof object with only two local states."""
    _check_sequence(sequence)
    _check_q(q)
    _validate_path(x,len(sequence))
    _validate_path(y,len(sequence))
    if len(x)!=len(sequence) or len(y)!=len(sequence) or x==y:
        raise ValueError("receipt needs two different full conformations")
    forward=proposal_row(x,policy).get(y,Fraction())
    reverse=proposal_row(y,policy).get(x,Fraction())
    if forward<=0 or reverse<=0:
        raise ValueError("unsupported directed pivot edge or missing reverse")
    pairs=_possible_hh_pairs(sequence)
    wx=q**_contacts(x,pairs)
    wy=q**_contacts(y,pairs)
    px=min(forward,wy*reverse/wx)
    py=min(reverse,wx*forward/wy)
    if wx*px!=wy*py:
        raise AssertionError("exact Hastings ratio failed detailed balance")
    return EdgeReceipt(x,y,forward,reverse,wx,wy,px,py,wx*px)


def verify_edge_receipt(
    sequence:str,q:Fraction,receipt:EdgeReceipt,policy:Policy,
)->Fraction:
    """Independently recompute the exact two-state balance receipt."""
    if not isinstance(receipt,EdgeReceipt):
        raise ValueError("invalid edge proof object")
    expected=propose_edge_receipt(
        sequence,q,receipt.from_path,receipt.to_path,policy,
    )
    if expected!=receipt:
        raise ValueError("tampered or unstable local Hastings receipt")
    if receipt.common_flow_mass!=receipt.boltzmann_to*receipt.transition_reverse:
        raise ValueError("edge detailed balance violated")
    return receipt.common_flow_mass


def mh_row(sequence:str,q:Fraction,path:Path,policy:Policy)->dict[Path,Fraction]:
    """Locally compute exact corrected transition row, O(polynomial in n)."""
    _check_sequence(sequence)
    _check_q(q)
    _validate_path(path,len(sequence))
    if len(path)!=len(sequence):
        raise ValueError("Hastings requires a full anchored conformation")
    pairs=_possible_hh_pairs(sequence)
    weight=q**_contacts(path,pairs)
    proposal=proposal_row(path,policy)
    transitions:dict[Path,Fraction]={}
    for y,forward in proposal.items():
        if y==path:
            continue
        reverse=proposal_row(y,policy).get(path,Fraction())
        if reverse>0:
            wy=q**_contacts(y,pairs)
            transitions[y]=min(forward,wy*reverse/weight)
    emitted=sum(transitions.values(),Fraction())
    if emitted>1:
        raise AssertionError("local Hastings probabilities exceed 1")
    transitions[path]=1-emitted
    return transitions


@dataclass(frozen=True)
class FiniteAudit:
    states:int
    partition:Fraction
    checked_edges:int
    event_mass:Fraction


def audit_finite_mh(
    sequence:str,q:Fraction,policy:Policy,*,max_states:int=2000,
    threshold:int=1,
)->FiniteAudit:
    """Complete small-model verification, including all directed transitions."""
    if type(threshold) is not int:
        raise ValueError("threshold must be integer")
    universe=all_states(sequence,max_states=max_states)
    lookup={path:i for i,path in enumerate(universe)}
    pairs=_possible_hh_pairs(sequence)
    weights=[q**_contacts(path,pairs) for path in universe]
    rows=[]
    for path in universe:
        row=mh_row(sequence,q,path,policy)
        if not set(row)<=set(lookup) or sum(row.values(),Fraction())!=1:
            raise ValueError("transition escaped the canonical finite state universe")
        rows.append(row)
    checked=0
    for i,path in enumerate(universe):
        for other,p in rows[i].items():
            j=lookup[other]
            if weights[i]*p!=weights[j]*rows[j].get(path,Fraction()):
                raise ValueError("state-dependent policy violated exact detailed balance")
            checked+=1
    return FiniteAudit(
        len(universe),sum(weights,Fraction()),checked,
        sum((w for path,w in zip(universe,weights)
             if _contacts(path,pairs)>=threshold),Fraction()),
    )


def exact_sample(row:dict[Path,Fraction],rng:random.Random)->Path:
    """Draw the next state using integer randomness, never float rounding."""
    denominator=lcm(*(p.denominator for p in row.values()))
    dart=rng.randrange(denominator)
    for path,p in sorted(row.items()):
        dart-=p.numerator*(denominator//p.denominator)
        if dart<0:
            return path
    raise AssertionError("stochastic row exhausted without a sample")


def sample_pivot_chain(
    sequence:str,q:Fraction,policy:Policy,*,steps:int,
    seed:int=0,start:Path|None=None,
)->tuple[Path,...]:
    """Sample locally verified steps without enumerating all conformations.

    The exact *mathematical* stationary guarantee is conditional on the
    proposal policy remaining deterministic and frozen for the whole run.
    The returned steps are not calibrated physical times.
    """
    if type(steps) is not int or steps<0:
        raise ValueError("steps must be nonnegative integer")
    _check_sequence(sequence)
    _check_q(q)
    if start is None:
        start=tuple((i,0) for i in range(len(sequence)))
    _validate_path(start,len(sequence))
    if len(start)!=len(sequence):
        raise ValueError("start must be a full canonical conformation")
    rng=random.Random(seed)
    states=[start]
    current=start
    for _ in range(steps):
        current=exact_sample(mh_row(sequence,q,current,policy),rng)
        states.append(current)
    return tuple(states)


def neural_policy(model,sequence:str,q:Fraction,threshold:int)->Policy:
    """Convert *frozen* JEPA model scores into exact state-dependent weights."""
    from .jepa_kinetic_teacher import conformation_features,pivot_class_features
    classes=move_classes(len(sequence))
    actions=tuple(pivot_class_features(sequence,c) for c in classes)

    def propose(path:Path)->tuple[Fraction,...]:
        state=conformation_features(sequence,q,threshold,path)
        scores=tuple(model.score(state,actions))
        return quantize_ranks(scores)

    return propose
