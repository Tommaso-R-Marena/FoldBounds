"""JEPA student of exact toy Markov dynamics and fail-closed proposal gating.

The student NEVER supplies certified rates, flux or Boltzmann probabilities:
only an independent exact Markov checker can authorize a proposed policy.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import random

from .jepa_guidance import STATE_DIM, ACTION_DIM, build_model
from .jepa_reversible_kinetics import (
    Kernel, MoveClass, _angles, build_kernel, event_flux, optimize_exact_flux,
    pivot, uniform_classes, verify_kernel,
)
from .lattice import Path, _contacts, _possible_hh_pairs


def conformation_features(sequence: str, q: Fraction, threshold: int,
                          path: Path) -> tuple[float, ...]:
    """Cheap 10-dimensional geometric + contact descriptors; no oracle labels."""
    n = len(sequence)
    pairs = _possible_hh_pairs(sequence)
    xs,ys = zip(*path)
    cx,cy = sum(xs)/n,sum(ys)/n
    rg2 = sum((x-cx)**2+(y-cy)**2 for x,y in path)/n
    value = (
        n/32,sequence.count("H")/n,threshold/max(1,n),
        math.log1p(float(q-1)),_contacts(path,pairs)/n,rg2/(n*n),
        (max(xs)-min(xs))/n,(max(ys)-min(ys))/n,
        math.hypot(*path[-1])/n,len(pairs)/max(1,n*n),
    )
    assert len(value)==STATE_DIM
    return value


def pivot_class_features(sequence: str, action: MoveClass) -> tuple[float, ...]:
    n=len(sequence)
    k,family=action
    prefix,suffix=sequence[:k+1],sequence[k+1:]
    h=lambda c: float(c=="H")
    result=(
        k/n,float(family=="quarter"),h(sequence[k]),h(sequence[k+1]),
        prefix.count("H")/len(prefix),suffix.count("H")/len(suffix),
        (n-k-1)/n,h(sequence[k-1]),h(sequence[k+2]) if k+2<n else 0.0,
    )
    assert len(result)==ACTION_DIM
    return result


def conditional_pivot_outcomes(kernel: Kernel, index: int, action: MoveClass):
    """Exact conditional transition P(y|x,action), including rejected moves."""
    if type(index) is not int or not 0<=index<len(kernel.states):
        raise ValueError("invalid state index")
    if action not in kernel.classes:
        raise ValueError("invalid pivot class")
    path=kernel.states[index]
    lookup={path:i for i,path in enumerate(kernel.states)}
    row={index:Fraction(1)}
    angles=_angles(action[1])
    for angle in angles:
        candidate=pivot(path,action[0],angle)
        if candidate is None or candidate==path:
            continue
        j=lookup[candidate]
        accepted=min(Fraction(1),kernel.weights[j]/kernel.weights[index])/len(angles)
        row[j]=row.get(j,Fraction())+accepted
        row[index]-=accepted
    result=tuple(sorted((j,p) for j,p in row.items() if p))
    if sum((p for _,p in result),Fraction())!=1:
        raise AssertionError("conditional transition not stochastic")
    return result


@dataclass(frozen=True)
class PivotExample:
    state: tuple[float,...]
    action: tuple[float,...]
    successors: tuple[tuple[tuple[float,...],Fraction],...]
    event_escape: Fraction
    expected_score_delta: Fraction


def collect_pivot_examples(kernel: Kernel, threshold: int, *,
                           max_states: int=32, seed: int=0) -> tuple[PivotExample,...]:
    """Generate training examples from frozen *training sequences only*."""
    if type(max_states) is not int or max_states<1:
        raise ValueError("max_states must be a positive integer")
    rng=random.Random(seed)
    selected=list(range(len(kernel.states)))
    if len(selected)>max_states:
        selected=rng.sample(selected,max_states)
    features=[conformation_features(kernel.sequence,kernel.q,threshold,path)
              for path in kernel.states]
    result=[]
    for i in selected:
        for action in kernel.classes:
            outcomes=conditional_pivot_outcomes(kernel,i,action)
            crossing=sum((p for j,p in outcomes
                          if kernel.scores[i]>=threshold and kernel.scores[j]<threshold),Fraction())
            delta=sum((p*(kernel.scores[j]-kernel.scores[i]) for j,p in outcomes),Fraction())
            result.append(PivotExample(
                features[i],pivot_class_features(kernel.sequence,action),
                tuple((features[j],p) for j,p in outcomes),
                crossing,delta/len(kernel.sequence),
            ))
    return tuple(result)


def fit_pivot_jepa(examples: tuple[PivotExample,...], *, epochs: int=30,
                   lr: float=1e-3, seed: int=0, model=None):
    """Train EMA-target JEPA on the exact *mixture of future embeddings*.

    Since the Markov process is stochastic, embedding a mean future vector
    is not the same as averaging future embeddings. We use the latter.
    The actor learns conditional A->not-A escape probabilities; the critic
    learns signed expected contact changes. Both are fallible regressions.
    """
    import torch
    from torch import nn
    if not examples or type(epochs) is not int or epochs<1 or lr<=0:
        raise ValueError("nonempty examples and positive epochs/lr required")
    torch.manual_seed(seed)
    model=build_model() if model is None else model
    optimizer=torch.optim.AdamW((p for p in model.parameters() if p.requires_grad),lr=lr)
    states=torch.tensor([x.state for x in examples],dtype=torch.float32)
    actions=torch.tensor([x.action for x in examples],dtype=torch.float32)
    escapes=torch.tensor([float(x.event_escape) for x in examples],dtype=torch.float32)
    deltas=torch.tensor([float(x.expected_score_delta) for x in examples],dtype=torch.float32)
    futures=[];owners=[];weights=[]
    for i,example in enumerate(examples):
        for vector,prob in example.successors:
            futures.append(vector)
            owners.append(i)
            weights.append(float(prob))
    future_tensor=torch.tensor(futures,dtype=torch.float32)
    owner_tensor=torch.tensor(owners,dtype=torch.int64)
    weight_tensor=torch.tensor(weights,dtype=torch.float32).unsqueeze(1)
    future_mean=torch.zeros_like(states).index_add(0,owner_tensor,
                                                   future_tensor*weight_tensor)
    losses=[]
    for _ in range(epochs):
        model.train()
        actor,critic,pred=model.heads(states,actions)
        with torch.no_grad():
            embedded=model.target_encoder(future_tensor)
            future_z=torch.zeros_like(pred).index_add(
                0,owner_tensor,embedded*weight_tensor,
            )
        loss=(
            nn.functional.smooth_l1_loss(pred,future_z)
            + nn.functional.mse_loss(model.decoder(pred),future_mean)
            + nn.functional.mse_loss(torch.sigmoid(actor),escapes)
            + nn.functional.mse_loss(critic,deltas)
        )
        if not torch.isfinite(loss):
            raise ArithmeticError("nonfinite JEPA kinetic loss")
        optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_((p for p in model.parameters() if p.requires_grad),2.0)
        optimizer.step()
        model.update_target()
        losses.append(float(loss.detach()))
    return model,tuple(losses)


def rank_classes_from_model(model, kernel: Kernel, threshold: int) -> tuple[float,...]:
    """Untrusted stationary-weighted class scores on a fully enumerated world.

    No transition-label lookups are used for class ranking, but this pilot
    still requires enumeration and is NOT a scalable protein folding method.
    """
    import torch
    verify_kernel(kernel)
    indices=[i for i,k in enumerate(kernel.scores) if k>=threshold]
    if not indices:
        return (0.0,)*len(kernel.classes)
    states=torch.tensor([
        conformation_features(kernel.sequence,kernel.q,threshold,kernel.states[i])
        for i in indices
    ],dtype=torch.float32)
    stationary=torch.tensor([
        float(kernel.weights[i]/kernel.partition) for i in indices
    ],dtype=torch.float32)
    model.eval()
    values=[]
    with torch.no_grad():
        for action in kernel.classes:
            features=torch.tensor([pivot_class_features(kernel.sequence,action)]*
                                  len(indices),dtype=torch.float32)
            predictions,_,_=model.heads(states,features)
            values.append(float((torch.sigmoid(predictions)*stationary).sum()))
    return tuple(values)


@dataclass(frozen=True)
class GatedProposal:
    accepted: bool
    chosen_class: MoveClass
    baseline_flux: Fraction
    achieved_flux: Fraction
    oracle_flux: Fraction
    actual_probs: tuple[Fraction,...]
    raw_candidate_flux: Fraction


def gate_model_proposal(kernel: Kernel, threshold: int,
                        scores: tuple[float,...], *,
                        floor: Fraction | None=None) -> GatedProposal:
    """Fail-closed: a learned proposal is used only if exact flux is no worse.

    This gate itself needs exhaustive exact enumeration; it is not an
    efficient large-protein certification procedure.
    """
    verify_kernel(kernel)
    if len(scores)!=len(kernel.classes) or any(
        isinstance(score,bool) or not isinstance(score,(int,float))
        or not math.isfinite(score) for score in scores
    ):
        raise ValueError("invalid untrusted JEPA score vector")
    baseline=uniform_classes(len(kernel.sequence))
    if floor is None:
        floor=min(baseline)/4
    if (isinstance(floor,bool) or not isinstance(floor,Fraction)
            or floor<=0 or floor>min(baseline)):
        raise ValueError("invalid rational probability floor")
    best=max(range(len(scores)),key=lambda i:(scores[i],-i))
    proposed=[floor]*len(scores)
    proposed[best]+=1-floor*len(scores)
    candidate=build_kernel(kernel.sequence,kernel.q,class_probs=tuple(proposed),
                           max_states=max(1,len(kernel.states)))
    reference=build_kernel(kernel.sequence,kernel.q,class_probs=baseline,
                           max_states=max(1,len(kernel.states)))
    verify_kernel(candidate)
    verify_kernel(reference)
    baseline_flux=event_flux(reference,threshold)[1]
    raw_flux=event_flux(candidate,threshold)[1]
    accepted=raw_flux>=baseline_flux
    chosen=candidate if accepted else reference
    oracle=optimize_exact_flux(reference,threshold).optimized_flux
    actual=event_flux(chosen,threshold)[1]
    if not baseline_flux<=actual<=oracle:
        raise AssertionError("checked kinetic proposal violates performance gate")
    return GatedProposal(accepted,kernel.classes[best],baseline_flux,
                         actual,oracle,chosen.class_probs,raw_flux)
