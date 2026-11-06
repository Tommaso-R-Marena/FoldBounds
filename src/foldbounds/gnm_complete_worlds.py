"""Exact complete classification of rank signs over finite optional-contact graphs.

The tree, claims and counterexample masks are UNTRUSTED proposals. This
checker validates every branch of the tree against exact Loewner/rank-two
bounds (with verified inverse matrices) and validates both opposite signs
for every purported non-invariant pair in actual complete contact graphs.

The exhaustive domain is ALL 2^m independently selected optional edges;
it is a SUPerset of networks achievable by continuous perturbed coordinates.
NO claim of physical 3D realizability of arbitrary Boolean edge masks.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .gnm_joint_rank import Edge, canonical_edges
from .gnm_reduced_woodbury import ReducedWorldCache,exact_joint_sign


@dataclass(frozen=True)
class CompleteWorldAudit:
    residues: int
    optional_contacts: int
    all_worlds: int
    invariant_pairs: int
    reversing_pairs: int
    original_certified: int
    extra_certified: int
    tree_nodes: int
    tree_leaves: int
    tree_max_depth: int
    counterexample_networks: int
    exact_inverse_graphs: int


def _strict_rank_sign(denlo: int, denhi: int, nlo, nhi,
                      i: int, j: int) -> int:
    """Check exact rank-two Loewner envelope for the strict sign of i-j.

    Returns 0 for abstention, +1 for Cii>Cjj, -1 for Cii<Cjj.
    Both input inverses must already be verified against graph endpoints.
    """
    a=nlo[i][i]*denhi-nhi[i][i]*denlo
    b=nlo[j][j]*denhi-nhi[j][j]*denlo
    c=nlo[i][j]*denhi-nhi[i][j]*denlo
    if a<0 or b<0 or a*b<c*c:
        raise ArithmeticError('endpoint inverse difference violates pairwise PSD')
    f=(nhi[i][i]-nhi[j][j])*denlo
    s=2*f+a-b
    rad=(a+b)**2-4*c*c
    if rad<0:raise ArithmeticError('negative discriminant')
    return (1 if s>0 and s*s>rad else
            -1 if s<0 and s*s>rad else 0)


def _parse_signed_orders(items: Any, n:int, label:str):
    if not isinstance(items,list):
        raise ValueError(label+' must be a JSON list')
    result={}
    for triple in items:
        if (not isinstance(triple,list) or len(triple)!=3
            or any(type(x) is not int for x in triple)):
            raise ValueError(label+' has a malformed signed pair')
        i,j,s=triple
        if not(0<=i<j<n) or s not in (-1,1) or (i,j) in result:
            raise ValueError(label+' contains duplicate/invalid signed pairs')
        result[(i,j)]=s
    return result


def verify_complete_contact_worlds(
    n:int, certain_edges:tuple[Edge,...], possible_edges:tuple[Edge,...],
    certificate:dict, *, selected_prior:tuple[Edge,...]=(),
    max_tree_nodes:int=10000, max_witness_masks:int=512,
    max_optional:int=18,
) -> CompleteWorldAudit:
    """Check all signs over EVERY Boolean optional-contact graph exactly.

    An untrusted numerical or symbolic producer may propose the tree and
    network masks. No floating-point calculation is used for acceptance.

    The initial selected_prior case split is itself recertified using
    the independent integer witnesses from PR #79. The tree then proves
    additional invariant signs on EVERY partition cell. For each
    remaining pair, two exact full-network witnesses show opposite signs.
    Validates a *complete* partition of n(n-1)/2 pair comparisons.
    """
    if type(n) is not int or not 2<=n<=128:
        raise ValueError('n out of bounds')
    certain=canonical_edges(n,certain_edges)
    possible=canonical_edges(n,possible_edges)
    if not set(certain)<=set(possible):
        raise ValueError('certain edges not contained in possible edges')
    optional=tuple(sorted(set(possible)-set(certain)))
    m=len(optional)
    if (type(max_optional) is not int or max_optional<0 or m>max_optional
        or type(max_tree_nodes) is not int or max_tree_nodes<1
        or type(max_witness_masks) is not int or max_witness_masks<1):
        raise ValueError('resource limit exceeded')
    if not isinstance(certificate,dict) or certificate.get('schema')!='foldbounds-gnm-full-contact-worlds-v1':
        raise ValueError('wrong certificate format')
    if certificate.get('optional_edges')!=[list(e) for e in optional]:
        raise ValueError('optional contact order / geometry mismatch')
    if (not isinstance(selected_prior,tuple)
            or len(set(selected_prior))!=len(selected_prior)
            or any(edge not in optional for edge in selected_prior)):
        raise ValueError('prior selected edges must be distinct optional contacts')
    # Exact verified base integer inverse, then exact verified k×k Woodbury
    # reductions, with k<=m. No numerical matrix inverses or stochastic grid.
    reduced=ReducedWorldCache(n,certain,possible,max_optional=max_optional)
    prior_indices=tuple(optional.index(edge) for edge in selected_prior)
    selected_mask=sum(1<<k for k in prior_indices)
    all_mask=(1<<m)-1
    previous=None
    for branch in range(1<<len(prior_indices)):
        lo_mask=sum((1<<k) for bit,k in enumerate(prior_indices)
                    if branch&(1<<bit))
        hi_mask=all_mask ^ (selected_mask ^ lo_mask)
        low=reduced.get(lo_mask);high=reduced.get(hi_mask)
        orders={}
        for i in range(n):
            for j in range(i+1,n):
                sign=exact_joint_sign(low,high,i,j)
                if sign:
                    orders[(i,j)]=sign
        if previous is None:
            previous=orders
        else:
            previous={pair:sign for pair,sign in previous.items()
                      if orders.get(pair)==sign}
    if previous is None:
        raise ArithmeticError('no prior graph cases were checked')
    extra=_parse_signed_orders(certificate.get('extra_stable_orders'),n,'extra_stable_orders')
    if set(previous)&set(extra):
        raise ValueError('additional signs overlap prior certificates')
    flips=certificate.get('changing_pairs')
    if not isinstance(flips,list):raise ValueError('changing_pairs must be a list')
    flip_set=set()
    for item in flips:
        if (not isinstance(item,list) or len(item)!=2
            or any(type(x) is not int for x in item)):
            raise ValueError('invalid changing pair')
        i,j=item
        if not 0<=i<j<n or (i,j) in flip_set:raise ValueError('duplicate/invalid changing pair')
        flip_set.add((i,j))
    if (set(previous)|set(extra)) & flip_set:
        raise ValueError('rank cannot be invariant and reversing')
    if len(previous)+len(extra)+len(flip_set)!=n*(n-1)//2:
        raise ValueError('incomplete classification: some pairs missing')

    tree=certificate.get('tree')
    nodes=0;leaves=0;max_depth=0

    def prove(subtree,yes:frozenset[Edge],no:frozenset[Edge],available:frozenset[Edge],depth:int):
        nonlocal nodes,leaves,max_depth
        nodes+=1;max_depth=max(max_depth,depth)
        if nodes>max_tree_nodes or depth>m:raise ValueError('tree budget exceeded')
        if not isinstance(subtree,dict):raise ValueError('tree node is not an object')
        if set(subtree)=={'leaf'} and subtree['leaf'] is True:
            leaves+=1
            low_mask=sum(1<<optional.index(e) for e in yes)
            high_mask=all_mask ^ sum(1<<optional.index(e) for e in no)
            low=reduced.get(low_mask);high=reduced.get(high_mask)
            for (i,j),sign in extra.items():
                actual=exact_joint_sign(low,high,i,j)
                if actual!=sign:
                    raise ArithmeticError(
                        f'unjustified cut at leaf for {(i,j)}: {actual} != {sign}')
            return
        if set(subtree)!={'split','zero','one'} or type(subtree['split']) is not int:
            raise ValueError('malformed split: both children required')
        index=subtree['split']
        if not 0<=index<m or optional[index] not in available:
            raise ValueError('split edge is invalid or already assigned')
        e=optional[index];rest=available-{e}
        # Both children are *always* checked, so there are no uncovered worlds.
        prove(subtree['zero'],yes,no|{e},rest,depth+1)
        prove(subtree['one'],yes|{e},no,rest,depth+1)

    prove(tree,frozenset(),frozenset(),frozenset(optional),0)

    mask_list=certificate.get('witness_masks')
    if (not isinstance(mask_list,list) or len(mask_list)>max_witness_masks
        or (not mask_list and flip_set)):
        raise ValueError('counterexample masks missing / over budget')
    if any(type(mask) is not int or not 0<=mask<(1<<m) for mask in mask_list):
        raise ValueError('counterexample mask out of range')
    if len(set(mask_list))!=len(mask_list):
        raise ValueError('duplicate counterexample masks')
    needed={pair:{-1,1} for pair in flip_set}
    for mask in mask_list:
        world=reduced.get(mask)
        diag=world.diagonal_numerators()
        for i,j in flip_set:
            diff=diag[i]-diag[j]
            sign=(diff>0)-(diff<0)
            needed[i,j].discard(sign)
    if any(signs for signs in needed.values()):
        remaining=[(pair,sorted(signs)) for pair,signs in needed.items() if signs]
        raise ArithmeticError(f'counterexample masks do not prove both signs: {remaining[:5]}')
    return CompleteWorldAudit(n,m,1<<m,len(previous)+len(extra),len(flip_set),
                              len(previous),len(extra),nodes,leaves,max_depth,
                              len(mask_list),reduced.count)
