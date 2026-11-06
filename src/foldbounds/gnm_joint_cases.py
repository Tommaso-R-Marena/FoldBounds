"""Source-only exhaustive case splits for exact joint GNM ranking certificates.

Selected optional contact edges are chosen deterministically using the exact
rank-one trace susceptibility of the *reference-only certain graph*. Every
assignment of the selected edge statuses is covered, and each case's two
endpoint inverse matrices are independently accepted using integer B*N=dI.
A pair's direction is certified globally only when *every* case certifies
that same strict sign. No NMR labels, target structures, or random sampling.

This is a finite Boolean partition and classical graph-update sensitivity,
not a new physical protein theory.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product

from .gnm_joint_rank import (
    Edge, IntegerInverse, RankCounts, canonical_edges, joint_rank_certificate,
    untrusted_propose_inverse, verify_inverse,
)


@dataclass(frozen=True)
class CaseConclusion:
    bits: tuple[int, ...]
    certain_nonlocal: int
    possible_nonlocal: int
    certified_pairs: int


@dataclass(frozen=True)
class GlobalCaseRanks:
    n: int
    selected_edges: tuple[Edge, ...]
    optional_edge_count: int
    baseline_pairs: int
    case_certified_pairs: int
    added_pairs: int
    inverse_witness_count: int
    conclusions: tuple[CaseConclusion, ...]
    signed_orders: tuple[tuple[int, int, int], ...]


def rank_source_only_optional_edges(
    n: int, certain_edges: tuple[Edge, ...], possible_edges: tuple[Edge, ...],
    certain_witness: IntegerInverse,
) -> tuple[tuple[Edge, Fraction], ...]:
    """Sort by the exact single-contact reduction in trace of GNM C.

    For b=e_u-e_v, adding spring bb^T changes C by -vv^T/(1+R),
    v=Cb, R=b^TCb. The trace decrement is ||v||²/(1+R).
    This is an *untrusted heuristic* for which edges to split: correctness
    of the final certificate is independent of its choice.
    """
    certain = canonical_edges(n, certain_edges)
    possible = canonical_edges(n, possible_edges)
    if not set(certain) <= set(possible):
        raise ValueError('certain edges must be a subset of possible edges')
    if not verify_inverse(n, certain, certain_witness):
        raise ValueError('invalid exact certain-graph inverse')
    opt = sorted(set(possible)-set(certain))
    N, d = certain_witness.numerator, certain_witness.denominator
    scored = []
    for u, v in opt:
        dif = [N[i][u]-N[i][v] for i in range(n)]
        resistance_num = n*(N[u][u]+N[v][v]-2*N[u][v])
        assert resistance_num >= 0
        trace_gain = Fraction(n*n*sum(x*x for x in dif),
                              d*(d+resistance_num))
        scored.append(((u,v),trace_gain))
    return tuple(sorted(scored,key=lambda x:(-x[1],x[0])))


def certify_split_cases(
    n: int, certain_edges: tuple[Edge, ...], possible_edges: tuple[Edge, ...],
    *, selected_count: int = 3, selected_override: tuple[Edge, ...] | None = None,
    max_cases: int = 32,
    update_cache=None,
) -> GlobalCaseRanks:
    """Recheck every Boolean world through selected-edge cases.

    Proposals are untrusted. Each endpoint inverse is proved by exact
    integer matrix multiplication before any order certificate is used.
    """
    certain = canonical_edges(n, certain_edges)
    possible = canonical_edges(n, possible_edges)
    opt = set(possible)-set(certain)
    if not set(certain)<=set(possible):
        raise ValueError('invalid contact interval')
    if type(selected_count) is not int or selected_count<0 or selected_count>len(opt):
        raise ValueError('invalid selected split depth')
    if type(max_cases) is not int or max_cases<1:
        raise ValueError('invalid case budget')
    if update_cache is None:
        initial = untrusted_propose_inverse(n,certain)
        baseline_possible = untrusted_propose_inverse(n,possible)
    else:
        if update_cache.n != n or update_cache.base_edges != certain:
            raise ValueError('inverse cache has a different base graph')
        initial = update_cache.get(certain)
        baseline_possible = update_cache.get(possible)
    baseline = joint_rank_certificate(n,certain,possible,initial,baseline_possible)
    if selected_override is None:
        chosen = tuple(edge for edge,_ in rank_source_only_optional_edges(n,certain,possible,initial)[:selected_count])
    else:
        if (not isinstance(selected_override,tuple) or len(selected_override)!=selected_count
                or len(set(selected_override))!=len(selected_override)
                or any(edge not in opt for edge in selected_override)):
            raise ValueError('split override must consist of distinct optional edges')
        chosen = selected_override
    if 2**len(chosen)>max_cases:
        raise ValueError('case enumeration exceeds budget')
    remainder = opt-set(chosen)
    # Cache by exact edge set to avoid redundant integer proposal generation.
    cache: dict[tuple[Edge,...], IntegerInverse] = {certain:initial,possible:baseline_possible}

    def inverse_for(edges):
        key = tuple(sorted(edges))
        if key not in cache:
            cache[key] = (untrusted_propose_inverse(n,key) if update_cache is None
                          else update_cache.get(key))
        return cache[key]

    orders_per_case = []
    conclusions = []
    for bits in product((0,1),repeat=len(chosen)):
        present = {edge for bit,edge in zip(bits,chosen) if bit}
        lower = tuple(sorted(set(certain)|present))
        upper = tuple(sorted(set(lower)|remainder))
        cert = joint_rank_certificate(n, lower,upper,
                                       inverse_for(lower),inverse_for(upper))
        orders_per_case.append({(i,j):sign for i,j,sign in cert.certified_orders})
        conclusions.append(CaseConclusion(
            bits, len(lower)-(n-1),len(upper)-(n-1),cert.joint_certified
        ))
    common = set(orders_per_case[0]) if orders_per_case else set()
    for d in orders_per_case[1:]:
        common &= set(d)
    signed = tuple((i,j,orders_per_case[0][i,j])
                   for i,j in sorted(common)
                   if all(d[i,j]==orders_per_case[0][i,j] for d in orders_per_case))
    # Refinement cannot logically lose previously certified global orders.
    initial_order = {(i,j):v for i,j,v in baseline.certified_orders}
    updated_order = {(i,j):v for i,j,v in signed}
    if any(updated_order.get(pair)!=sign for pair,sign in initial_order.items()):
        raise ArithmeticError('exact case split violates certificate refinement')
    return GlobalCaseRanks(
        n,chosen,len(opt),baseline.joint_certified,len(signed),
        len(signed)-baseline.joint_certified,len(cache),tuple(conclusions),signed
    )
