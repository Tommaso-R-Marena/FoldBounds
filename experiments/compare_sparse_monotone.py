"""Controlled leaf-budget comparisons and exact independent probability enclosure checks."""
from __future__ import annotations
from fractions import Fraction
from time import perf_counter
import json
from foldbounds.lattice import ROOT, build_frontier, exact_spectrum
from foldbounds.shell_portfolio import optimize_shell_portfolio
from foldbounds.sparse_monotone import optimize_sparse_monotone,verify_sparse_monotone
from foldbounds.lookahead import verify_lookahead

def benchmark(seq,q,t,budget=25,maxdepth=9):
    old=optimize_shell_portfolio(seq,ROOT,q,t,leaf_budget=budget,
                                  max_core_length=maxdepth,objective='event')
    t0=perf_counter()
    new=optimize_sparse_monotone(seq,ROOT,q,t,leaf_budget=budget,
                                 max_core_length=maxdepth,objective='event')
    proposer_wall=perf_counter()-t0
    t1=perf_counter()
    mass=verify_sparse_monotone(new.certificate)
    checker_wall=perf_counter()-t1
    assert mass==new.mass
    assert mass.event>=old.mass.event
    t_union=perf_counter()
    union=optimize_sparse_monotone(seq,ROOT,q,t,leaf_budget=budget,
                                   max_core_length=maxdepth,
                                   objective='event',include_nonmonotone_turns=True)
    union_planner_wall=perf_counter()-t_union
    t_checker=perf_counter()
    checked_union=verify_sparse_monotone(union.certificate)
    union_checker_wall=perf_counter()-t_checker
    assert checked_union==union.mass
    assert union.mass.event>=mass.event
    case=dict(sequence=seq,q=str(q),threshold=t,leaf_budget=budget,
              depth_cap=maxdepth,candidates=new.candidate_prefixes,
              old_event_lower=str(old.mass.event),
              sparse_event_lower=str(new.mass.event),
              turn_union_event_lower=str(union.mass.event),
              turn_union_complement_lower=str(union.mass.complement),
              turn_union_walks=union.mass.distinct_paths,
              turn_union_candidate_prefixes=union.candidate_prefixes,
              turn_union_planner_seconds=round(union_planner_wall,5),
              turn_union_checker_seconds=round(union_checker_wall,5),
              sparse_complement_lower=str(new.mass.complement),
              sparse_walks=mass.distinct_paths,
              sparse_leaves=mass.leaf_count,
              planner_seconds=round(proposer_wall,5),
              independent_checker_seconds=round(checker_wall,5))
    cert=build_frontier(seq,q,t,2)
    kw=dict(lookahead=0,use_matching=True,use_transfer=True,
            transfer_memory=4,use_monotone=True,use_turn_family=True,
            use_coupled_polytope=True)
    baseline=verify_lookahead(cert,**kw,shell_portfolios=(old.certificate,))
    improved=verify_lookahead(cert,**kw,shell_portfolios=(old.certificate,),
                               sparse_certificates=(new.certificate,))
    improved_union=verify_lookahead(cert,**kw,shell_portfolios=(old.certificate,),
                                     sparse_certificates=(new.certificate,union.certificate))
    assert improved.event_lower>=baseline.event_lower
    assert improved.complement_lower>=baseline.complement_lower
    assert improved.partition_upper==baseline.partition_upper
    assert improved.probability_lower>=baseline.probability_lower
    assert improved.probability_upper<=baseline.probability_upper
    assert improved_union.probability_lower>=improved.probability_lower
    assert improved_union.probability_upper<=improved.probability_upper
    case.update(previous_probability_lower=str(baseline.probability_lower),
                sparse_probability_lower=str(improved.probability_lower),
                previous_probability_upper=str(baseline.probability_upper),
                sparse_probability_upper=str(improved.probability_upper),
                union_probability_lower=str(improved_union.probability_lower),
                union_probability_upper=str(improved_union.probability_upper))
    if len(seq)<=12:
        spectrum=exact_spectrum(seq)
        z=sum((count*q**score for score,count in spectrum.items()),Fraction(0))
        e=sum((count*q**score for score,count in spectrum.items() if score>=t),Fraction(0))
        assert improved.probability_lower<=e/z<=improved.probability_upper
        case['exact_probability']=str(e/z)
    print(json.dumps(case,sort_keys=True),flush=True)

if __name__=='__main__':
    benchmark('HPHHPHHPHHHH',Fraction(3,2),2)
    benchmark('H'*12,Fraction(2),2)
    benchmark('H'*16,Fraction(2),2)
    # This is NOT a full model-mass calculation: the restricted family is
    # represented by only 25 prefix leaves and integer DP, with no n=64 oracle.
    t0=perf_counter()
    large=optimize_sparse_monotone('H'*64,ROOT,Fraction(2),1,
                                   leaf_budget=25,max_core_length=6)
    checked=verify_sparse_monotone(large.certificate)
    assert checked==large.mass
    assert checked.distinct_paths==20032011142543966146
    print(json.dumps(dict(sequence='H'*64,q='2',threshold=1,
                          type='scalability_restricted_family_NOT_exact_partition',
                          selected_walks=checked.distinct_paths,
                          event_lower=str(checked.event),
                          complement_lower=str(checked.complement),
                          selected_leaves=checked.leaf_count,
                          prefixes_examined=large.candidate_prefixes,
                          total_seconds=round(perf_counter()-t0,5))),flush=True)

    # Controlled comparison against the strongest preexisting prefix endpoint
    # upper bounds and fixed depth-10 shell lower masses. This uses the SAME
    # upper-bound settings for the off/on probability enclosures.
    seq='HPHHPHHPHHHH'
    q=Fraction(3,2)
    t=2
    candidates=tuple(optimize_sparse_monotone(
        seq,ROOT,q,t,leaf_budget=200,max_core_length=8,objective=objective,
        include_nonmonotone_turns=True).certificate
        for objective in ('event','complement'))
    assert candidates[0].leaves and candidates[1].leaves
    m_event=verify_sparse_monotone(candidates[0])
    m_complement=verify_sparse_monotone(candidates[1])
    assert m_event.event==Fraction(421407,32)
    assert m_complement.complement==Fraction(19426)
    cert=build_frontier(seq,q,t,2)
    controls=dict(lookahead=0,use_matching=True,use_transfer=True,
                  transfer_memory=4,terminal_box=True,use_monotone=True,
                  use_turn_family=True,shell_core=10,use_coupled_polytope=True)
    old=verify_lookahead(cert,**controls)
    new=verify_lookahead(cert,**controls,sparse_certificates=candidates)
    assert new.probability_lower>old.probability_lower
    assert new.probability_upper<=old.probability_upper
    assert new.partition_upper==old.partition_upper
    exact=exact_spectrum(seq)
    z=sum((count*q**k for k,count in exact.items()),Fraction(0))
    e=sum((count*q**k for k,count in exact.items() if k>=t),Fraction(0))
    assert new.probability_lower<=e/z<=new.probability_upper
    print(json.dumps(dict(type='200_leaf_vs_strongest_terminal_box',
                          sequence=seq,q=str(q),threshold=t,
                          event_lower=str(m_event.event),
                          complement_lower=str(m_complement.complement),
                          old_probability_lower=str(old.probability_lower),
                          new_probability_lower=str(new.probability_lower),
                          old_probability_upper=str(old.probability_upper),
                          new_probability_upper=str(new.probability_upper),
                          exact_probability=str(e/z))),flush=True)
