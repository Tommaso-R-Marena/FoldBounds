"""Exact matched-leaf-budget controlled ablation for optimal shell portfolios.

No uncontrolled merging of certificates and no claims of matching proposer
computational budget: the DP explores the full candidate prefix tree to depth.
"""
from __future__ import annotations

import json
from fractions import Fraction
from time import perf_counter

from foldbounds.lattice import ROOT
from foldbounds.shell_family import extremal_shell_mass
from foldbounds.shell_portfolio import optimize_shell_portfolio, verify_shell_portfolio

CASES=(
    ('HPHHPHHPHHHH',Fraction(3,2),2),
    ('H'*12,Fraction(2),2),
    ('H'*16,Fraction(2),2),
)

def comparisons():
    for seq,q,threshold in CASES:
        for core,budget in ((5,25),(6,50)):
            start=perf_counter()
            uniform=extremal_shell_mass(seq,ROOT,q,threshold,core_length=core)
            uniform_seconds=perf_counter()-start
            start=perf_counter()
            plan=optimize_shell_portfolio(seq,ROOT,q,threshold,
                   leaf_budget=budget,max_core_length=9,objective='event')
            planner_seconds=perf_counter()-start
            start=perf_counter()
            checked=verify_shell_portfolio(plan.certificate)
            checker_seconds=perf_counter()-start
            assert checked==plan.mass
            assert checked.leaf_count<=budget
            if uniform.seed_count <= budget:
                assert checked.event >= uniform.event
            yield {
                'sequence':seq,'q':str(q),'threshold':threshold,
                'uniform_core':core,'uniform_leaves':uniform.seed_count,
                'uniform_event_lower':str(uniform.event),
                'uniform_complement_lower':str(uniform.complement),
                'portfolio_leaf_budget':budget,
                'portfolio_actual_leaves':checked.leaf_count,
                'portfolio_event_lower':str(checked.event),
                'portfolio_complement_lower':str(checked.complement),
                'candidate_prefixes_examined':plan.candidate_nodes,
                'planner_wall_s':round(planner_seconds,5),
                'checker_wall_s':round(checker_seconds,5),
                'uniform_wall_s':round(uniform_seconds,5),
            }

if __name__=='__main__':
    for entry in comparisons():
        print(json.dumps(entry,sort_keys=True))