"""Exact-rational matched ablation for disjoint extremal quadrant shells."""
from __future__ import annotations

import json
from fractions import Fraction
from time import perf_counter
from foldbounds.lattice import ROOT, build_frontier
from foldbounds.lookahead import verify_lookahead
from foldbounds.shell_family import extremal_shell_mass

CASES=(
    ('HPHHPHHPHHHH',Fraction(3,2),2,10,4),
    ('H'*12,Fraction(2),2,10,4),
    ('H'*18,Fraction(2),1,12,6),
)


def comparison():
    for seq,q,t,core,memory in CASES:
        cert=build_frontier(seq,q,t,2)
        options=dict(lookahead=0,use_matching=True,use_transfer=True,
                     use_monotone=True,use_turn_family=True,
                     terminal_box=True,transfer_memory=memory)
        start=perf_counter()
        base=verify_lookahead(cert,**options)
        base_time=perf_counter()-start
        start=perf_counter()
        candidate=verify_lookahead(cert,shell_core=core,**options)
        new_time=perf_counter()-start
        shell=extremal_shell_mass(seq,ROOT,q,t,core_length=core)
        assert candidate.event_lower>=base.event_lower
        assert candidate.complement_lower>=base.complement_lower
        assert candidate.partition_upper==base.partition_upper
        assert candidate.probability_lower>=base.probability_lower
        assert candidate.probability_upper<=base.probability_upper
        assert candidate.event_lower>=shell.event
        assert candidate.complement_lower>=shell.complement
        yield {
            'sequence':seq, 'q':str(q), 'threshold':t,'core_length':core,
            'memory':memory,'seed_count':shell.seed_count,
            'counted_distinct_complete_paths':shell.distinct_paths,
            'old_event_lower':str(base.event_lower),
            'new_event_lower':str(candidate.event_lower),
            'old_complement_lower':str(base.complement_lower),
            'new_complement_lower':str(candidate.complement_lower),
            'old_partition_upper':str(base.partition_upper),
            'new_partition_upper':str(candidate.partition_upper),
            'old_probability_lower':str(base.probability_lower),
            'new_probability_lower':str(candidate.probability_lower),
            'old_probability_upper':str(base.probability_upper),
            'new_probability_upper':str(candidate.probability_upper),
            'old_wall_seconds':round(base_time,5),
            'new_wall_seconds':round(new_time,5),
        }


if __name__=='__main__':
    for entry in comparison():
        print(json.dumps(entry,sort_keys=True))