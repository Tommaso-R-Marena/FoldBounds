"""Controlled ablation of 0..2-turn constructive lower masses.

Upper-bound settings are identical across each pair, and the optional
turn family is combined with monotone/ray witnesses via max, never sum.
No published baseline or real-protein accuracy claim is made.
"""
from __future__ import annotations

import csv
import sys
from fractions import Fraction
from time import perf_counter

from foldbounds.lattice import ROOT, build_frontier
from foldbounds.lookahead import verify_lookahead
from foldbounds.turn_family import turn_family_mass


def main():
    fields = ('sequence', 'n', 'q', 'threshold', 'memory', 'turn_paths',
              'turn_event_lower', 'turn_complement_lower',
              'baseline_event_lower', 'enhanced_event_lower',
              'baseline_probability_lower', 'enhanced_probability_lower',
              'baseline_partition_upper', 'enhanced_partition_upper',
              'baseline_ms', 'enhanced_ms')
    writer = csv.DictWriter(sys.stdout, fieldnames=fields)
    writer.writeheader()
    cases = (
        ('H'*8, Fraction(2), 1, 2),
        ('H'*12, Fraction(2), 1, 2),
        ('H'*16, Fraction(2), 1, 2),
        ('HPHHPHHPHHHH', Fraction(3, 2), 2, 4),
        ('HHPHHHPHHH', Fraction(2), 2, 4),
    )
    for sequence, q, threshold, memory in cases:
        cert = build_frontier(sequence, q, threshold, 2)
        family = turn_family_mass(sequence, ROOT, q, threshold)
        kw = dict(lookahead=0, use_matching=True, use_transfer=True,
                  use_monotone=True, transfer_memory=memory)
        start = perf_counter()
        old = verify_lookahead(cert, use_turn_family=False, **kw)
        old_ms = (perf_counter() - start) * 1000
        start = perf_counter()
        new = verify_lookahead(cert, use_turn_family=True, **kw)
        new_ms = (perf_counter() - start) * 1000
        assert new.event_lower >= old.event_lower
        assert new.complement_lower >= old.complement_lower
        assert new.partition_upper == old.partition_upper
        assert new.event_upper == old.event_upper
        assert new.complement_upper == old.complement_upper
        assert new.probability_lower >= old.probability_lower
        assert new.probability_upper <= old.probability_upper
        writer.writerow(dict(
            sequence=sequence, n=len(sequence), q=q, threshold=threshold,
            memory=memory, turn_paths=family.distinct_paths,
            turn_event_lower=family.event,
            turn_complement_lower=family.complement,
            baseline_event_lower=old.event_lower,
            enhanced_event_lower=new.event_lower,
            baseline_probability_lower=old.probability_lower,
            enhanced_probability_lower=new.probability_lower,
            baseline_partition_upper=old.partition_upper,
            enhanced_partition_upper=new.partition_upper,
            baseline_ms=f'{old_ms:.3f}', enhanced_ms=f'{new_ms:.3f}'))


if __name__ == '__main__':
    main()
