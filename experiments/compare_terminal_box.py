"""Matched geometry-on/off endpoint transfer bound ablation.

The optional terminal-conditioned isoperimetric cap is the only changed
mathematical bound. Both cases count the same finite-memory relaxed paths;
retain all other upper/lower bound settings. Single-run times are noisy.
"""
from __future__ import annotations

import csv
import sys
from fractions import Fraction
from time import perf_counter

from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead


def run():
    fields = (
        'sequence', 'n', 'q', 'event_at_least', 'memory',
        'box_enabled', 'elapsed_ms', 'partition_lower', 'partition_upper',
        'event_lower', 'event_upper', 'complement_lower', 'complement_upper',
        'probability_lower', 'probability_upper', 'probability_width',
        'exact_z_if_known', 'exact_probability_if_known',
    )
    out = csv.DictWriter(sys.stdout, fieldnames=fields)
    out.writeheader()
    cases = (
        ('H' * 8, Fraction(2), 1, 2),
        ('H' * 12, Fraction(2), 1, 4),
        ('H' * 16, Fraction(2), 1, 4),
        ('H' * 18, Fraction(2), 1, 6),
        ('HPHHPHHPHHHH', Fraction(3, 2), 2, 4),
        ('HHPHHHPHHH', Fraction(2), 2, 4),
    )
    for sequence, q, threshold, memory in cases:
        cert = build_frontier(sequence, q, threshold, 2)
        exact_z = exact_p = None
        if len(sequence) <= 8:
            spectrum = exact_spectrum(sequence)
            exact_z = sum((multiplicity * q ** score
                           for score, multiplicity in spectrum.items()), Fraction(0))
            event = sum((multiplicity * q ** score for score, multiplicity
                         in spectrum.items() if score >= threshold), Fraction(0))
            exact_p = event / exact_z
        previous = None
        for enabled in (False, True):
            start = perf_counter()
            b = verify_lookahead(
                cert, lookahead=0, use_matching=True, use_transfer=True,
                transfer_memory=memory, use_monotone=True, use_turn_family=True,
                terminal_box=enabled)
            elapsed_ms = 1000 * (perf_counter() - start)
            if previous:
                assert b.partition_lower == previous.partition_lower
                assert b.partition_upper <= previous.partition_upper
                assert b.event_upper <= previous.event_upper
                assert b.complement_upper <= previous.complement_upper
                assert b.probability_lower >= previous.probability_lower
                assert b.probability_upper <= previous.probability_upper
            if exact_z is not None:
                assert b.partition_lower <= exact_z <= b.partition_upper
                assert b.probability_lower <= exact_p <= b.probability_upper
            previous = b
            out.writerow(dict(
                sequence=sequence, n=len(sequence), q=q,
                event_at_least=threshold, memory=memory,
                box_enabled=enabled, elapsed_ms=f'{elapsed_ms:.3f}',
                partition_lower=b.partition_lower, partition_upper=b.partition_upper,
                event_lower=b.event_lower, event_upper=b.event_upper,
                complement_lower=b.complement_lower,
                complement_upper=b.complement_upper,
                probability_lower=b.probability_lower,
                probability_upper=b.probability_upper,
                probability_width=b.probability_upper-b.probability_lower,
                exact_z_if_known=exact_z if exact_z is not None else '',
                exact_probability_if_known=exact_p if exact_p is not None else '',
            ))


if __name__ == '__main__':
    run()
