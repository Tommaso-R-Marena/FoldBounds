"""Ablate guaranteed recent suffix contacts in complement upper bounds.

All comparisons keep identical finite-memory relaxed paths and identical
Boltzmann upper weights. The only difference is whether visible contacts
between recent suffix H residues are used to reject complement membership.

This is an internal soundness/interval diagnostic, NOT a matched runtime
comparison to another transfer implementation or published algorithm.
"""
from __future__ import annotations

import csv
import sys
from fractions import Fraction

from foldbounds.isoperimetric import prefix_contact_cap
from foldbounds.lattice import ROOT, _contacts, _possible_hh_pairs
from foldbounds.monotone import monotone_family_mass
from foldbounds.transfer_bound import (
    future_future_cap, nonbacktracking_mass_upper,
    recent_memory_cross_spectrum,
)


def historical_complement_upper(seq, path, q, threshold, memory):
    """Reconstruct the pre-change complement classifier exactly."""
    fixed = _contacts(path, _possible_hh_pairs(seq))
    flow = future_future_cap(seq, len(path))
    iso = prefix_contact_cap(path, len(seq))
    total = Fraction(0)
    for cross, multiplicity in recent_memory_cross_spectrum(seq, path, memory):
        cap = min(flow, max(0, iso - fixed - cross))
        if fixed + cross < threshold:
            total += multiplicity * q ** (fixed + cross + cap)
    return total


def ratio_lower(event_lower, complement_upper):
    if complement_upper == 0:
        return Fraction(1)
    if event_lower == 0:
        return Fraction(0)
    return event_lower / (event_lower + complement_upper)


def main():
    columns = (
        'sequence', 'n', 'q', 'threshold', 'memory',
        'event_lower', 'old_complement_upper', 'new_complement_upper',
        'relative_complement_reduction', 'old_probability_lower',
        'new_probability_lower', 'partition_upper_unchanged',
    )
    writer = csv.DictWriter(sys.stdout, fieldnames=columns)
    writer.writeheader()
    cases = (
        ('H' * 8, Fraction(2), 1, 4),
        ('H' * 12, Fraction(2), 1, 4),
        ('H' * 16, Fraction(2), 1, 4),
        ('H' * 18, Fraction(2), 1, 6),
        ('HPHHPHHPHHHH', Fraction(3, 2), 2, 4),
        ('HPHHPHHPHHHH', Fraction(3, 2), 2, 6),
    )
    for seq, q, threshold, memory in cases:
        old_comp = historical_complement_upper(seq, ROOT, q, threshold, memory)
        new = nonbacktracking_mass_upper(seq, ROOT, q, threshold,
                                        suffix_memory=memory)
        lower = monotone_family_mass(seq, ROOT, q, threshold)
        assert new.complement <= old_comp
        old_p = ratio_lower(lower.event, old_comp)
        new_p = ratio_lower(lower.event, new.complement)
        assert new_p >= old_p
        writer.writerow(dict(
            sequence=seq, n=len(seq), q=q, threshold=threshold, memory=memory,
            event_lower=lower.event,
            old_complement_upper=old_comp,
            new_complement_upper=new.complement,
            relative_complement_reduction=(
                str((old_comp - new.complement) / old_comp)
                if old_comp > 0 else 'N/A'),
            old_probability_lower=old_p,
            new_probability_lower=new_p,
            partition_upper_unchanged='yes (same weights and relaxed paths)',
        ))


if __name__ == '__main__':
    main()
