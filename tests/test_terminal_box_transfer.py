"""Independent test oracles for terminal-conditioned finite-memory transfer.

All exhaustive search is confined to small toy models in tests. Production
terminal-box DP remains a relaxed-path counting method, not full SAW enumeration.
"""
from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.adaptive import adaptive_refine, verify_refinement_run
from foldbounds.isoperimetric import minimum_final_box_sum
from foldbounds.lattice import ROOT, build_frontier
from foldbounds.lookahead import verify_lookahead
from foldbounds.transfer_bound import (
    nonbacktracking_mass_upper, recent_memory_contact_spectrum,
)

STEPS = ((1, 0), (-1, 0), (0, 1), (0, -1))


def full_walks(prefix, n):
    paths = [prefix]
    for _ in range(len(prefix), n):
        next_paths = []
        for path in paths:
            x, y = path[-1]
            for dx, dy in STEPS:
                site = (x + dx, y + dy)
                if site not in path:
                    next_paths.append(path + (site,))
        paths = next_paths
    return paths


def relaxed_walks(prefix, n, memory):
    paths = [prefix]
    occupied = set(prefix)
    for _ in range(len(prefix), n):
        next_paths = []
        for path in paths:
            x, y = path[-1]
            for dx, dy in STEPS:
                site = (x + dx, y + dy)
                if site not in occupied and site not in path[-memory:]:
                    next_paths.append(path + (site,))
        paths = next_paths
    return paths


def all_contact_score(seq, path):
    return sum(
        seq[i] == seq[j] == 'H'
        and abs(path[i][0] - path[j][0])
        + abs(path[i][1] - path[j][1]) == 1
        for j in range(len(seq)) for i in range(j - 1)
    )


def prefix_cross_local(seq, path, m, memory):
    cross = local = 0
    for j in range(m, len(path)):
        for i in range(j - 1):
            if seq[i] != 'H' or seq[j] != 'H':
                continue
            if abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) != 1:
                continue
            if i < m:
                cross += 1
            elif j - i <= memory:
                local += 1
    return cross, local


@pytest.mark.parametrize('sequence', ('HHHHHHH', 'HPHHPHPH', 'PPPPPPPP', 'HHHHHHHH'))
@pytest.mark.parametrize('memory', (2, 4, 6))
@pytest.mark.parametrize('prefix', (ROOT, ROOT + ((2, 0),), ROOT + ((2, 0), (2, 1))))
@pytest.mark.parametrize('saturation', (0, 1, 2))
def test_terminal_joint_histogram_agrees_with_full_trajectory_reference(
        sequence, memory, prefix, saturation):
    counts = Counter()
    for path in relaxed_walks(prefix, len(sequence), memory):
        cross, local = prefix_cross_local(sequence, path, len(prefix), memory)
        counts[(cross, min(local, saturation), path[-1])] += 1
    actual = dict(recent_memory_contact_spectrum(
        sequence, prefix, memory, local_cap=saturation, with_endpoint=True))
    assert actual == dict(counts)
    marginal = Counter()
    for (cross, local, _), value in counts.items():
        marginal[(cross, local)] += value
    plain = dict(recent_memory_contact_spectrum(
        sequence, prefix, memory, local_cap=saturation))
    assert plain == dict(marginal)


@pytest.mark.parametrize('sequence', ('HHHHHHH', 'HPHHPHPH', 'HPPHPPHH', 'PPPPPPPP'))
@pytest.mark.parametrize('q', (Fraction(1), Fraction(3, 2), Fraction(2)))
@pytest.mark.parametrize('threshold', (-1, 0, 1, 2, 4))
@pytest.mark.parametrize('memory', (2, 4, 6))
def test_endpoint_conditioning_preserves_exhaustive_mass_enclosure(
        sequence, q, threshold, memory):
    walks = full_walks(ROOT, len(sequence))
    z = event = other = Fraction(0)
    for path in walks:
        score = all_contact_score(sequence, path)
        weight = q ** score
        z += weight
        if score >= threshold:
            event += weight
        else:
            other += weight
    old = nonbacktracking_mass_upper(sequence, ROOT, q, threshold,
                                    suffix_memory=memory)
    new = nonbacktracking_mass_upper(sequence, ROOT, q, threshold,
                                    suffix_memory=memory, terminal_box=True)
    assert z <= new.partition <= old.partition
    assert event <= new.event <= old.event
    assert other <= new.complement <= old.complement
    assert new.relaxed_walk_count == old.relaxed_walk_count


@pytest.mark.parametrize('sequence', ('HHHHHHHH', 'HHPHHHPH', 'HPHHPHHP'))
@pytest.mark.parametrize('memory', (2, 4))
def test_integrated_interval_nesting_and_adaptive_replay(sequence, memory):
    cert = build_frontier(sequence, Fraction(3, 2), 2, 2)
    settings = dict(lookahead=0, use_matching=True, use_transfer=True,
                    use_monotone=True, transfer_memory=memory)
    old = verify_lookahead(cert, terminal_box=False, **settings)
    new = verify_lookahead(cert, terminal_box=True, **settings)
    assert new.partition_lower == old.partition_lower
    assert new.partition_upper <= old.partition_upper
    assert new.event_upper <= old.event_upper
    assert new.complement_upper <= old.complement_upper
    assert new.probability_lower >= old.probability_lower
    assert new.probability_upper <= old.probability_upper
    run = adaptive_refine(cert, steps=2, policy='mass', terminal_box=True,
                          use_matching=True, use_transfer=True,
                          transfer_memory=memory)
    assert run.terminal_box
    assert verify_refinement_run(run) == run.bounds
    with pytest.raises(ValueError):
        verify_refinement_run(replace(run, terminal_box=False))


def test_terminal_isoperimetric_bound_contains_all_complete_small_walk_scores():
    for n in range(2, 9):
        for walk in full_walks(ROOT, n):
            # Count all nonconsecutive lattice contacts, not just H-H.
            contacts = sum(
                abs(walk[i][0] - walk[j][0]) +
                abs(walk[i][1] - walk[j][1]) == 1
                for j in range(n) for i in range(j-1))
            x0, x1 = min(p[0] for p in ROOT), max(p[0] for p in ROOT)
            y0, y1 = min(p[1] for p in ROOT), max(p[1] for p in ROOT)
            x, y = walk[-1]
            width = max(x1, x) - min(x0, x) + 1
            height = max(y1, y) - min(y0, y) + 1
            cap = max(0, n + 1 - minimum_final_box_sum(n, width, height))
            assert contacts <= cap


@pytest.mark.parametrize('bad', (0, 'yes', None, 3.5))
def test_invalid_terminal_geometry_controls_fail_closed(bad):
    with pytest.raises(ValueError):
        nonbacktracking_mass_upper('HHHHH', ROOT, Fraction(2), 1,
                                   terminal_box=bad)
    with pytest.raises(ValueError):
        recent_memory_contact_spectrum('HHHHH', ROOT, 4,
                                       with_endpoint=bad)
    with pytest.raises(ValueError):
        verify_lookahead(build_frontier('HHHHH', Fraction(2), 1, 2),
                         use_transfer=True, terminal_box=bad)


def test_terminal_geometry_requires_enabled_transfer():
    with pytest.raises(ValueError):
        verify_lookahead(build_frontier('HHHHH', Fraction(2), 1, 2),
                         terminal_box=True)
