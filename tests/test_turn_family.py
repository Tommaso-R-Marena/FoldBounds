"""Independent complete-walk oracle for restricted two-turn witness families."""
from collections import Counter
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.adaptive import adaptive_refine, verify_refinement_run
from foldbounds.lattice import ROOT, build_frontier
from foldbounds.lookahead import verify_lookahead
from foldbounds.turn_family import turn_family_mass

MOVES = ((1, 0), (0, 1), (-1, 0), (0, -1))


def all_self_avoiding_extensions(prefix, n):
    paths = [prefix]
    for _ in range(len(prefix), n):
        next_paths = []
        for path in paths:
            x, y = path[-1]
            for dx, dy in MOVES:
                nxt = (x + dx, y + dy)
                if nxt not in path:
                    next_paths.append(path + (nxt,))
        paths = next_paths
    return paths


def independent_contacts(seq, path):
    return sum(
        seq[i] == 'H' and seq[j] == 'H'
        and abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) == 1
        for i in range(len(path))
        for j in range(i + 2, len(path))
    )


def independent_turns(path, m):
    suffix = [
        (path[i][0] - path[i - 1][0], path[i][1] - path[i - 1][1])
        for i in range(m, len(path))
    ]
    return sum(a != b for a, b in zip(suffix, suffix[1:]))


@pytest.mark.parametrize("sequence", (
    "HHHHHHHH", "HPHHPHHP", "HPPHHHHP", "PPPPPPPP",
))
@pytest.mark.parametrize("prefix", (ROOT, ROOT + ((2, 0), (2, 1))))
@pytest.mark.parametrize("max_turns", (0, 1, 2))
@pytest.mark.parametrize("threshold", (0, 1, 2))
def test_low_turn_mass_equals_independent_complete_walk_oracle(
        sequence, prefix, max_turns, threshold):
    q = Fraction(3, 2)
    event = other = Fraction(0)
    count = 0
    for full in all_self_avoiding_extensions(prefix, len(sequence)):
        if independent_turns(full, len(prefix)) > max_turns:
            continue
        count += 1
        weight = q ** independent_contacts(sequence, full)
        if independent_contacts(sequence, full) >= threshold:
            event += weight
        else:
            other += weight
    result = turn_family_mass(
        sequence, prefix, q, threshold, max_turns=max_turns)
    assert result.event == event
    assert result.complement == other
    assert result.distinct_paths == count


def test_mixed_sequence_event_now_has_strict_positive_constructive_mass():
    sequence = "HPHHPHHPHHHH"
    turn = turn_family_mass(sequence, ROOT, Fraction(3, 2), 2)
    assert turn.event == Fraction(459, 8)
    cert = build_frontier(sequence, Fraction(3, 2), 2, 2)
    old = verify_lookahead(
        cert, lookahead=0, use_matching=True, use_transfer=True,
        use_monotone=True, transfer_memory=4)
    new = verify_lookahead(
        cert, lookahead=0, use_matching=True, use_transfer=True,
        use_monotone=True, use_turn_family=True, transfer_memory=4)
    assert old.event_lower == 0
    assert new.event_lower >= turn.event > 0
    assert new.probability_lower > old.probability_lower
    assert new.partition_upper == old.partition_upper
    assert new.event_upper == old.event_upper
    assert new.complement_upper == old.complement_upper
    assert new.probability_upper <= old.probability_upper


@pytest.mark.parametrize("sequence", ("H"*8, "HPHHPHHP", "H"*12))
def test_integrated_turn_family_never_loosens_mass_enclosures(sequence):
    cert = build_frontier(sequence, Fraction(2), 1, 2)
    old = verify_lookahead(cert, lookahead=0, use_matching=True,
                           use_transfer=True, use_monotone=True)
    new = verify_lookahead(cert, lookahead=0, use_matching=True,
                           use_transfer=True, use_monotone=True,
                           use_turn_family=True)
    assert new.event_lower >= old.event_lower
    assert new.complement_lower >= old.complement_lower
    assert new.partition_lower >= old.partition_lower
    assert new.partition_upper == old.partition_upper
    assert new.probability_lower >= old.probability_lower
    assert new.probability_upper <= old.probability_upper


def test_adaptive_turn_option_is_replayed_and_tampering_rejected():
    cert = build_frontier("HPHHPHHPHHHH", Fraction(3, 2), 2, 2)
    run = adaptive_refine(cert, steps=2, use_matching=True,
                          use_transfer=True, use_monotone=True,
                          use_turn_family=True, policy="mass")
    assert run.use_turn_family
    assert run.bounds == verify_refinement_run(run)
    with pytest.raises(ValueError):
        verify_refinement_run(replace(run, use_turn_family=False))


@pytest.mark.parametrize("bad", (-1, 3, True, False, 2.5, "2", None))
def test_invalid_turn_budget_rejected(bad):
    with pytest.raises(ValueError):
        turn_family_mass("HHHH", ROOT, Fraction(2), 1, max_turns=bad)
