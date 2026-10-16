"""Independent finite-memory relaxed-walk and soundness oracles.

The reference enumerator below keeps whole trajectories and scores every
pair from scratch. The production DP only retains K sites and aggregated
(integer) multiplicities. These tests do not replace a Lean refinement proof.
"""
from collections import Counter
from fractions import Fraction

import pytest

from foldbounds.lattice import ROOT
from foldbounds.transfer_bound import (
    nonbacktracking_mass_upper,
    recent_memory_contact_spectrum,
    recent_memory_cross_spectrum,
)

DIRECTIONS = ((1, 0), (0, 1), (-1, 0), (0, -1))


def relaxed_paths(prefix, n, memory):
    """Reference enumeration independent of transfer-state aggregation."""
    candidates = [prefix]
    occupied = set(prefix)
    for _ in range(len(prefix), n):
        following = []
        for path in candidates:
            x, y = path[-1]
            for dx, dy in DIRECTIONS:
                new = (x + dx, y + dy)
                if new in occupied or new in path[-memory:]:
                    continue
                following.append(path + (new,))
        candidates = following
    return candidates


def true_paths(prefix, n):
    candidates = [prefix]
    for _ in range(len(prefix), n):
        following = []
        for path in candidates:
            x, y = path[-1]
            for dx, dy in DIRECTIONS:
                new = (x + dx, y + dy)
                if new not in path:
                    following.append(path + (new,))
        candidates = following
    return candidates


def independent_scores(sequence, path, m, memory):
    """Calculate all eligible contacts from the complete coordinate word."""
    cross = local = total = 0
    for j in range(len(path)):
        if sequence[j] != 'H':
            continue
        for i in range(j - 1):
            if sequence[i] != 'H':
                continue
            if abs(path[j][0] - path[i][0]) + abs(path[j][1] - path[i][1]) != 1:
                continue
            total += 1
            if i < m <= j:
                cross += 1
            elif m <= i and j - i <= memory:
                local += 1
    return cross, local, total


@pytest.mark.parametrize("sequence", (
    "HHHHHHHH", "HPHHPHHP", "HPPHHHHP", "PPPPPPPP"
))
@pytest.mark.parametrize("memory", (2, 3, 4, 6))
@pytest.mark.parametrize("prefix", (
    ROOT,
    ROOT + ((2, 0),),
    ROOT + ((2, 0), (2, 1)),
))
def test_joint_spectrum_matches_full_trajectory_enumerator(sequence, memory, prefix):
    counts = Counter()
    for full in relaxed_paths(prefix, len(sequence), memory):
        cross, local, _ = independent_scores(sequence, full, len(prefix), memory)
        counts[(cross, local)] += 1
    joint = dict(recent_memory_contact_spectrum(sequence, prefix, memory))
    assert joint == dict(counts)
    marginal = Counter()
    for (cross, _), count in counts.items():
        marginal[cross] += count
    assert dict(recent_memory_cross_spectrum(sequence, prefix, memory)) == dict(marginal)


@pytest.mark.parametrize("sequence", ("HHHHHHHH", "HPHHPHHP", "HPPHHHHP"))
@pytest.mark.parametrize("memory", (2, 3, 4, 6))
@pytest.mark.parametrize("threshold", (0, 1, 2, 3))
def test_real_conformations_always_enclosed(sequence, memory, threshold):
    q = Fraction(3, 2)
    prefix = ROOT
    z = event = complement = Fraction(0)
    for full in true_paths(prefix, len(sequence)):
        _, _, score = independent_scores(sequence, full, len(prefix), memory)
        weight = q ** score
        z += weight
        if score >= threshold:
            event += weight
        else:
            complement += weight
    bounds = nonbacktracking_mass_upper(sequence, prefix, q, threshold,
                                       suffix_memory=memory)
    assert z <= bounds.partition
    assert event <= bounds.event
    assert complement <= bounds.complement


def test_recent_contacts_strictly_prune_complement_without_reweighting():
    seq = 'H' * 8
    q = Fraction(2)
    memory = 4
    threshold = 1
    bounds = nonbacktracking_mass_upper(seq, ROOT, q, threshold,
                                       suffix_memory=memory)
    # The previous classifier used prefix/suffix contacts alone, ignoring
    # genuine local suffix/suffix contacts when classifying complement mass.
    # Compare equal upper weights under the same fixed memory.
    from foldbounds.transfer_bound import future_future_cap
    from foldbounds.isoperimetric import prefix_contact_cap
    ff = future_future_cap(seq, len(ROOT))
    perimeter = prefix_contact_cap(ROOT, len(seq))
    previous_complement = Fraction(0)
    for (cross, _local), multiplicity in recent_memory_contact_spectrum(
            seq, ROOT, memory):
        cap = min(ff, max(0, perimeter - cross))
        if cross < threshold:
            previous_complement += multiplicity * q ** (cross + cap)
    assert bounds.complement < previous_complement
    assert bounds.event >= 0
    assert bounds.partition >= 0


@pytest.mark.parametrize("memory", (2, 3, 4, 6))
def test_contact_score_floor_never_exceeds_true_score(memory):
    for sequence in ("HHHHHHHH", "HPHHPHHP"):
        for path in true_paths(ROOT, len(sequence)):
            cross, local, score = independent_scores(
                sequence, path, len(ROOT), memory)
            assert cross + local <= score
