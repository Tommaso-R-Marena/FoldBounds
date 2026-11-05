"""Independent oracle tests for partially directed disjoint SAW spectra."""
from collections import Counter
from fractions import Fraction
from itertools import product
from random import Random

import pytest

from foldbounds.lattice import build_frontier
from foldbounds.partial_directed import (
    partially_directed_spectrum, is_partially_directed_path,
)


def enumerate_all_walks(n):
    walks = [((0, 0), (1, 0))]
    for _ in range(2, n):
        next_walks = []
        for walk in walks:
            x, y = walk[-1]
            for dest in ((x+1, y), (x-1, y), (x, y+1), (x, y-1)):
                if dest not in walk:
                    next_walks.append(walk + (dest,))
        walks = next_walks
    return walks


def independent_contact_score(seq, walk):
    return sum(
        seq[i] == seq[j] == "H"
        and abs(walk[i][0] - walk[j][0]) + abs(walk[i][1] - walk[j][1]) == 1
        for i in range(len(seq)) for j in range(i+2, len(seq))
    )


def independently_classify(walk):
    """Reference implementation using literal direction words and segments.

    The first move along +x is already present in the root. No DP geometry,
    no relative-coordinate contact scoring, and no production classification.
    """
    directions = []
    for prev, current in zip(walk, walk[1:]):
        d = (current[0]-prev[0], current[1]-prev[1])
        directions.append({(1, 0): "E", (-1, 0): "W",
                           (0, 1): "N", (0, -1): "S"}[d])
    words = ''.join(directions)

    def no_reverse_within_layers(word, forward, side1, side2):
        if any(letter not in (forward, side1, side2) for letter in word):
            return False
        for layer in word.split(forward):
            if side1 in layer and side2 in layer:
                return False
        return True

    e = no_reverse_within_layers(words[1:], 'E', 'N', 'S')
    n = 'W' in words and no_reverse_within_layers(words, 'N', 'E', 'W')
    s = 'W' in words and no_reverse_within_layers(words, 'S', 'E', 'W')
    assert int(e) + int(n) + int(s) <= 1
    return bool(e or n or s)


def test_exhaustive_short_all_sequences_and_arbitrary_prefixes():
    for n in range(2, 9):
        walks = enumerate_all_walks(n)
        prefixes = {p[:depth] for p in walks
                    for depth in (2, min(4,n), min(5,n), n)}
        for walk in walks:
            assert is_partially_directed_path(walk) == independently_classify(walk)
        sequences = ('H'*n, 'P'*n, ('HP'*n)[:n], ('PH'*n)[:n])
        for seq in sequences:
            for prefix in prefixes:
                got = Counter(dict(partially_directed_spectrum(seq, prefix)))
                expect = Counter(independent_contact_score(seq, walk)
                                 for walk in walks
                                 if walk[:len(prefix)] == prefix
                                 and independently_classify(walk))
                assert got == expect, (seq, prefix, got, expect)


def test_terminal_prefix_scores_and_nonmatching_yield_empty_spectra():
    n = 8
    walks = enumerate_all_walks(n)
    seq = 'HHHPPHHH'
    for walk in walks:
        hist = dict(partially_directed_spectrum(seq, walk))
        assert sum(hist.values()) == int(independently_classify(walk))
        if independently_classify(walk):
            assert hist == {independent_contact_score(seq, walk):1}


def test_disjoint_directions_grow_exponentially_without_enumeration():
    for n in (30, 40, 60):
        hist = dict(partially_directed_spectrum('H'*n))
        assert len(hist) > n // 5
        assert sum(hist.values()) >= 2**(n-2)
        assert min(hist) == 0
        assert max(hist) > 0


def test_adversarial_prefix_inputs_rejected():
    seq='HHHHH'
    invalid = [((0, 0), (0, 1)),
               ((0, 0), (1, 0), (0, 0)),
               ((0, 0), (1, 0), (2, 0), (4, 0)),
               ((0, 0), (1, 0), (2., 0.)),
               ((0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0))]
    for prefix in invalid:
        with pytest.raises((ValueError, TypeError)):
            partially_directed_spectrum(seq, prefix)


def test_random_exact_mass_lower_and_root_concordance():
    rng=Random(38842)
    for _ in range(60):
        n=rng.randint(3,9)
        seq=''.join(rng.choice('HP') for _ in range(n))
        q=rng.choice((Fraction(1),Fraction(3,2),Fraction(2)))
        exact=sum((q**independent_contact_score(seq,w) for w in enumerate_all_walks(n)), Fraction(0))
        hist=dict(partially_directed_spectrum(seq))
        known=sum((count*q**score for score,count in hist.items()),Fraction(0))
        assert Fraction(0)<known<=exact


def test_east_family_matches_independent_exact_generating_function():
    """A(z)=(1+z)/(1-2z-z²) gives 1,3,7,17,... exactly."""
    from foldbounds.partial_directed import _rooted_axis_spectrum
    expected=[1,3]
    for _ in range(2,43):
        expected.append(2*expected[-1]+expected[-2])
    root=((0,0),(1,0))
    for r in range(43):
        seq='P'*(r+2)
        got=sum(_rooted_axis_spectrum(seq,root,'E').values())
        assert got==expected[r]


def test_analytic_rooted_partition_lower_bound_scales_and_is_sound():
    from foldbounds.partial_directed import (
        east_directed_path_count, universal_hp_partition_lower,
    )
    assert [east_directed_path_count(r) for r in range(7)] == [1,3,7,17,41,99,239]
    for n in range(2,10):
        seq='H' * n
        for q in (Fraction(1),Fraction(3,2),Fraction(2)):
            full=sum((q**independent_contact_score(seq,p) for p in enumerate_all_walks(n)),Fraction(0))
            assert universal_hp_partition_lower(seq,q) <= full
    assert east_directed_path_count(1000) > 2**1000
    for bad in (-1,True,2.5):
        with pytest.raises(ValueError):
            east_directed_path_count(bad)
    with pytest.raises(ValueError):
        universal_hp_partition_lower('HHH', Fraction(1,2))
