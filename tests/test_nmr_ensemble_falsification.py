"""Real NMR multi-conformer falsification: exact coordinate-model contact drift."""
from __future__ import annotations

from math import isqrt
from pathlib import Path
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data"))
from experimental_ca_packed import STRUCTURES, decode
from experimental_ca_1LCD_ensembles import MODEL_PACKED, SOURCE_SHA
from foldbounds.sparse_distance_certificate import make_certificate

def squared(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b))

@pytest.mark.parametrize("model,expected_changes,expected_losses", [
    (2, 12, 9),
    (3, 15, 11),
])
def test_exact_nmr_contact_change_and_anchor_assumption_failure(
        model, expected_changes, expected_losses):
    base = decode(STRUCTURES["1LCD"]["packed"])
    other = decode(MODEL_PACKED[model])
    assert len(base) == len(other) == 51
    assert SOURCE_SHA == STRUCTURES["1LCD"]["sha"]

    changed = lost = 0
    for i in range(51):
        for j in range(i + 4, 51):
            was = squared(base[i], base[j]) <= 8000 ** 2
            now = squared(other[i], other[j]) <= 8000 ** 2
            changed += int(was != now)
            lost += int(was and not now)
    assert (changed, lost) == (expected_changes, expected_losses)

    selected = make_certificate(base).anchor_indices
    # Certify without floats that at least one mutual anchor distance
    # changes by >1 Angstrom between these deposited NMR conformers.
    def differs_by_more_than_one_angstrom(s, t):
        # |sqrt(s)-sqrt(t)| > 1000 if
        # |s-t| > 1000*(sqrt(s)+sqrt(t)).
        # isqrt(s)+isqrt(t)+2 strictly bounds the denominator.
        return abs(s-t) > 1000 * (isqrt(s) + isqrt(t) + 2)

    assert any(
        differs_by_more_than_one_angstrom(
            squared(base[i],base[j]), squared(other[i],other[j]))
        for x, i in enumerate(selected) for j in selected[x + 1:]
    )

def test_nmr_dataset_has_no_independent_holdout_claim():
    assert set(MODEL_PACKED) == {2, 3}
    assert all(len(decode(v)) == 51 for v in MODEL_PACKED.values())
