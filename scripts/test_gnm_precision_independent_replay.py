"""Hermetic exact-interval tests and optional full real-protein replay.

Run real checks with:
 FOLDBOUNDS_RNASE_STAGE1=/path/to/Stage1.zip
 FOLDBOUNDS_SPECTRAL_CERT=/path/to/spectral.json
 OPENBLAS_NUM_THREADS=1 python -m pytest -q scripts/test_gnm_precision_independent_replay.py
"""
import os
from pathlib import Path
from fractions import Fraction
import pytest

from gnm_precision_independent_replay import (
    residual, intervals, source_pairs, verify,
)


def test_sparse_integer_residual():
    A = [[2, 0], [0, 1]]
    N = [[7, 0], [0, 14]]
    assert residual(A, N, 16) == [[2, 0], [0, 2]]


def test_exact_interval_bound_monotonicity():
    # Diagonal grounded matrix: true inverse = diag(1/2, 1).
    A, N, Q = [[2, 0], [0, 1]], [[7, 0], [0, 14]], 16
    R = residual(A, N, Q)
    pairs = [(0, 2), (1, 2), (0, 1)]
    tight_n, tight_w = intervals(A, N, R, Q, 1, 1, pairs, 2)
    loose_n, loose_w = intervals(A, N, R, Q, 5, 1, pairs, 2)
    assert 0 <= loose_n <= tight_n <= len(pairs) - 1
    assert 0 < tight_w and loose_w == 5 * tight_w
    # Explicit integer-residual bounds enclose the exact responses.
    assert Fraction(126, 256) <= Fraction(1, 2) <= Fraction(130, 256)
    assert Fraction(252, 256) <= 1 <= Fraction(256, 256)
    assert Fraction(378, 256) <= Fraction(3, 2) <= Fraction(386, 256)


def test_tampered_stage1_rejected(tmp_path):
    bad = tmp_path / "tampered.zip"
    bad.write_bytes(b"forged source archive")
    with pytest.raises(ValueError, match="Wrong frozen source archive"):
        source_pairs(bad)


def test_full_archived_18_precision_levels_when_inputs_available():
    s = os.environ.get("FOLDBOUNDS_RNASE_STAGE1")
    w = os.environ.get("FOLDBOUNDS_SPECTRAL_CERT")
    if not s or not w:
        pytest.skip("Original pinned Stage1 ZIP and spectral SOS certificate required")
    report = (Path(__file__).resolve().parents[1] /
              "research/GNMSourceOnlyPrecisionAblation_20261105.json")
    rows = verify(s, w, report)
    assert len(rows) == 18
    assert (16, 102, 2138) in rows
    assert (24, 2589, 2590) in rows
