from fractions import Fraction

import pytest

from foldbounds.lattice import build_frontier
from foldbounds.lean_export import (certificate_to_lean_text, path_to_word, raw_lean_text,
                                    word_to_path)


@pytest.mark.parametrize("word", [".", "0", "02", "0213", "2210", "3300"])
def test_word_path_roundtrip(word):
    assert path_to_word(word_to_path(word)) == word


def test_rejects_wrong_root():
    with pytest.raises(ValueError):
        path_to_word(((0, 0), (0, 1)))


def test_rejects_non_unit_step():
    with pytest.raises(ValueError):
        path_to_word(((0, 0), (1, 0), (3, 0)))


@pytest.mark.parametrize("seq,depth", [("HHHH", 2), ("HPHPH", 3), ("HHPPHH", 4)])
def test_certificate_text_shape(seq, depth):
    cert = build_frontier(seq, Fraction(3, 2), 1, depth)
    text = certificate_to_lean_text(cert, transfer=True)
    lines = text.splitlines()
    assert lines[:4] == ["FOLDBOUNDS-FRONTIER 1", f"sequence {seq}", "q 3/2", "threshold 1"]
    assert lines[-1] == "end"
    assert len(lines) == 5 + len(cert.frontier)
    assert all(line.split(" ")[0] in {"E", "T", "D"} for line in lines[4:-1])


def test_raw_text():
    assert raw_lean_text("HH", Fraction(2), 0, ["E ."]) == (
        "FOLDBOUNDS-FRONTIER 1\nsequence HH\nq 2\nthreshold 0\nE .\nend\n")
