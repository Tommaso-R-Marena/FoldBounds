from fractions import Fraction

import pytest

from foldbounds.lattice import build_frontier
from foldbounds.lean_export_v2 import certificate_to_lean_text_v2, raw_lean_text_v2


@pytest.mark.parametrize("seq,depth,leaf,memory", [
    ("HHHH", 2, "M", 0), ("HPHPH", 3, "M", 2), ("HHPPHH", 4, "M", 4),
    ("HHPPHH", 4, "T", 0), ("HHPPHH", 3, "B", 0)])
def test_v2_text_shape(seq, depth, leaf, memory):
    cert = build_frontier(seq, Fraction(3, 2), 1, depth)
    text = certificate_to_lean_text_v2(cert, leaf=leaf, memory=memory)
    lines = text.splitlines()
    assert lines[:4] == ["FOLDBOUNDS-FRONTIER 2", f"sequence {seq}", "q 3/2", "threshold 1"]
    assert lines[-1] == "end"
    assert len(lines) == 5 + len(cert.frontier)
    kinds = {"E", "D", f"M{memory}"} if leaf == "M" else {"E", "D", leaf}
    assert all(line.split(" ")[0] in kinds for line in lines[4:-1])


def test_v2_enum_threshold():
    cert = build_frontier("HHHHHH", Fraction(2), 2, 3)
    text = certificate_to_lean_text_v2(cert, leaf="M", memory=2, enum_remaining=10)
    body = text.splitlines()[4:-1]
    assert body and all(line.split(" ")[0] in {"X", "E", "D"} for line in body)


@pytest.mark.parametrize("leaf,memory", [("Q", 0), ("M", -1)])
def test_v2_rejects_bad_options(leaf, memory):
    cert = build_frontier("HHHH", Fraction(2), 1, 2)
    with pytest.raises(ValueError):
        certificate_to_lean_text_v2(cert, leaf=leaf, memory=memory)


def test_v2_witness_must_extend_node():
    cert = build_frontier("HHHHH", Fraction(2), 1, 3)
    node = next(n for n in cert.frontier if n.kind == "bound")
    bogus = ((0, 0), (1, 0), (1, -1), (2, -1), (3, -1))
    if bogus[:len(node.path)] == node.path:
        bogus = ((0, 0), (1, 0), (1, 1), (2, 1), (3, 1))
    with pytest.raises(ValueError):
        certificate_to_lean_text_v2(cert, leaf="M", memory=1, witnesses={node.path: [bogus]})


def test_v2_raw_text():
    assert raw_lean_text_v2("HH", Fraction(2), 0, ["X ."]) == (
        "FOLDBOUNDS-FRONTIER 2\nsequence HH\nq 2\nthreshold 0\nX .\nend\n")
