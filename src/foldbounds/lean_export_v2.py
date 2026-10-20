"""Export frontier certificates to the version-2 Lean text format.

UNTRUSTED producer.  The Lean checker ``FoldBounds.FrontierNext`` (runner
``foldcheck2``) re-validates every node; soundness of accepted intervals is
the Lean theorem ``FoldBounds.FrontierNext.checkFlat2_sound``.  Only the
echoed model parameters (sequence, q, threshold) must be compared by a human
with the intended ones.

Leaf kinds (version 2):

* ``E w``           exact complete conformation
* ``D w``           dead prefix
* ``X w``           exhaustive local evaluation inside the verifier
* ``B w [wit..]``   legacy basic bound          (version-1 semantics)
* ``T w [wit..]``   legacy transfer bound       (version-1 semantics)
* ``M<k> w [wit..]`` memory-k transfer + geometric caps + PD families

This module never modifies ``lean_export.py`` (version 1 stays available).
"""
from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

from .lattice import Certificate, Path
from .lean_export import path_to_word, _fraction

HEADER_V2 = "FOLDBOUNDS-FRONTIER 2"


def _suffixes(prefix_word: str, full_paths: Iterable[Path]) -> list[str]:
    prefix = prefix_word if prefix_word != "." else ""
    out = []
    for full in full_paths:
        fw = path_to_word(full)
        fw = fw if fw != "." else ""
        if not fw.startswith(prefix):
            raise ValueError("witness does not extend its node")
        out.append(fw[len(prefix):] or ".")
    return out


def leaf_line(kind: str, word: str, suffixes: Sequence[str] = ()) -> str:
    return " ".join([kind, word, *suffixes])


def certificate_to_lean_text_v2(cert: Certificate, *, leaf: str = "M", memory: int = 2,
                                witnesses: dict[Path, Sequence[Path]] | None = None,
                                enum_remaining: int = -1) -> str:
    """Serialise a Python frontier certificate in the version-2 format.

    ``leaf`` selects the kind used for bounded nodes: ``"B"``, ``"T"`` or ``"M"``
    (memory ``memory``).  Bounded nodes with at most ``enum_remaining`` free
    directions left are emitted as ``X`` (exhaustive local evaluation).
    """
    if leaf not in ("B", "T", "M"):
        raise ValueError("leaf must be B, T or M")
    if type(memory) is not int or memory < 0:
        raise ValueError("memory must be a natural number")
    witnesses = witnesses or {}
    n = len(cert.sequence)
    lines = [HEADER_V2, f"sequence {cert.sequence}", f"q {_fraction(cert.q)}",
             f"threshold {cert.event_at_least}"]
    for node in cert.frontier:
        word = path_to_word(node.path)
        if node.kind == "exact":
            lines.append(f"E {word}")
        elif node.kind == "dead":
            lines.append(f"D {word}")
        elif node.kind == "bound":
            remaining = n - len(node.path)
            if remaining <= enum_remaining:
                lines.append(f"X {word}")
                continue
            kind = f"M{memory}" if leaf == "M" else leaf
            lines.append(leaf_line(kind, word, _suffixes(word, witnesses.get(node.path, ()))))
        else:
            raise ValueError(f"unknown node kind {node.kind!r}")
    lines.append("end")
    return "\n".join(lines) + "\n"


def raw_lean_text_v2(sequence: str, q: Fraction, threshold: int, nodes: Iterable[str]) -> str:
    """Build a version-2 text certificate from raw node lines (adversarial tests)."""
    return "\n".join([HEADER_V2, f"sequence {sequence}", f"q {_fraction(q)}",
                      f"threshold {threshold}", *nodes, "end"]) + "\n"
