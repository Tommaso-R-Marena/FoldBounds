"""Export Python frontier certificates to the Lean frontier-checker text format.

The exporter is UNTRUSTED: the Lean checker (``FoldBounds.FrontierCheck``)
re-validates every node.  Only the model parameters (sequence, q, threshold)
must be checked by a human against the echo printed by the Lean runner.

Lean direction encoding: 0 = +x, 1 = -x, 2 = +y, 3 = -y.  (Python's
``DIRECTIONS`` tuple uses a different order; paths are converted step by step
from coordinates, so the order never matters.)
"""
from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

from .lattice import ROOT, Certificate, Path

_STEP_TO_LEAN = {(1, 0): "0", (-1, 0): "1", (0, 1): "2", (0, -1): "3"}
_LEAN_TO_STEP = {v: k for k, v in _STEP_TO_LEAN.items()}


def path_to_word(path: Path) -> str:
    """Direction word (after the fixed first bond) of an anchored path."""
    if tuple(path[:2]) != ROOT:
        raise ValueError("path must start with the canonical root (0,0),(1,0)")
    out = []
    for (x0, y0), (x1, y1) in zip(path[1:], path[2:]):
        step = (x1 - x0, y1 - y0)
        if step not in _STEP_TO_LEAN:
            raise ValueError(f"non-unit lattice step {step}")
        out.append(_STEP_TO_LEAN[step])
    return "".join(out) or "."


def word_to_path(word: str) -> Path:
    """Inverse of :func:`path_to_word` (no self-avoidance check)."""
    path = list(ROOT)
    if word != ".":
        for ch in word:
            dx, dy = _LEAN_TO_STEP[ch]
            x, y = path[-1]
            path.append((x + dx, y + dy))
    return tuple(path)


def _fraction(q: Fraction) -> str:
    return f"{q.numerator}" if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def certificate_to_lean_text(cert: Certificate, *, transfer: bool = False,
                             witnesses: dict[Path, Sequence[Path]] | None = None) -> str:
    """Serialise a Python frontier certificate.

    ``transfer`` marks every bounded node for the verified transfer bound.
    ``witnesses`` optionally maps a bounded node path to complete paths
    (proposed by any untrusted search) extending it.
    """
    witnesses = witnesses or {}
    lines = ["FOLDBOUNDS-FRONTIER 1", f"sequence {cert.sequence}",
             f"q {_fraction(cert.q)}", f"threshold {cert.event_at_least}"]
    for node in cert.frontier:
        word = path_to_word(node.path)
        if node.kind == "exact":
            lines.append(f"E {word}")
        elif node.kind == "dead":
            lines.append(f"D {word}")
        elif node.kind == "bound":
            kind = "T" if transfer else "B"
            prefix = word if word != "." else ""
            sufs = []
            for full in witnesses.get(node.path, ()):
                fw = path_to_word(full)
                fw = fw if fw != "." else ""
                if not fw.startswith(prefix):
                    raise ValueError("witness does not extend its node")
                sufs.append(fw[len(prefix):] or ".")
            lines.append(" ".join([kind, word, *sufs]))
        else:
            raise ValueError(f"unknown node kind {node.kind!r}")
    lines.append("end")
    return "\n".join(lines) + "\n"


def raw_lean_text(sequence: str, q: Fraction, threshold: int,
                  nodes: Iterable[str]) -> str:
    """Build a text certificate from raw node lines (for adversarial tests)."""
    return "\n".join(["FOLDBOUNDS-FRONTIER 1", f"sequence {sequence}",
                      f"q {_fraction(q)}", f"threshold {threshold}", *nodes, "end"]) + "\n"
