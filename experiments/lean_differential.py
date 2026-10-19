#!/usr/bin/env python3
"""Differential and adversarial tests: Python reference vs. the Lean checker.

Runs the *Lean* frontier checker (``scripts/foldcheck.lean``) on certificates
exported from Python and compares against

1. an independent brute-force enumerator written here (no shared code with
   ``foldbounds.lattice``), and
2. the Python reference verifier ``foldbounds.lattice.verify_certificate``.

It also checks that the Lean decoder/scorer enumerates exactly the same
conformations with the same contact counts, and that mutated (malformed or
dishonest) certificates are rejected.

THIS IS TESTING, NOT PROOF.  Agreement on these cases does not establish a
cross-language refinement theorem.

Usage (from the Lake project root that builds ``FoldBounds``):

    python3 FoldBounds-main/experiments/lean_differential.py \
        --lake-root . --runner FoldBounds-main/scripts/foldcheck.lean
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from fractions import Fraction
from itertools import product
from pathlib import Path as FsPath

HERE = FsPath(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

from foldbounds.lattice import (Certificate, Node, build_frontier,  # noqa: E402
                                exact_spectrum, verify_certificate, CertificateError)
from foldbounds.lean_export import (certificate_to_lean_text, path_to_word,  # noqa: E402
                                    raw_lean_text)
from foldbounds.witnesses import propose_dfs_witnesses  # noqa: E402

STEPS = {"0": (1, 0), "1": (-1, 0), "2": (0, 1), "3": (0, -1)}


# ---------------------------------------------------------------- oracle ---
def brute_force(sequence: str):
    """Independent enumerator: all direction words, naive all-pairs contacts."""
    n = len(sequence)
    out = []
    for word in product("0123", repeat=n - 2):
        path = [(0, 0), (1, 0)]
        ok = True
        for ch in word:
            dx, dy = STEPS[ch]
            x, y = path[-1]
            nxt = (x + dx, y + dy)
            if nxt in path:
                ok = False
                break
            path.append(nxt)
        if not ok:
            continue
        k = sum(1 for i in range(n) for j in range(i + 2, n)
                if sequence[i] == "H" and sequence[j] == "H"
                and abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) == 1)
        out.append((tuple(path), k))
    return out


def exact_values(sequence: str, q: Fraction, thr: int):
    confs = brute_force(sequence)
    za = sum((q ** k for _, k in confs if k >= thr), Fraction(0))
    zn = sum((q ** k for _, k in confs if k < thr), Fraction(0))
    return za + zn, za, zn, za / (za + zn)


# ------------------------------------------------------------ lean runner ---
class Lean:
    def __init__(self, root: str, runner: str):
        self.root, self.runner = root, runner

    def _run(self, args):
        res = subprocess.run(["lake", "env", "lean", "--run", self.runner, *args],
                             cwd=self.root, capture_output=True, text=True, check=True)
        return res.stdout

    def check_texts(self, texts):
        with tempfile.TemporaryDirectory() as d:
            files = []
            for i, t in enumerate(texts):
                f = os.path.join(d, f"c{i}.txt")
                FsPath(f).write_text(t)
                files.append(f)
            t0 = time.time()
            out = self._run(["check", *files])
            wall = time.time() - t0
        rows = [json.loads(line) for line in out.splitlines() if line.strip()]
        assert len(rows) == len(texts)
        return rows, wall

    def enumerate(self, sequence):
        return json.loads(self._run(["enumerate", sequence]))


def frac(s: str) -> Fraction:
    a, b = s.split("/")
    return Fraction(int(a), int(b))


# ------------------------------------------------------------------ tests ---
SEQUENCES = ["HH", "HP", "HHH", "HHHH", "HPPH", "PPPPP", "HHHHH", "HPHPH", "HHPHH",
             "HPHPHH", "HHHHHH", "PHHPHP", "HPHPHPH", "HHPPHHP", "HHHHHHH", "HPPHHPPH",
             "HHHHHHHH", "PPPPPPPP"]
QS = [Fraction(1), Fraction(3, 2), Fraction(2)]


def test_semantics(lean: Lean, report):
    mismatches = 0
    for seq in SEQUENCES:
        lean_rows = lean.enumerate(seq)
        lean_set = {(tuple(tuple(s) for s in r["sites"]), r["contacts"]) for r in lean_rows}
        bf = set(brute_force(seq))
        py_spec = exact_spectrum(seq)
        lean_spec = {}
        for _, k in lean_set:
            lean_spec[k] = lean_spec.get(k, 0) + 1
        ok = lean_set == bf and dict(sorted(lean_spec.items())) == py_spec
        mismatches += not ok
        report.append(f"semantics {seq}: {len(lean_set)} conformations, spectrum "
                      f"{dict(sorted(lean_spec.items()))} -> {'OK' if ok else 'MISMATCH'}")
    return mismatches


def certificate_cases():
    for seq in SEQUENCES:
        if len(seq) < 4:
            depths = [2, len(seq)]
        else:
            depths = sorted({2, 3, len(seq) // 2 + 1, len(seq)})
        for q in QS:
            for thr in [0, 1, 2, 4]:
                for depth in depths:
                    yield seq, q, thr, depth


def test_certificates(lean: Lean, report):
    texts, meta = [], []
    for seq, q, thr, depth in certificate_cases():
        cert = build_frontier(seq, q, thr, depth)
        py = verify_certificate(cert)
        wit = {}
        for node in cert.frontier:
            if node.kind == "bound":
                wit[node.path] = propose_dfs_witnesses(node.path, len(seq), node_budget=60)
        for mode in ("B", "T", "TW"):
            texts.append(certificate_to_lean_text(cert, transfer=(mode != "B"),
                                                  witnesses=wit if mode == "TW" else None))
            meta.append((seq, q, thr, depth, mode, py))
    rows, wall = lean.check_texts(texts)
    failures = 0
    for (seq, q, thr, depth, mode, py), r in zip(meta, rows):
        z, za, zn, p = exact_values(seq, q, thr)
        if not r.get("accepted"):
            failures += 1
            report.append(f"REJECTED honest cert {seq} q={q} thr={thr} depth={depth} {mode}")
            continue
        lo = {k: frac(r[k]) for k in ["zLo", "zHi", "eLo", "eHi", "nLo", "nHi", "pLo", "pHi"]}
        encl = (lo["zLo"] <= z <= lo["zHi"] and lo["eLo"] <= za <= lo["eHi"]
                and lo["nLo"] <= zn <= lo["nHi"] and lo["pLo"] <= p <= lo["pHi"])
        same_as_python = True
        if mode == "B":
            same_as_python = (lo["zLo"], lo["zHi"], lo["eLo"], lo["eHi"], lo["nLo"], lo["nHi"],
                              lo["pLo"], lo["pHi"]) == (
                py.partition_lower, py.partition_upper, py.event_lower, py.event_upper,
                py.complement_lower, py.complement_upper, py.probability_lower,
                py.probability_upper)
        if not (encl and same_as_python):
            failures += 1
            report.append(f"FAIL {seq} q={q} thr={thr} depth={depth} {mode}: encl={encl} "
                          f"python_match={same_as_python} lean={r}")
    report.append(f"certificates: {len(rows)} Lean checks ({wall:.1f}s wall incl. startup), "
                  f"{failures} failures")
    return failures


def mutations(seq, q, thr, depth):
    cert = build_frontier(seq, q, thr, depth)
    base = certificate_to_lean_text(cert).splitlines()[4:-1]
    n = len(seq)
    out = []
    for i in range(len(base)):
        out.append(("drop node %d" % i, base[:i] + base[i + 1:]))
    out.append(("duplicate node", base + [base[0]]))
    for i, line in enumerate(base):
        k, w = line.split(" ")[:2]
        if k == "B":
            out.append(("bound->exact", base[:i] + [f"E {w}"] + base[i + 1:]))
            out.append(("bound->dead", base[:i] + [f"D {w}"] + base[i + 1:]))
            out.append(("hidden descendant", base + [f"B {w}0" if w != '.' else "B 0"]))
            out.append(("self-intersecting witness",
                        base[:i] + [f"B {w} " + ("1" * (n - 2 - (0 if w == '.' else len(w))))]
                        + base[i + 1:]))
            break
    for i, line in enumerate(base):
        k, w = line.split(" ")[:2]
        if k == "E":
            out.append(("exact->dead", base[:i] + [f"D {w}"] + base[i + 1:]))
            out.append(("reversal step", base[:i] + [f"E {w[:-1]}{'1' if w[-1] == '0' else '0'}"]
                        + base[i + 1:]) if w != '.' else None)
            break
    out.append(("invalid geometry (reverse first step)", base + ["E " + "1" * (n - 2)]))
    return [(name, raw_lean_text(seq, q, thr, nodes)) for name, nodes in
            [m for m in out if m is not None]]


def test_adversarial(lean: Lean, report):
    texts, names = [], []
    for seq, depth in [("HPHPHH", 3), ("HHHHHH", 4), ("HHPHHPH", 3), ("HPPHHPPH", 4)]:
        for name, text in mutations(seq, Fraction(2), 2, depth):
            texts.append(text)
            names.append(f"{seq}/{name}")
    # q < 1 and too-short sequence and duplicate witnesses
    texts.append(raw_lean_text("HHHH", Fraction(1, 2), 1, ["B ."]))
    names.append("q<1")
    texts.append(raw_lean_text("HHHH", Fraction(2), 1, ["B . 00 00"]))
    names.append("duplicate witness")
    texts.append(raw_lean_text("HHHH", Fraction(2), 1, ["B . 0"]))
    names.append("short witness")
    rows, _ = lean.check_texts(texts)
    accepted_bad = 0
    for name, r in zip(names, rows):
        if r.get("accepted"):
            # dropping/adding a node may legitimately still be a valid partition only if the
            # dropped node was a collision dead node; record and verify enclosure separately
            accepted_bad += 1
            report.append(f"ACCEPTED mutated certificate: {name} -> {r}")
    report.append(f"adversarial: {len(rows)} mutated certificates, {accepted_bad} accepted")
    return accepted_bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lake-root", default=str(HERE.parent))
    ap.add_argument("--runner", default="scripts/foldcheck.lean")
    args = ap.parse_args()
    lean = Lean(args.lake_root, args.runner)
    report = []
    bad = test_semantics(lean, report)
    bad += test_certificates(lean, report)
    bad += test_adversarial(lean, report)
    print("\n".join(report))
    print("TOTAL PROBLEMS:", bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
