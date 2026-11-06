#!/usr/bin/env python3
"""Independent exact-integer replay of source-only GNM precision comparison.

Does not import the ablation producer. Recomputes the untrusted floating
inverse; each response interval is then checked with Python unbounded integers.
The exact spectral SOS lower bound is independently verified before use.

Run with the original SHA-pinned reference-only Stage1 archive:
 OPENBLAS_NUM_THREADS=1 python scripts/gnm_precision_independent_replay.py \
   STAGE1.zip spectral.json research/GNMSourceOnlyPrecisionAblation_20261105.json
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from pathlib import Path
import json
import zipfile
import io
import hashlib
import numpy as np

from gnm_spectral_sos_independent_replay import (
    verify as verify_spectral,
    frozen_source,
    certain_edges_independent,
    grounded_laplacian_integer,
)

EXPECTED_SOURCE = "ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff"
Q_BITS = (10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 44, 48, 52, 55)


def source_pairs(path):
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != EXPECTED_SOURCE:
        raise ValueError("Wrong frozen source archive")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        manifest = json.loads(archive.read("SHA256SUMS.json"))
        if (archive.testzip() is not None or
            set(manifest) != set(archive.namelist()) - {"SHA256SUMS.json"} or
            any(hashlib.sha256(archive.read(k)).hexdigest() != digest
                for k, digest in manifest.items())):
            raise ValueError("Frozen reference archive manifest mismatch")
        with np.load(io.BytesIO(archive.read("FROZEN_REFERENCE_ONLY.npz")),
                     allow_pickle=False) as data:
            return list(zip([int(x) for x in data["pair_i"]],
                            [int(x) for x in data["pair_j"]]))


def residual(A, N, Q):
    n = len(A)
    sparse = [[(j, x) for j, x in enumerate(row) if x] for row in A]
    return [[int(Q * (i == j) -
                 sum(a * N[k][j] for k, a in sparse[i])) for j in range(n)]
            for i in range(n)]


def intervals(A, N, R, Q, B_num, B_den, pairs, root):
    n = len(A)
    idx = {k: (k if k < root else k - 1)
           for k in range(n + 1) if k != root}
    bounds = []
    for i, j in pairs:
        ai, bj = idx.get(i), idx.get(j)
        # The pair force has only two nonzero entries; column differences
        # avoid dense O(n^2) matrix/vector multiplication.
        y = [(N[k][ai] if ai is not None else 0) -
             (N[k][bj] if bj is not None else 0) for k in range(n)]
        r = [(R[k][ai] if ai is not None else 0) -
             (R[k][bj] if bj is not None else 0) for k in range(n)]
        b_y = (y[ai] if ai is not None else 0) - (y[bj] if bj is not None else 0)
        lower = Q * b_y + sum(yy * rr for yy, rr in zip(y, r))
        bounds.append((lower * B_den,
                       lower * B_den + B_num * sum(rr * rr for rr in r)))
    order = sorted(range(len(bounds)), key=lambda k: (bounds[k][0], k))
    count = sum(bounds[order[i]][1] < bounds[order[i + 1]][0]
                for i in range(len(order) - 1))
    width = max(hi - lo for lo, hi in bounds)
    return count, Fraction(width, B_den * Q * Q)


def verify(stage1, spectral, ablation, selected_bits=None):
    report = verify_spectral(stage1, spectral)
    if report["B_upper"] != "1000000/169821":
        raise ValueError("Wrong spectral lower bound")
    coords, _ = frozen_source(stage1)
    certain, optional = certain_edges_independent(coords)
    if (len(certain), len(optional)) != (972, 108):
        raise ValueError("Unexpected source topology")
    root = 45
    A = grounded_laplacian_integer(certain, 124, root)
    pairs = source_pairs(stage1)
    if (len(pairs) != 2614 or len(set(pairs)) != 2614 or
        any(i < 0 or i >= j or j >= 124 for i, j in pairs)):
        raise ValueError("Not 2614 canonical frozen pairs")
    rows_doc = json.loads(Path(ablation).read_text())
    if (rows_doc.get("format") != "GNM_SOURCE_ONLY_PRECISION_ABLATION_V1" or
        rows_doc.get("source_sha256") != EXPECTED_SOURCE):
        raise ValueError("Wrong ablation report identity")
    rows = rows_doc["results"]
    if [row["q_bits"] for row in rows] != list(Q_BITS):
        raise ValueError("Unexpected witness precisions")
    required = {"q_bits", "old_adjacent_orders", "new_adjacent_orders",
                "old_max_score_width", "new_max_score_width"}
    if any(set(row) != required or not
           (0 <= row["old_adjacent_orders"] <= row["new_adjacent_orders"] <= 2613)
           for row in rows):
        raise ValueError("Malformed precision report or violated monotonicity")
    inverse = np.linalg.inv(np.asarray(A, dtype=np.float64))  # untrusted proposal
    out = []
    for row in rows:
        bits = row["q_bits"]
        if selected_bits is not None and bits not in selected_bits:
            continue
        Q = 1 << bits
        N = [[int(round(float(x) * Q)) for x in line] for line in inverse]
        R = residual(A, N, Q)
        old_n, old_w = intervals(A, N, R, Q, 262, 1, pairs, root)
        new_n, new_w = intervals(A, N, R, Q, 1000000, 169821, pairs, root)
        if (old_n, new_n) != (row["old_adjacent_orders"],
                              row["new_adjacent_orders"]):
            raise AssertionError("Published adjacent counts not reproduced at Q=2^" +
                                 str(bits))
        for field, width in (("old_max_score_width", old_w),
                             ("new_max_score_width", new_w)):
            expected = row[field]
            if abs(float(width) - expected) > 1e-13 * max(abs(expected), 1e-25):
                raise AssertionError("Published width not reproduced: " +
                                     field + ", Q=2^" + str(bits))
        out.append((bits, old_n, new_n))
    print("INDEPENDENT PRECISION ABLATION VERIFIED:", out)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage1")
    ap.add_argument("spectral")
    ap.add_argument("ablation")
    ap.add_argument("--bits", type=int, nargs="+", default=None)
    x = ap.parse_args()
    verify(x.stage1, x.spectral, x.ablation, x.bits)
