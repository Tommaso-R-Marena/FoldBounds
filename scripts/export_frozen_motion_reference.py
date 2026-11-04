#!/usr/bin/env python3
"""Fail-closed stdlib exporter for the frozen FoldBounds protein references.

Converts archived IEEE-754 coordinates into exact rational numerator/denominator
pairs, checks every archive SHA256 and frozen-prediction hash, and never reads NMR
evaluation arrays into the certified model input.

This is an independent Python input adapter, NOT a Lean-verified decoder.
"""
import argparse
import ast
import hashlib
import io
import json
import math
import struct
import sys
import zipfile

MAX_MEMBER = 20_000_000


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_npy(blob, dtype):
    if len(blob) < 12 or blob[:6] != b"\x93NUMPY":
        raise ValueError("bad NPY magic")
    major, minor = blob[6:8]
    if (major, minor) == (1, 0):
        offset, nhead = 10, struct.unpack("<H", blob[8:10])[0]
    elif (major, minor) in ((2, 0), (3, 0)):
        offset, nhead = 12, struct.unpack("<I", blob[8:12])[0]
    else:
        raise ValueError("unsupported NPY version")
    if nhead > 16000:
        raise ValueError("excessive NPY header")
    meta = ast.literal_eval(blob[offset:offset + nhead].decode("latin1"))
    if meta.get("descr") != dtype or meta.get("fortran_order") is not False:
        raise ValueError("wrong dtype or order")
    shape = meta.get("shape")
    if (not isinstance(shape, tuple) or len(shape) not in (1, 2)
            or not all(isinstance(k, int) and 0 < k < 100000 for k in shape)):
        raise ValueError("invalid shape")
    count = math.prod(shape)
    typ = {"<f8": "d", "<i4": "i"}[dtype]
    payload = blob[offset + nhead:]
    if len(payload) != count * struct.calcsize(typ):
        raise ValueError("wrong data length")
    return shape, struct.unpack("<" + str(count) + typ, payload)


def parse_reference_npz(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names = z.namelist()
        if len(names) != len(set(names)) or len(names) > 30:
            raise ValueError("duplicate or excessive NPZ entries")
        wanted = {"reference_coords.npy", "residues.npy", "pair_i.npy", "pair_j.npy"}
        if not wanted <= set(names):
            raise ValueError("missing reference array")
        for name in wanted:
            if z.getinfo(name).file_size > MAX_MEMBER:
                raise ValueError("oversize array")
        shape, coords = read_npy(z.read("reference_coords.npy"), "<f8")
        ids_shape, ids = read_npy(z.read("residues.npy"), "<i4")
        pair_shape, pi = read_npy(z.read("pair_i.npy"), "<i4")
        other_shape, pj = read_npy(z.read("pair_j.npy"), "<i4")
    n = shape[0]
    if shape != (n, 3) or ids_shape != (n,) or pair_shape != other_shape:
        raise ValueError("inconsistent shapes")
    if not all(math.isfinite(value) for value in coords):
        raise ValueError("nonfinite coordinates")
    ids_set = set(ids)
    pairs = list(zip(pi, pj))
    if len(ids_set) != n or len(set(pairs)) != len(pairs) or any(
            i >= j or i not in ids_set or j not in ids_set for i, j in pairs):
        raise ValueError("invalid residue indices or duplicate pair")
    rational_coords = []
    for i in range(n):
        rational_coords.append([
            list(coords[3*i + axis].as_integer_ratio()) for axis in range(3)])
    return {"residue_ids": list(ids),
            "coordinate_rationals": rational_coords,
            "evaluation_pair_residue_ids": [[i, j] for i, j in pairs]}


def export_archive(path):
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        if len(set(names)) != len(names) or any(
                not name or name.startswith("/") or "\\" in name
                or ".." in name.split("/") for name in names):
            raise ValueError("unsafe outer archive members")
        if any(item.file_size > MAX_MEMBER for item in z.infolist()):
            raise ValueError("oversize outer archive member")
        sums = json.loads(z.read("SHA256SUMS.json"))
        if set(sums) != set(names) - {"SHA256SUMS.json"}:
            raise ValueError("incomplete SHA256 manifest")
        for name, claimed in sums.items():
            if sha(z.read(name)) != claimed:
                raise ValueError("SHA256 mismatch " + name)
        frozen = json.loads(z.read("ALL_PREDICTION_SHA256_BEFORE_NMR.json"))
        protocol = json.loads(z.read("FROZEN_PROTOCOL.json"))
        if (protocol.get("contact_A") != 10.0
                or protocol.get("ridge") != 0.2
                or protocol.get("distance_scale_A") != 7.0):
            raise ValueError("unexpected frozen protocol")
        result = {
            "format": "FOLDREF-RAT-1",
            "archive_sha256": sha(open(path, "rb").read()),
            "model": {"contact_squared_A2": [100, 1],
                      "ridge": [1, 5],
                      "distance_scale_squared_A2": [49, 1]},
            "inputs": {}}
        for protein in ("CspB", "Crambin"):
            name = protein + "_FROZEN_BEFORE_TARGET.npz"
            if frozen[protein] != sums[name]:
                raise ValueError("unmatched frozen hash " + protein)
            out = parse_reference_npz(z.read(name))
            out["frozen_npz_sha256"] = sums[name]
            result["inputs"][protein] = out
        return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive")
    parser.add_argument("--output")
    args = parser.parse_args()
    result = export_archive(args.archive)
    txt = json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="\n") as f:
            f.write(txt)
    else:
        sys.stdout.write(txt)
    for name, values in result["inputs"].items():
        print(name, len(values["residue_ids"]), "residues,",
              len(values["evaluation_pair_residue_ids"]), "pairs", file=sys.stderr)


if __name__ == "__main__":
    main()
