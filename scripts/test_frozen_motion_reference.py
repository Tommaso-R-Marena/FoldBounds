#!/usr/bin/env python3
"""Synthetic, adversarial tests for the independent frozen NPZ byte adapter."""
import io
import json
import os
import struct
import tempfile
import zipfile
from export_frozen_motion_reference import export_archive, sha


def npy(shape, dtype, values):
    magic = b"\x93NUMPY\x01\x00"
    header = str({"descr": dtype, "fortran_order": False, "shape": shape})
    n = len(magic) + 2 + len(header) + 1
    header += " " * ((-n) % 16) + "\n"
    hdr = magic + struct.pack("<H", len(header)) + header.encode("ascii")
    code = {"<f8": "d", "<i4": "i"}[dtype]
    return hdr + struct.pack("<" + str(len(values)) + code, *values)


def frozen_npz(make_nan=False, bad_pair=False):
    # Seven residues, noncontiguous PDB IDs and one eligible 0-based pair.
    # This catches the former off-by-one position-versus-residue-ID bug.
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as z:
        xs = [1000.0, 0.0, 200.0, 300.0, 400.0, 8.0, 500.0]
        vals = [v for xx in xs for v in (xx, 0.0, 0.0)]
        if make_nan:
            vals[16] = float("nan")
        z.writestr("reference_coords.npy", npy((7, 3), "<f8", vals))
        z.writestr("residues.npy",
                   npy((7,), "<i4", [11, 17, 20, 23, 29, 35, 42]))
        z.writestr("pair_i.npy", npy((1,), "<i4", [1]))
        z.writestr("pair_j.npy", npy((1,), "<i4", [4 if bad_pair else 5]))
    return buffer.getvalue()


def archive(path, make_nan=False, tamper=False, bad_pair=False):
    a = frozen_npz(make_nan, bad_pair)
    b = frozen_npz()
    files = {
        "CspB_FROZEN_BEFORE_TARGET.npz": a,
        "Crambin_FROZEN_BEFORE_TARGET.npz": b,
        "FROZEN_PROTOCOL.json": json.dumps(
            {"contact_A": 10, "ridge": 0.2, "distance_scale_A": 7}).encode(),
        "ALL_PREDICTION_SHA256_BEFORE_NMR.json":
            json.dumps({"CspB": sha(a), "Crambin": sha(b)}).encode(),
    }
    sums = {k: sha(v) for k, v in files.items()}
    if tamper:
        files["FROZEN_PROTOCOL.json"] += b" "
    files["SHA256SUMS.json"] = json.dumps(sums).encode()
    with zipfile.ZipFile(path, "w") as z:
        for key, val in files.items():
            z.writestr(key, val)


def expect_bad(path, **kwargs):
    archive(path, **kwargs)
    try:
        export_archive(path)
    except ValueError:
        return
    raise AssertionError("unsafe input accepted")


def main():
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "fixture.zip")
        archive(path)
        result = export_archive(path)
        assert result["format"] == "FOLDREF-RAT-2"
        assert result["inputs"]["CspB"]["coordinate_rationals"][5][0] == [8, 1]
        assert result["inputs"]["CspB"]["evaluation_pair_positions"] == [[1, 5]]
        assert result["inputs"]["CspB"]["evaluation_pair_residue_ids"] == [[17, 35]]
        assert result["inputs"]["Crambin"]["evaluation_pair_residue_ids"] == [[17, 35]]
        expect_bad(path, tamper=True)
        expect_bad(path, make_nan=True)
        expect_bad(path, bad_pair=True)
    print("PASS: frozen reference exact-rational parser, hash-tamper rejection and NaN rejection")


if __name__ == "__main__":
    main()
