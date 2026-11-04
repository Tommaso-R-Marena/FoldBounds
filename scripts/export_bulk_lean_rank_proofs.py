#!/usr/bin/env python3
"""Rebuild a complete Lean integer-rank replay from a frozen FoldBounds audit ZIP.

Exact endpoint ordering ONLY. The score intervals must be separately checked
with MotionScaledCertificates to establish actual physical-model score bounds.
The ZIP and member SHA-256 hashes and corrected position/ID mappings are
checked, but these byte checks are Python, not Lean-kernel claims.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile

PREFIX = "FoldBounds_Bulk_Rational_Audit_20261104/"
EXPECTED = {"CspB": 1170, "Crambin": 547}


def sha(blob):
    return hashlib.sha256(blob).hexdigest()


def load_audit(path):
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if not names or len(names) != len(set(names)):
            raise ValueError("duplicate/empty ZIP members")
        if any(not n.startswith(PREFIX) or "\\" in n or ".." in n.split("/")
               for n in names):
            raise ValueError("unsafe ZIP name")
        if any(i.file_size > 15_000_000 for i in archive.infolist()):
            raise ValueError("oversized ZIP member")
        man = PREFIX + "SHA256SUMS.json"
        hashes = json.loads(archive.read(man))
        if set(hashes) != {n[len(PREFIX):] for n in names if n != man}:
            raise ValueError("incomplete SHA256 manifest")
        for name, expected in hashes.items():
            if not re.fullmatch("[a-f0-9]{64}", expected):
                raise ValueError("bad digest")
            if sha(archive.read(PREFIX + name)) != expected:
                raise ValueError("archive SHA256 mismatch: " + name)
        inp = json.loads(archive.read(
            PREFIX + "FoldBounds_Frozen_Exact_Rational_Inputs_CORRECTED_V2.json"))
        if inp.get("format") != "FOLDREF-RAT-2":
            raise ValueError("uncorrected reference identity format")
        result = {}
        for protein, count in EXPECTED.items():
            filename = f"FoldBounds_Real_Bulk_{protein}_Certificate_Intervals.json"
            cert = json.loads(archive.read(PREFIX + filename))
            refs = inp["inputs"][protein]
            if (cert.get("format") != "FOLDBULK-INTERVAL-1"
                or cert.get("protein") != protein
                or cert.get("frozen_npz_sha256") != refs["frozen_npz_sha256"]
                or cert.get("fixed_point_denominator") != "2**96"
                or cert.get("proposal_denominator") != "2**30"):
                raise ValueError("incompatible fixed-point data")
            entries = cert["scores"]
            positions = refs["evaluation_pair_positions"]
            ids = refs["evaluation_pair_residue_ids"]
            residues = refs["residue_ids"]
            if len(entries) != count or len(positions) != count:
                raise ValueError("wrong number of pairs")
            checked = []
            seen = set()
            for t, (entry, pair, rid) in enumerate(zip(entries, positions, ids)):
                i, j = pair
                if (not 0 <= i < j < len(residues) or tuple(pair) in seen
                    or [residues[i], residues[j]] != rid
                    or [entry["position_i"], entry["position_j"]] != pair
                    or [entry["residue_id_i"], entry["residue_id_j"]] != rid):
                    raise ValueError("incorrect position/ID mapping")
                seen.add(tuple(pair))
                lo = int(entry["lower_scaled_2pow96"])
                hi = int(entry["upper_scaled_2pow96"])
                if lo > hi:
                    raise ValueError("reversed interval")
                checked.append((lo, hi, t))
            checked.sort()
            if any(checked[i][1] >= checked[i+1][0]
                   for i in range(count-1)):
                raise ValueError("overlap or tie in strict adjacent ranking")
            result[protein] = (checked, hashes[filename])
    return result, sha(Path(path).read_bytes())


def overlapping_blocks(entries, width):
    if width < 2:
        raise ValueError("need at least two intervals per block")
    start = 0
    while start + 1 < len(entries):
        end = min(start + width, len(entries))
        yield entries[start:end]
        start = end - 1


def generate_source(proteins, archive_sha, block_size=32):
    lines = [
        "import FoldBounds.MotionBulkChain", "",
        "/-! Bulk frozen dyadic rank arithmetic.",
        "Generated from SHA256-validated reference source: " + archive_sha,
        "Each theorem checks a block of PROPOSED intervals, not",
        "the correctness of all physical-model score endpoints. -/", "",
        "open FoldBounds.MotionScaledCertificates",
        "open FoldBounds.MotionBulkChain",
        "namespace FoldBounds.FrozenRankReplay", ""
    ]
    stats = {}
    for protein, (entries, source_sha) in proteins.items():
        blocks = list(overlapping_blocks(entries, block_size))
        assert sum(len(b) - 1 for b in blocks) == len(entries) - 1
        stats[protein] = {
            "score_intervals": len(entries),
            "adjacent_comparisons": len(entries)-1,
            "blocks": len(blocks),
            "certificate_file_sha256": source_sha,
            "min_integer_gap": str(min(entries[i+1][0]-entries[i][1]
                                       for i in range(len(entries)-1)))
        }
        for k, block in enumerate(blocks):
            tag = f"{protein}_block_{k:03d}"
            lines.append(f"private def {tag} : List ScaledInterval := [")
            for i, (lo, hi, original_pair_index) in enumerate(block):
                comma = "," if i + 1 != len(block) else ""
                lines.append(f"  ⟨{lo}, {hi}⟩{comma} -- frozen pair #{original_pair_index}")
            lines += ["]", f"theorem {tag}_ok : checkAdjacent {tag} = true := by native_decide",
                      f"theorem {tag}_strict : StrictAdjacent {tag} :=",
                      f"  checkAdjacent_sound {tag} {tag}_ok", ""]
    lines += ["end FoldBounds.FrozenRankReplay", ""]
    return "\n".join(lines), stats


def main():
    p = argparse.ArgumentParser()
    p.add_argument("archive", type=Path)
    p.add_argument("--output", type=Path, default=Path("MotionFrozenRanksGenerated.lean"))
    p.add_argument("--manifest", type=Path,
                   default=Path("MotionFrozenRanksGenerated.manifest.json"))
    p.add_argument("--block-size", type=int, default=32)
    args = p.parse_args()
    data, archive_sha = load_audit(args.archive)
    content, stats = generate_source(data, archive_sha, args.block_size)
    args.output.write_text(content, encoding="utf-8")
    metadata = {"format": "FOLDBULK-RANK-LEAN-1",
                "outer_archive_sha256": archive_sha,
                "generated_lean_sha256": sha(content.encode("utf-8")),
                "chunk_size": args.block_size, "proteins": stats}
    args.manifest.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
