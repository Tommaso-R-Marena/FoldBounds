#!/usr/bin/env python3
"""Bulk replay of all 1,717 frozen FoldBounds fixed-reference pairs.

Requires numpy and scipy for UNTRUSTED witness generation. All score
enclosures use the standalone integer-outward interval kernel in
motion_fixedpoint_core.py, not the floating-point solver.

Usage:
  python3 scripts/check_frozen_bulk_intervals.py \
    FoldBounds_Prospective_SEND_BACK.zip --output-dir audit

The frozen NPZ arrays contain ZERO-BASED COORDINATE POSITIONS; they are
NOT residue IDs. Every source-byte hash and eligibility pair is verified.
This is reproducible exact-integer Python arithmetic, not a Lean kernel
replay and not evidence that the biological motion model is correct.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from motion_fixedpoint_core import B, D, Y_SCALE, hessian_boxes, bound_score


def sha(data):
    return hashlib.sha256(data).hexdigest()


def frozen_reference(z, label):
    raw = z.read(label + "_FROZEN_BEFORE_TARGET.npz")
    with np.load(io.BytesIO(raw), allow_pickle=False) as d:
        xyz = np.array(d["reference_coords"], dtype=np.float64)
        ids = [int(i) for i in d["residues"]]
        pairs = np.stack((d["pair_i"],d["pair_j"]),axis=1).tolist()
        archive_prediction = np.array(
            d["prediction_ANM_weighted_inverse_square"], dtype=np.float64)
    n = len(ids)
    if xyz.shape != (n, 3) or not np.isfinite(xyz).all():
        raise ValueError("invalid coordinates")
    if len(set(ids)) != n or any(i <= 0 for i in ids):
        raise ValueError("bad residue IDs")
    ci = np.rint(xyz*D).astype(np.int64)
    if not np.array_equal(ci.astype(np.float64)/D, xyz):
        raise ValueError("coordinates not exactly represented at 2**24")
    ci = ci.tolist()
    required = []
    for i in range(1, n-1):
        for j in range(i+4,n-1):
            s = sum((ci[i][k]-ci[j][k])**2 for k in range(3))
            if 25*D*D <= s <= 256*D*D:
                required.append([i,j])
    if pairs != required:
        raise ValueError("frozen pair positions violate the complete 0-based eligibility rule")
    if len(pairs) != len(archive_prediction):
        raise ValueError("archived pair/prediction length mismatch")
    return raw, xyz, ci, ids, pairs, archive_prediction


def float_proposals(xyz, ci, pairs):
    """Untrusted Cholesky proposals constructed from the same rational geometry."""
    n = len(ci)
    H = np.eye(3*n, dtype=float)/5
    for a in range(n):
        for b in range(a+1, n):
            vec = [ci[a][k]-ci[b][k] for k in range(3)]
            s = sum(v*v for v in vec)
            if s == 0:
                raise ValueError("colliding C-alpha coordinates")
            if s > 100*D*D:
                continue
            v = np.array(vec,dtype=float)/D
            K = 49*np.outer(v,v)/(float(s)/D**2)**2
            A = slice(3*a, 3*a+3)
            C = slice(3*b, 3*b+3)
            H[A,A]+=K
            H[C,C]+=K
            H[A,C]-=K
            H[C,A]-=K
    rhs = np.zeros((3*n,len(pairs)),dtype=float)
    for t,(a,b) in enumerate(pairs):
        v = xyz[a]-xyz[b]
        rhs[3*a:3*a+3,t]=v
        rhs[3*b:3*b+3,t]=-v
    solve = cho_solve(cho_factor(H,lower=True),rhs)
    return rhs,solve


def check_one(z,label,outdir):
    raw,xyz,ci,ids,pairs,archived = frozen_reference(z,label)
    intervals, active = hessian_boxes(ci)
    rhs,solutions = float_proposals(xyz,ci,pairs)
    bounds=[]
    proposed=[]
    for t,pair in enumerate(pairs):
        proposal=[int(round(float(v)*Y_SCALE)) for v in solutions[:,t]]
        lo,hi=bound_score(intervals,pair,ci,proposal)
        bounds.append((lo,hi))
        proposed.append(proposal)
    scores = np.einsum("ij,ij->j",rhs,solutions)
    dsq = np.array([sum((xyz[i][k]-xyz[j][k])**2
                         for k in range(3)) for i,j in pairs])
    scores /= dsq
    err = float(np.max(np.abs(scores-archived)))
    if err > 5e-10:
        raise ValueError(f"archived score mismatch: {label}, {err}")
    order = sorted(range(len(pairs)),key=lambda t:archived[t])
    separated = sum(bounds[order[t]][1] < bounds[order[t+1]][0]
                    for t in range(len(order)-1))
    records=[{
        "position_i":i,"position_j":j,
        "residue_id_i":ids[i],"residue_id_j":ids[j],
        "lower_scaled_2pow96":str(bounds[t][0]),
        "upper_scaled_2pow96":str(bounds[t][1])
    } for t,(i,j) in enumerate(pairs)]
    interval_file=outdir / (label+"_Interval_Certificates.json")
    interval_file.write_text(json.dumps({
        "format":"FOLDBULK-INTERVAL-1","protein":label,
        "frozen_npz_sha256":sha(raw),
        "interval_denominator":"2**96",
        "proposal_denominator":"2**30","records":records},
        separators=(",",":"),sort_keys=True)+"\n",encoding="utf8")
    np.savez_compressed(outdir/(label+"_Untrusted_Proposals.npz"),
        proposal_y_numerators=np.asarray(proposed,dtype=np.int64),
        zero_based_pair_positions=np.asarray(pairs,dtype=np.int32),
        denominator=np.asarray([Y_SCALE],dtype=np.int64))
    widths=[(hi-lo)/B for lo,hi in bounds]
    return {
        "protein":label,
        "residues":len(ids),
        "frozen_npz_sha256":sha(raw),
        "eligible_pairs":len(pairs),
        "active_contacts":active,
        "intervals_completed":len(bounds),
        "strict_adjacent_ranks":separated,
        "adjacent_rank_total":len(pairs)-1,
        "max_interval_width":max(widths),
        "maximum_float_vs_archived_difference":err
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("archive")
    p.add_argument("--output-dir",default="motion_bulk_audit")
    args=p.parse_args()
    folder=Path(args.output_dir)
    folder.mkdir(parents=True,exist_ok=True)
    archive=Path(args.archive)
    with zipfile.ZipFile(archive) as z:
        names=z.namelist()
        if len(names)!=len(set(names)) or any(
                name.startswith("/") or ".." in name.split("/")
                or "\\" in name for name in names):
            raise ValueError("unsafe archive members")
        hashes=json.loads(z.read("SHA256SUMS.json"))
        if set(hashes)!=set(names)-{"SHA256SUMS.json"}:
            raise ValueError("incomplete SHA256 manifest")
        for name,expected in hashes.items():
            if sha(z.read(name)) != expected:
                raise ValueError("SHA256 mismatch "+name)
        frozen=json.loads(z.read("ALL_PREDICTION_SHA256_BEFORE_NMR.json"))
        protocol=json.loads(z.read("FROZEN_PROTOCOL.json"))
        if (protocol["contact_A"],protocol["ridge"],
                protocol["distance_scale_A"])!=(10.0,0.2,7.0):
            raise ValueError("unexpected model parameters")
        for label in ("CspB","Crambin"):
            fn=label+"_FROZEN_BEFORE_TARGET.npz"
            if frozen[label]!=hashes[fn]:
                raise ValueError("invalid frozen-reference hash")
        results=[check_one(z,label,folder) for label in ("CspB","Crambin")]
    summary={"archive_sha256":sha(archive.read_bytes()),
             "trust_boundary":"Python exact outward integer intervals; NOT Lean replay",
             "result":results}
    (folder/"SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
