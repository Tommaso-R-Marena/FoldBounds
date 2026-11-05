#!/usr/bin/env python3
"""Exact integer residual certificates for connected, unweighted GNM networks.

This verifies the *specified* GNM effective resistance without a floating-point
trust assumption.  It is not a proof that an NMR ensemble is a physical dynamic
process, and it is not an end-to-end Lean certificate for parsed NPZ files.

For a connected contact graph and a grounded Laplacian A, with anchor vertex a,
we use a reference-only BFS tree.  Let B = sum_v distance_T(v,a).
Cauchy--Schwarz on each tree path implies ||u||² <= B*u'Au for ALL grounded u.
For proposed integer N approximating Q*A^{-1}, define r=Q*b-A*(N*b), u=N*b.
Then with exact integer numerators:
    l=(Q*b'u+u'r)/Q²,
    h=(Q*b'u+u'r + B*r'r)/Q²,
we have l <= b'A^{-1}b <= h for every target pair b=e_i-e_j.
The proposed matrix N is completely untrusted; all arithmetic after reading
N is exact integer arithmetic. Graph threshold is checked with exact binary64
rational coordinates; no rounded Euclidean distance decides contact status.

Computational detail: 0<distance²<=cutoff²; threshold is inclusive.  Unrelated
pairs at exact coincident coordinates cause fail-closed. An approximate matrix
inverse requires floating point only for the *proposal*, not the verifier.
"""
from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile
import numpy as np

STAGE1_SHA = "ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff"
SOURCE_NAME = "FROZEN_REFERENCE_ONLY.npz"
Q_DEFAULT = 1 << 55


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def open_sealed_stage1(path: str | Path, expected_sha: str = STAGE1_SHA):
    data = Path(path).read_bytes()
    if sha256(data) != expected_sha:
        raise ValueError("ABSTAIN: Stage1 byte SHA256 mismatch")
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        if z.testzip() is not None:
            raise ValueError("ABSTAIN: ZIP CRC failure")
        manifest = json.loads(z.read("SHA256SUMS.json"))
        if set(manifest) != set(z.namelist()) - {"SHA256SUMS.json"}:
            raise ValueError("ABSTAIN: incomplete ZIP SHA manifest")
        for name, digest in manifest.items():
            if sha256(z.read(name)) != digest:
                raise ValueError("ABSTAIN: internal SHA256 mismatch: " + name)
        npz = z.read(SOURCE_NAME)
    with np.load(io.BytesIO(npz), allow_pickle=False) as a:
        return {key: a[key].copy() for key in a.files}, sha256(npz)


def _validate_inputs(coords, pairs_i, pairs_j):
    x = np.asarray(coords, dtype=np.float64)
    i = np.asarray(pairs_i)
    j = np.asarray(pairs_j)
    if x.ndim != 2 or x.shape[1] != 3 or x.shape[0] < 3:
        raise ValueError("ABSTAIN: invalid C-alpha coordinate array")
    if not np.isfinite(x).all() or not np.issubdtype(i.dtype, np.integer) or not np.issubdtype(j.dtype, np.integer):
        raise ValueError("ABSTAIN: invalid coordinates or pair index types")
    if i.ndim != 1 or j.shape != i.shape or len(i) == 0:
        raise ValueError("ABSTAIN: invalid pair index arrays")
    if np.any(i < 0) or np.any(j >= len(x)) or np.any(i >= j):
        raise ValueError("ABSTAIN: pair indices not canonical 0 <= i < j < n")
    pairs = list(zip(map(int,i),map(int,j)))
    if len(pairs) != len(set(pairs)):
        raise ValueError("ABSTAIN: duplicate requested score pair")
    return x, pairs


def exact_contact_graph(coords, cutoff_sq=Fraction(100)):
    x = np.asarray(coords, dtype=np.float64)
    if not np.isfinite(x).all():
        raise ValueError("ABSTAIN: NaN or infinite coordinate")
    if cutoff_sq <= 0:
        raise ValueError("ABSTAIN: nonpositive cutoff")
    # Binary64 values are interpreted as exactly their stored dyadic rationals.
    xyz = [[Fraction.from_float(float(v)) for v in row] for row in x]
    n = len(x)
    neighbors = [set() for _ in range(n)]
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            d2 = sum((xyz[i][k] - xyz[j][k])**2 for k in range(3))
            if not d2:
                raise ValueError("ABSTAIN: coincident C-alpha coordinates")
            if d2 <= cutoff_sq:
                edges.append((i,j))
                neighbors[i].add(j); neighbors[j].add(i)
    return edges, neighbors


def choose_tree(neighbors):
    n = len(neighbors)
    def bfs(root):
        depths = [-1]*n
        parent = [-1]*n
        depths[root] = 0
        queue = deque([root])
        while queue:
            i = queue.popleft()
            for j in sorted(neighbors[i]):
                if depths[j] == -1:
                    depths[j] = depths[i] + 1
                    parent[j] = i
                    queue.append(j)
        if min(depths) < 0:
            raise ValueError("ABSTAIN: disconnected GNM graph")
        return depths, parent
    # Fixed, reference-only anchor with smallest total tree depth.
    best = min(((sum(bfs(root)[0]), root) for root in range(n)))
    B, anchor = best
    d,p = bfs(anchor)
    if sum(d) != B or B <= 0:
        raise AssertionError("Bad BFS witness")
    return anchor, B, d, p


def grounded_laplacian(edges, n, anchor):
    remaining=[i for i in range(n) if i != anchor]
    idx={v:i for i,v in enumerate(remaining)}
    A=[[0]*(n-1) for _ in range(n-1)]
    for a,b in edges:
        if a != anchor: A[idx[a]][idx[a]] += 1
        if b != anchor: A[idx[b]][idx[b]] += 1
        if a != anchor and b != anchor:
            A[idx[a]][idx[b]] -= 1
            A[idx[b]][idx[a]] -= 1
    return A,idx


def untrusted_inverse_proposal(A, Q=Q_DEFAULT):
    # Numerical step, deliberately outside the TCB.
    numeric_inverse=np.linalg.inv(np.asarray(A,dtype=float))
    if not np.isfinite(numeric_inverse).all():
        raise ValueError("ABSTAIN: inverse proposal invalid")
    return [[int(round(float(v)*Q)) for v in row] for row in numeric_inverse]


def verify_integer_witness(A, N, Q):
    m=len(A)
    if not isinstance(Q,int) or isinstance(Q,bool) or not (1 <= Q <= 1<<128):
        raise ValueError("ABSTAIN: invalid exact rational denominator")
    if len(N)!=m or any(len(row)!=m for row in N):
        raise ValueError("ABSTAIN: wrong inverse witness shape")
    if any(type(z) is not int or z.bit_length() > 256 for row in N for z in row):
        raise ValueError("ABSTAIN: nonintegral or huge inverse witness")
    # Recompute R=Q I-A*N using only Python unbounded integers.
    R=[]
    for i in range(m):
        sparse=[(k,int(v)) for k,v in enumerate(A[i]) if v]
        R.append([Q*int(i==j)-sum(v*N[k][j] for k,v in sparse)
                  for j in range(m)])
    return R


def score_numerators(pair, N, R, Q, B, idx):
    i,j=pair
    mi,mj=idx.get(i),idx.get(j)
    m=len(N)
    v=[(N[k][mi] if mi is not None else 0) - (N[k][mj] if mj is not None else 0)
       for k in range(m)]
    r=[(R[k][mi] if mi is not None else 0) - (R[k][mj] if mj is not None else 0)
       for k in range(m)]
    bv=(v[mi] if mi is not None else 0)-(v[mj] if mj is not None else 0)
    lower=Q*bv+sum(vi*ri for vi,ri in zip(v,r))
    upper=lower+B*sum(ri*ri for ri in r)
    if upper < lower:
        raise AssertionError("Invalid bound orientation")
    return lower, upper


def verify_and_score(coords, pairs_i, pairs_j, Q, N, cutoff_sq=Fraction(100)):
    x,pairs = _validate_inputs(coords,pairs_i,pairs_j)
    edges,neighbors=exact_contact_graph(x, cutoff_sq)
    anchor,B,depths,parent=choose_tree(neighbors)
    A,idx=grounded_laplacian(edges,len(x),anchor)
    R=verify_integer_witness(A,N,Q)
    bounds=[score_numerators(pair,N,R,Q,B,idx) for pair in pairs]
    # Global exact score order cert: if all successive intervals are separated,
    # every pairwise rank is fixed in every realization of the model.
    order=sorted(range(len(bounds)),key=lambda k:(bounds[k][0],k))
    strict_adjacent=sum(bounds[order[i]][1]<bounds[order[i+1]][0]
                        for i in range(len(bounds)-1))
    score_mid=np.array([((lo+hi)/2)/(Q*Q) for lo,hi in bounds])
    widths=[hi-lo for lo,hi in bounds]
    return {
        "score_numerators":bounds,"ordered_indices":order,
        "strict_adjacent_count":strict_adjacent,"max_width_num":max(widths),
        "anchor":anchor,"B":B,"depths":depths,"parent":parent,
        "contacts":len(edges),"edges":edges,"score_midpoints":score_mid,
        "denominator":Q*Q,
        "max_residual_integer":max(abs(v) for row in R for v in row),
    }


def make_certificate(stage1_path, output_path, Q=Q_DEFAULT):
    array, npz_sha=open_sealed_stage1(stage1_path)
    coords=array["reference_coords"]
    pairs_i=array["pair_i"];pairs_j=array["pair_j"]
    # Choosing the anchor and preparing the matrix never reads any target.
    edges,neighbors=exact_contact_graph(coords)
    anchor,B,depths,parent=choose_tree(neighbors)
    A,idx=grounded_laplacian(edges,len(coords),anchor)
    N=untrusted_inverse_proposal(A,Q)
    checks=verify_and_score(coords,pairs_i,pairs_j,Q,N)
    # Numerical equivalence to the frozen original GNM, not a mathematical
    # proof about the old NumPy eigensolver.
    target=array["prediction_GNM_resistance"]
    if np.max(np.abs(checks["score_midpoints"]-target)) > 1e-11:
        raise ValueError("ABSTAIN: exact graph model doesn't reproduce frozen GNM")
    certificate={
        "format":"FOLDBOUNDS_GNM_INTEGER_RESIDUAL_V1",
        "source_archive_sha256":STAGE1_SHA,"source_npz_sha256":npz_sha,
        "network":"all unweighted undirected C-alpha pairs with exact dyadic separation squared <= 100",
        "pair_index_base":0,"n_residues":len(coords),"n_pairs":len(pairs_i),
        "anchor":anchor,"tree_parent":parent,"tree_depth":depths,
        "sum_tree_depths_B":B,"denominator_Q":Q,
        "integer_inverse_proposal":N,
        "integer_witness_sha256":sha256(json.dumps(N,separators=(',',':')).encode("ascii")),
        "score_numerators":checks["score_numerators"],
        "ordered_indices":checks["ordered_indices"],
        "strict_adjacent_count":checks["strict_adjacent_count"],
        "n_contacts":checks["contacts"],
        "max_residual_integer":checks["max_residual_integer"],
        "limitations":[
            "Exact Python integer arithmetic; no Lean replay of this dataset yet",
            "Structural model certificate, not physically calibrated NMR prediction",
            "Stored binary64 C-alpha coordinates treated as exact dyadic rationals",
            "Untrusted floating inverse only proposes an integer witness",
        ],
    }
    path=Path(output_path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(certificate,sort_keys=True,separators=(',',':'))+'\n',encoding='utf8')
    print(json.dumps({"status":"PASS","archive":str(path),"sha256":sha256(path.read_bytes()),
                      "n":len(coords),"contacts":checks["contacts"],
                      "pairs":len(pairs_i),"adjacent_strict":checks["strict_adjacent_count"],
                      "adjacent_total":len(pairs_i)-1,"anchor":anchor,
                      "B":B,"Q":Q,"max_float_discrepancy":float(np.max(abs(checks["score_midpoints"]-target)))},indent=2))
    return certificate


def verify_certificate(stage1_path, certificate_path):
    array,npz_sha=open_sealed_stage1(stage1_path)
    raw=Path(certificate_path).read_bytes()
    c=json.loads(raw)
    if (c.get("source_archive_sha256")!=STAGE1_SHA or
        c.get("source_npz_sha256")!=npz_sha or
        c.get("format")!="FOLDBOUNDS_GNM_INTEGER_RESIDUAL_V1"):
        raise ValueError("ABSTAIN: wrong frozen experiment or certificate version")
    if c.get("integer_witness_sha256") != sha256(json.dumps(c["integer_inverse_proposal"],separators=(',',':')).encode("ascii")):
        raise ValueError("ABSTAIN: integer witness SHA256 mismatch")
    check=verify_and_score(array["reference_coords"],array["pair_i"],
                           array["pair_j"],c["denominator_Q"],c["integer_inverse_proposal"])
    expected={"score_numerators":[list(b) for b in check["score_numerators"]],
              "ordered_indices":check["ordered_indices"],
              "strict_adjacent_count":check["strict_adjacent_count"],
              "anchor":check["anchor"],"tree_parent":check["parent"],
              "tree_depth":check["depths"],"sum_tree_depths_B":check["B"],
              "n_contacts":check["contacts"],
              "n_pairs":len(array["pair_i"]),
              "n_residues":len(array["reference_coords"]),
              "max_residual_integer":check["max_residual_integer"]}
    for key,v in expected.items():
        if c.get(key)!=v:
            raise ValueError("ABSTAIN: certificate mismatch: "+key)
    frozen=array["prediction_GNM_resistance"]
    if np.max(abs(check["score_midpoints"]-frozen))>1e-11:
        raise ValueError("ABSTAIN: model mismatch with frozen baseline")
    print("VERIFIED exact GNM all-pair residual certificates; adjacent orders:",
          check["strict_adjacent_count"],"/",len(frozen)-1)
    return check


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("stage1_zip")
    ap.add_argument("certificate_json")
    ap.add_argument("--verify",action="store_true")
    args=ap.parse_args()
    if args.verify:verify_certificate(args.stage1_zip,args.certificate_json)
    else:make_certificate(args.stage1_zip,args.certificate_json)


if __name__=="__main__":main()