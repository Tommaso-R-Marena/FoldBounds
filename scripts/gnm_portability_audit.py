#!/usr/bin/env python3
"""Exact scaled-integer GNM inverse test on previously frozen CspB/Crambin references.

Inputs are extracted ONLY from FROZEN_BEFORE_TARGET.npz members of the
previously recorded FoldBounds_Prospective_SEND_BACK.zip, never from NMR data.
The proposed inverses use SymPy; a separate integer checker verifies A*N=D*I.
"""
from __future__ import annotations
import argparse
import hashlib
from io import BytesIO
import json
import time
from zipfile import ZipFile
import numpy as np
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

ARCHIVE_SHA = '666a5d734fbd61be0f182c954474841bb552b71b3a1e9156d6c21b8c2f095f95'
FROZEN = {
    'CspB': ('4972bc4b7f580460607a1f6f27a5c44be0e4d675c39b9098a171529975d9f48c', 67, 1170, 543),
    'Crambin': ('a273deecac065246c5a891fe4fbeff7aaa380e26835a06f39bddeaf3b7253afd', 46, 547, 334),
}

def source(zip_path, tag):
    blob = open(zip_path, 'rb').read()
    if hashlib.sha256(blob).hexdigest() != ARCHIVE_SHA:
        raise ValueError('ABSTAIN: prospective outer ZIP differs from original freeze')
    with ZipFile(BytesIO(blob)) as z:
        manifest = json.loads(z.read('SHA256SUMS.json'))
        if set(manifest) != set(z.namelist()) - {'SHA256SUMS.json'}:
            raise ValueError('ABSTAIN: prospective ZIP manifest invalid')
        # Do NOT read target evaluation members, even to recompute hashes.
        # Exact outer archive SHA authenticates all bytes without opening them.
        filename = f'{tag}_FROZEN_BEFORE_TARGET.npz'
        b = z.read(filename)
        if manifest.get(filename) != FROZEN[tag][0] or hashlib.sha256(b).hexdigest() != FROZEN[tag][0]:
            raise ValueError('ABSTAIN: reference-only npz differs from frozen bytes')
    with np.load(BytesIO(b), allow_pickle=False) as ar:
        coords = ar['reference_coords'].copy()
        pairs = list(zip(map(int, ar['pair_i']),map(int, ar['pair_j'])))
        frozen_gnm = ar['prediction_GNM_resistance'].copy()
    _, n, m, _ = FROZEN[tag]
    if coords.shape != (n,3) or len(pairs) != m or any(i<0 or j<=i or j>=n for i,j in pairs):
        raise ValueError('ABSTAIN: wrong frozen reference geometry/pairs')
    if len(set(pairs)) != m or not np.isfinite(coords).all() or not np.isfinite(frozen_gnm).all():
        raise ValueError('ABSTAIN: duplicate pair or invalid coordinates/predictions')
    return coords,pairs,frozen_gnm


def contacts(coords):
    n = len(coords)
    # Exact binary64 rational distances, no numerical cutoff ambiguity.
    xyz = [[float(v).as_integer_ratio() for v in row] for row in coords]
    def within(i,j):
        from fractions import Fraction
        return sum((Fraction(*xyz[i][k]) - Fraction(*xyz[j][k]))**2 for k in range(3)) <= 100
    return [(i,j) for i in range(n) for j in range(i+1,n) if within(i,j)]


def integer_graph(edges,n,anchor):
    ix={k:(k if k<anchor else k-1) for k in range(n) if k!=anchor}
    A=[[0]*(n-1) for _ in range(n-1)]
    for a,b in edges:
        if a!=anchor:A[ix[a]][ix[a]]+=1
        if b!=anchor:A[ix[b]][ix[b]]+=1
        if a!=anchor and b!=anchor:
            A[ix[a]][ix[b]]-=1; A[ix[b]][ix[a]]-=1
    return A,ix


def audit(zip_path,tag):
    t=time.monotonic()
    coords,pairs,old=source(zip_path,tag)
    edges=contacts(coords)
    if len(edges) != FROZEN[tag][3]:
        raise ValueError(f'ABSTAIN: changed contact count {len(edges)} for {tag}')
    n=len(coords)
    anchor=n-1
    A,ix=integer_graph(edges,n,anchor)
    NN,DD=DomainMatrix.from_list(A,ZZ).inv_den()
    N=[[int(v) for v in row] for row in NN.to_list()]
    D=int(DD)
    if D<0: N=[[-v for v in row] for row in N];D=-D
    # Checked independently of SymPy via integer sparse multiply; generic verifier.
    # validate_matrix in RNase module hard-codes 124: check a native variable-size identity here.
    if D<=0 or any(N[i][j]!=N[j][i] for i in range(n-1) for j in range(i+1,n-1)):
        raise ValueError('ABSTAIN: denominator/sign/asymmetry')
    sparse=[[(j,c) for j,c in enumerate(row) if c] for row in A]
    for i in range(n-1):
        for j in range(n-1):
            if sum(v*N[k][j] for k,v in sparse[i])!=(D if i==j else 0):
                raise ValueError('ABSTAIN: integer identity A*N=D*I failed')
    nums=[]
    for i,j in pairs:
        a,b=ix.get(i),ix.get(j)
        v=(N[a][a] if a is not None else 0)+(N[b][b] if b is not None else 0)-2*(N[a][b] if a is not None and b is not None else 0)
        if v<=0:raise ValueError('ABSTAIN: nonpositive exact resistance')
        nums.append(v)
    from fractions import Fraction
    exact_order=sorted(range(len(nums)),key=lambda k:(nums[k],k))
    floating=[float(Fraction(v,D)) for v in nums]
    float_order=list(np.argsort(old,kind='stable'))
    rank_same=exact_order==float_order
    maxdiff=float(np.max(np.abs(np.asarray(floating)-old)))
    # The exact CspB scores have nine ties; near-equal frozen floats can
    # permute those ties. Compare only mathematically strict orderings and
    # allow 1e-12 for the frozen binary64 approximation.
    tolerance=1e-12
    checked=0; violations=0; tied_pairs=0
    for a_pos in range(len(exact_order)):
        i=exact_order[a_pos]
        for b_pos in range(a_pos+1,len(exact_order)):
            j=exact_order[b_pos]
            if nums[i]==nums[j]:
                tied_pairs+=1
                continue
            checked+=1
            if old[i]>old[j]+tolerance:violations+=1
    if violations or maxdiff>1e-10:
        raise ValueError('ABSTAIN: frozen strict orders disagree with exact GNM')
    return dict(protein=tag,residues=n,contact_edges=len(edges),pair_scores=len(nums),
                exact_ties=len(nums)-len(set(nums)),
                exact_total_sort_equal_to_floating=rank_same,
                frozen_strict_order_violations=violations,
                frozen_strict_order_comparisons=checked,
                exact_tied_pair_comparisons=tied_pairs,
                max_float_difference=maxdiff,
                denominator_digits=len(str(D)),time_seconds=round(time.monotonic()-t,3),
                frozen_reference_sha256=FROZEN[tag][0],
                matrix_identity_checked=True, target_data_accessed=False)

if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('prospective_zip')
    a=ap.parse_args()
    result=[audit(a.prospective_zip,x) for x in FROZEN]
    print(json.dumps(result,indent=2))
    if any(r['frozen_strict_order_violations'] or r['max_float_difference']>1e-10 for r in result):
        raise SystemExit('ABSTAIN: frozen baseline mismatch')