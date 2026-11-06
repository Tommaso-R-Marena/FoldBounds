# Exact sparse-distance reconstruction on experimental protein coordinate models

**Research status:** exact integer mathematics and Python-checked certificates,
not a new theorem of Euclidean distance geometry, Lean-kernel proof, de novo
fold prediction, or calibration to actual instrument uncertainty.

## Two complementary certificates

**Four known, affinely independent reference sites** a0,a1,a2,a3 fix the
reference frame. Write B as the 3x3 matrix of rows a_i-a0, i=1..3. If
s_i=||x-a_i||^2, subtraction of s0 from each other squared distance yields

    2 B (x-a0) = (||a_i-a0||^2+s0-s_i)_(i=1..3).

Invertible B makes every remaining coordinate unique. In a frame-free
description, the four-anchor distance graph has six distances within the
tetrahedron and four from each of the n-4 remaining sites: **4n-10 distances**.
The executable instead accepts the four *Cartesian reference coordinates* and
4(n-4) additional distances. It is not a claim to compress 3n coordinates
into fewer bits. The verifier computes determinant, adjugate and exact
divisibility in arbitrary-precision integers, and checks all anchor distances,
coordinate digest, and typed inputs.

**Three known, noncollinear reference sites** require only 3n-6 distances
(three among anchors plus 3(n-3) from nonanchors) and n-3 orientation signs.
Let u=a1-a0, v=a2-a0, N=u cross v, G=||N||^2, and
t1=(||u||^2+s0-s1)/2, t2=(||v||^2+s0-s2)/2. The squared signed
triple product obeys the exact Gram determinant identity

    V^2=G*s0-(||v||^2*t1^2+||u||^2*t2^2-2(u.v)*t1*t2).

The sign of V chooses the point above/below the reference plane.
Solve the linear system with rows u,v,N and right-hand side t1,t2,V.
The verifier rejects negative/non-square V^2, nonintegral solutions, wrong
distance squares, inconsistent signs, and source-digest mismatches.
On these four deposited structures, all orientation signs are ±1.
This tradeoff is classical trilateration with chirality information.

## Universal adversarial error certificate (four anchors)

Assume the *anchor coordinates are exact*. For every other residue i, each
of its four observed anchor distances may be changed independently by any
real amount of absolute value at most eps. Let D_i be a conservative integer
upper bound on the four original distances and

    eta_i = 2 D_i eps + eps^2,
    R = ||B^{-1}||_infinity, h_i=R eta_i.

Each squared distance changes by <=eta_i. Subtracting the equations gives
a right-hand-side error bounded by eta_i, so every recovered site satisfies
||x'_i-x_i||_infinity <=h_i. A pair distance can therefore change by at most
sqrt(3)(h_i+h_j), safely relaxed to 2(h_i+h_j). At cutoff r, mark:

- guaranteed contact if original distance <= r-2(h_i+h_j);
- guaranteed noncontact if original distance > r+2(h_i+h_j);
- otherwise **undecided**, never silently assert certainty.

The implementation uses exact **integer squared-distance inequalities**
and cross multiplication of rational radii, with no floating-point
roundoff in the final decision. Anchor selection is a deterministic greedy
farthest/triangle-area/tetrahedral-volume heuristic, **not** a global optimum.

The alternate three-anchor signed certificate does not inherit this linear
noise guarantee near its reference plane.

## Frozen, real deposited coordinate models

Public provenance: biopython/biopython \`Tests/PDB\`, originally deposited
in the Protein Data Bank. Parse first MODEL, chain A, unique ATOM C-alpha
records in file order, coordinates at the deposited 0.001 Angstrom precision.
All source Git blob SHA values are in \`data/experimental_ca_packed.py\` and
\`data/experimental_ca_7DDO.py\`. The source can be independently re-extracted.

| PDB | Experimental technique | C-alpha sites | Eligible pairs (sequence separation>=4) | Nominal 8A contacts | Guaranteed contacts | Undecided pairs at eps=0.01A |
|---|---|---:|---:|---:|---:|---:|
| 1A8O | X-ray | 66 | 1,953 | 88 | 78 | 17 |
| 1LCD | solution NMR | 51 | 1,128 | 94 | 87 | 18 |
| 2BEG | solution NMR | 26 | 253 | 1 | 0 | 2 |
| 7DDO | electron microscopy | 597 | 176,121 | 1,235 | 1,123 | 229 |

At eps=0.01A, naive first-four-residue reference tetrahedra left
411,117,13,5693 eligible pairs undecided; the farthest-volume rule left
17,18,2,229, respectively. The improvement is **in certified labels for
these deposited coordinates**, not prediction accuracy or natural
uncertainty calibration. Very large fractions of certifiable pairs are
noncontacts, so the contact-specific counts are provided separately.

Exact 3-anchor reconstructions succeeded for all 740 CA sites; the largest
requires 1,785 anchor-related distances (including the three
anchor-anchor distances), plus 594 orientation sign bits. Exact reference quadruple decoding also succeeded.
A separately implemented BigInt geometry checker directly consumed all
four pinned source PDBs and reproduced every four-anchor count, exact
coordinate reconstruction and three-anchor Gram perfect-square test.

## Limitations and reproducibility

Source coordinates were extracted from the deposited structure models;
the method does **not** derive structures from experimental NMR restraints,
X-ray structure factors, EM density, or sequence alone. The 0.01A and
0.05A perturbations are hypothetical, NOT instrument error estimates.
NMR ensemble model dispersion is not an equilibrium trajectory. C-alpha
contacts omit chemistry, solvation, sterics, torsion and chirality.
The selected public PDB test fixtures form a convenience sample, NOT a
representative external benchmark. These standard geometric identities
have established prior art; no mathematical priority is claimed.

    python -m pip install -e '.[test]'
    python -m pytest tests/test_sparse_distance_realdata.py -q
    python experiments/sparse_geometry_public_demo.py

No Lean proof of these new Python functions exists. Keep this PR draft until
the independent source audit, full CI, and scientific review pass.


## Local minimality of the signed metric data

For n generic labeled points in 3D, configurations modulo the six
continuous degrees of Euclidean rigid motion have dimension 3n−6.
Any system of m smooth scalar distance measurements has Jacobian
rank at most m. If m<3n−6, it cannot give a full-rank *local*
coordinate chart on this configuration space. This is the ordinary
rank/inverse-function theorem, not a new lower-bound theorem.

The signed three-anchor construction uses exactly 3n−6 pair distances
(three among the anchors plus three per nonanchor). For every site
x with nonzero signed volume V=(u×v)·(x−a0), the Jacobian of its
three squared distances to fixed anchors is the 3×3 matrix with
rows 2(x−a0),2(x−a1),2(x−a2), with exact determinant **8V**.
Hence each nonanchor block is nonsingular. The three mutual anchor
distances fix the reference triangle modulo rigid motion, and the
discrete orientation signs choose one of the mirror branches.
This is locally full-rank using the minimum number of independent
*continuous* metric scalars, conditional on the sign side-information.

A new regression computes the determinant using integer operations
on every nonanchor site in all four deposited structures: none has
V=0 for the chosen anchors. This **does not** imply that arbitrary
distance graphs with 3n−6 edges are globally rigid, nor that one can
discard the orientation signs, nor that bit-level compression is optimal.
