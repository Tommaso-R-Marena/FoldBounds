# Exact inverse contact geometry and 3D ground-state entropy

**Status:** independent, additive FoldBounds reference. Python/C++ software assurance only, not Lean-kernel proof or physical protein prediction. The cube grid automorphism theorem and lattice isoperimetric theorem are classical. Do not merge without review.

## 1. Mathematical statement: complete contact graph determines embedding modulo rigid symmetries

For a Hamiltonian walk through every vertex of a square or cubic grid of side m, the index-labeled backbone together with *every* nonconsecutive nearest-neighbor contact is exactly the graph P_m^d. There is no unknown adjacency: every lattice-neighbor pair is either a backbone edge or a contact edge.

The inverse algorithm starts with a vertex c of degree d (a corner), chooses an ordered axis basis from its d neighbors, then follows straight grid axes. A straight successor two steps beyond c shares exactly one neighbor with c, whereas a right-angle turn shares two. This yields opposite axis endpoint a_i. Shortest-path distances reconstruct each lattice coordinate:

    x_i(v) = (dist(c,v) + m-1 - dist(a_i,v))/2.

The checker refuses nonintegers, out-of-box coordinates, collisions, and any disagreement between the supplied contact+backbone graph and all grid adjacencies computed from the reconstructed coordinates. No graph-isomorphism oracle is trusted. Every accepted candidate is a full integer-vertex cube embedding.

An automorphism of P_m^d permutes its d Cartesian axes and independently flips them. The group has 2^d*d! elements. A canonical rooted self-avoiding chain begins at origin with its first directed bond +x, so exactly

    2^(d-1)*(d-1)!

distinct anchored embeddings realize every accepted complete labeled contact map. This is 2 in 2D and 8 in 3D. In 3D, four are proper rotations around the oriented x axis and four are reflected equivalents. The latter divide into two proper-rotation orbits, which is an identifiability limitation for contact-only Hamiltonians, not a claim about biochemical enantiomerism.

Even if the pairwise contact energy E is a separate arbitrary residue-index-specific rational coefficient for every nonconsecutive pair, E takes the exact same value on all eight map-equivalent states. Torsions, hydrogen bond orientation, handed backbone stereochemistry, and solvation are *not* assumed to obey this invariance.

## 2. Exact all-hydrophobic 27-site ground-state census

For a filled m-by-m-by-m cube, the number of occupied adjacent site-pairs is 3*m*m*(m-1). Subtract n-1 backbone edges to get nonbonded contacts:

    K_cube = 2*m^3 - 3*m^2 + 1.

The classical edge-isoperimetric equality case at n=m^3 states this is the global maximum for any connected n-site lattice set, attained only by a translated filled cube. All Hamiltonian paths of such a cube are exactly maximal-contact all-H SAWs.

Two standalone C++17 enumerators reconstruct the n=27 anchored census without sharing traversal logic:

- **826768** first-East Hamiltonian paths on a completely occupied 3x3x3 cube. This agrees with OEIS A308586: 2480304 undirected Hamiltonian paths, times two for direction, divided by six for first-bond rotational normalization.
- **103346** distinct index-labeled nonbackbone contact masks, each realized exactly **8** times.
- For 2x2x2: **24** anchored Hamiltonian paths and **3** contact masks, each realized 8 times.
- Every mask is a full 192-bit value with equality checks; hashing is used for storage but cannot merge unequal maps by hash collision.

The connectivity tests prune only unvisited regions that are disconnected from the current endpoint. The separate map enumerator recomputes contacts from coordinate adjacency rather than relying on a precomputed score from a Python checker.

## 3. Thermodynamic occupancy and rigorous no-go

Write C=826768 exact ground states at maximum score K=28 and let the positive statewise correction a_s range independently within [l,u], where l>0. If A_t is the event that a conformation has at least t contacts, t<=K, then every ground state is in A_t. Every state outside A_t has score at most t-1.

A bound W on the **total** number of anchored 27-site walks gives

    P_q(A_t) >= l*C*q^(K-t+1) / (l*C*q^(K-t+1)+u*W).

This is valid for q>=1 and monotone in q, so the rational checkpoint at q0 proves the same lower bound for every **real** q>=q0, even if per-state correction factors vary with q within [l,u].

To certify W, a walk's n-2 steps beyond the anchored first edge are cut into blocks. For a fixed block length k, the number of legal next k-step continuations is never greater than the independently enumerated free-space rooted count c[k], because discarding older occupied sites only permits *more* extensions. Therefore W<=product c[k_i] for every partition sum(k_i)=n-2. A finite exact-integer DP picks the smallest product. The separate C++ free-space audit checks

    c[0..10] = (1,5,25,121,589,2821,13565,64661,308981,1468313,6989025).

For n=27, the minimizing block partition is (5,10,10), yielding

    W=137795893141213125.

For q=10 and t=16, this yields **P(A_16)>=98.3606413855% nominally** and **>=93.7499478473%** with a_s in [1/2,2], each uniform for every q>=10.

For any ONE exact ground microstate and any ONE full-ground contact map, the presence of other equal-energy ground states gives respectively

    P(one ground microstate) <= u/[u+l*(C-1)],
    P(one full-ground map) <= 8*u/[8*u+l*(C-8)].

These no-go ceilings are valid for all q>=1, and for the nominal case simplify to 1/826768 and 1/103346. Under [1/2,2], the latter is 32/(C+24). This bound is sharp if the ensemble includes only the C ground states with adversarial factors at the corners; the actual full ensemble also has excited states.

This is a matched-model distinction between compactness as a *bulk event* and identifying a *specific geometric state* from its contacts. Rotations and translations should be quotient-equivalent when asking about physically distinct structures; the remaining pair of mirror classes is the contact-only ambiguity. Do not interpret all eight rooted lattice realizations as eight different naturally occurring proteins.

## 4. Reproduction and assurance

The PR provides two standalone Python modules, independent C++17 full-cube enumerators, a separate C++ free-space walk counter, a demo, and pytest checks.

    python -m pip install -e '.[test]'
    python -m pytest tests/test_cubic_ground_rigidity.py -q
    python experiments/cubic_ground_rigidity_demo.py --audit-cpp

Local result: 16 pytest parameterizations pass. Tests include complete independent 2x2 square, 3x3 square, and 2x2x2 cube direct enumeration; one full 3x3x3 C++ recount for every map; exact n=8 unrestricted Gibbs probabilities; forged/tampered inputs, out-of-range indices, impossible maps, and rational pairwise energy invariance. g++ is needed to replay the external audits.

No Python or C++ theorem is Lean-refined. The classical grid-product automorphism and cubic edge-isoperimetric equality results are mathematical dependencies. The count tables are fixed after the independent C++ enumeration and rechecked on full CI when the compiler is available. These are deliberately simplified HP model demonstrations, without experimental temperature calibration, folding kinetics, native-structure accuracy, or theorem novelty claims.

**Prior art:** OEIS A308586 (https://oeis.org/A308586); the classical edge-isoperimetric problem on cubic lattices and grid graph automorphisms; existing lattice-polymer work on contact-map degeneracy and mirror symmetry.