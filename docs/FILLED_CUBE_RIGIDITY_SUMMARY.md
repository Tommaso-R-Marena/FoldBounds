# Filled-cube contact-map rigidity and exact 3D ground-state entropy

Independent, nonconflicting research contribution. Full local Python/C++ implementation and test overlay is available as a downloadable archive in the associated ChatGPT conversation. This branch starts from main and leaves other agents' branches unchanged.

## Mathematical result

For a completely occupied d-dimensional m-by-m grid, the backbone edges together with all labeled nonconsecutive lattice-contact pairs recover the complete Cartesian-product grid graph. Its graph automorphisms are signed coordinate permutations; after fixing a rooted first +x edge, each full-grid contact map has exactly 2^(d-1)*(d-1)! rooted geometric realizations: 2 in 2D and **8 in 3D**.

An inverse geometry algorithm reconstructs the coordinates by choosing a corner, following uniquely recognized straight axes, and applying shortest-path distance landmarks:

    x_i(v) = (dist(c,v) + m-1 - dist(a_i,v))/2.

The algorithm checks a complete bijection with every grid cell and replays every backbone/contact edge, failing closed on incomplete or malformed graphs. In 3D, the orbit splits into 4 proper rotations and 4 reflected counterparts. Arbitrary residue-index-specific contact-only energy coefficients therefore cannot distinguish the mirror classes.

## Fully audited 3x3x3 all-H example

- 27-site maximum nonbackbone contact count: 28 (classical cubic edge-isoperimetric result).
- **826,768** anchored maximum-contact Hamiltonian paths (published crosscheck: OEIS A308586, directed count / 6).
- **103,346** distinct index-labeled maximum-contact maps.
- **8** geometric paths for every distinct map.
- Independent complete C++ path-count and contact-mask enumerators agree on all values; every complete contact map has multiplicity exactly 8.
- Exact free-space cubic rooted block counts for 0..10 steps: 1,5,25,121,589,2821,13565,64661,308981,1468313,6989025.
- The optimal block-product cap on all 27-site anchored walks is 137795893141213125, attained by suffix blocks (5,10,10).
- At q=10, the event of at least 16 contacts has certified probability >=98.3606413855% nominally and >=93.7499478473% under adversarial per-state factors in [1/2,2], uniformly for every real q>=10.
- Any one exact ground contact map has nominal population <=1/103346, and corrected population <=4/103349, for any q>=1.

The underlying graph automorphisms, Hamiltonian path counts, chirality limitation, and discrete edge-isoperimetry are established mathematics. The proposed software is **not Lean-verified** and the all-H cubic lattice is **not experimentally predictive biophysics**. Do not merge without source and scientific review.

## Local reproducibility

The complete source overlay contains two Python modules, three independent C++ reference programs, tests, a detailed proof argument, and a benchmark/demo:

    python -m pytest tests/test_cubic_ground_rigidity.py -q
    python experiments/cubic_ground_rigidity_demo.py --audit-cpp

Local outcome: 16/16 test functions passed, including exhaustive short-space oracles, full 3x3x3 C++ audits, large-box inverse reconstruction, and adversarial invalid-certificate tests.

Research provenance: https://oeis.org/A308586 ; Schram & Schiessel, 2013, Corrigendum on exact 4x4x4 Hamiltonian walks; cubic-lattice edge isoperimetry (https://pmc.ncbi.nlm.nih.gov/articles/PMC6733839/).
