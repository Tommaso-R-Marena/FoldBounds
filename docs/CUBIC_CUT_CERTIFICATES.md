# Independently checked 3D HP cut-tree certificates and all-q native-mass bounds

**Scope:** additive, opt-in research code stacked on draft PR #64 (itself #61). It does not change any existing Lean source, production model solver, other research branch, blind structural experiment, or CI configuration. The new producer, independent verifier, unit tests and demonstration live in separate files.

**Research integrity:** This is *not* a new protein-folding law, original discrete Loomis–Whitney theorem, or experimental validation. Branch-and-bound and HP lattice ground-state searches have prior art. The new result is a small, reproducible, independently rechecked **3D exact-tail reference oracle** and an exact continuum bound for its thermal population, without enumerating the full 3D sequence-conformation space. The C++ executable code is **not Lean-kernel-refined**. Passing any preexisting Lean audit cannot change that.

## 1. Mathematical search certificate

Fix a string `s` of n H/P beads and the 3D simple-cubic polymer lattice. The first directed chain edge is always `(0,0,0)->(1,0,0)`. Every later vertex occupies a free lattice neighbor of the chain tip. A nonbonded H–H contact is a lattice-neighbor pair of H beads with sequence-index separation ≥3. The exact contact-score of a completed walk is K.

The question is to count **all** rooted walks with `K >= threshold`, without visiting every complete walk. At each partial walk of d vertices, the producer proposes a certificate bit:

- `1` = prune the entire subtree; an independently recomputed **integer upper bound** must be strictly smaller than `threshold`.
- `0` = expand every legal neighboring lattice step in the fixed six-direction order; at completed leaves, the verifier directly rescores the entire coordinate path.

The fixed-direction traversal and self-avoidance rules determine the tree, so the single-bit-per-visited-node **preorder cut stream** is sufficient to specify *every* search branch. The independent verifier rejects missing/extra bits, a false summary count, malformed headers/padding, excessive resource claims, or any pruning decision that fails its independently computed bound. It recomputes partial geometry and contact scores from scratch rather than importing the producer's incremental state. No committed certificate JSON claims are trusted on their own.

### A. Per-residue contact capacity and parity

Every possible *future* H–H contact involves two hydrophobic residue indices of opposite index parity (because Z³ is bipartite), with separation ≥3. Already-placed/placed contacts are either known or impossible to add later. A remaining pair joining a placed H at index i to a future H at j is possible only if its place can be reached in time: 

`ManhattanDistance(tip, position[i]) <= (j - d + 1) + 1`.

For a placed H, count its *currently empty* lattice-neighbor positions. If that residue is the current tip and completion remains, reserve one such position for the mandatory **next backbone bond**. For an unplaced H, permit at most four nonbonded neighbors at an internal chain position and five at the terminal chain position. These deliberately optimistic capacities are valid for every geometrically legal completion.

Build the list E of candidate future H–H residue pairs meeting all *necessary* conditions above. For each H index i let `eligibleDegree[i]` be its degree in that candidate-pair graph and `cap[i]` its remaining nonbonded contact-site capacity. Every newly formed H–H contact uses one even-indexed and one odd-indexed H. Consequently:

`remainingHH <= min(|E|, sum_even min(cap[i],eligibleDegree[i]), sum_odd min(cap[i],eligibleDegree[i]))`.

This is an upper relaxation, not a constructive matching proof. A realizable completion can never exceed it.

### B. 3D discrete Loomis–Whitney projection upper bound

For any finite complete occupied vertex set V of size n, let `p_xy`, `p_xz`, `p_yz` denote the numbers of distinct two-coordinate projected positions. Each x-parallel line with r occupied vertices contains at most r−1 lattice edges, so

`E_lattice(V) <= 3n − (p_xy+p_xz+p_yz)`.

The **discrete Loomis–Whitney inequality** states `p_xy * p_xz * p_yz >= n²`. Every completed walk also has projection cardinalities at least as large as those already realized by its prefix. Define `Smin` as the minimum **integer** `a+b+c` over `prefixPxy<=a<=n`, `prefixPxz<=b<=n`, `prefixPyz<=c<=n`, with `a*b*c>=n²`. This is a finite integer problem, so no floating-point geometry is used.

Every completed anchored walk contains exactly n−1 **backbone edges**, hence its total number of nonconsecutive lattice-neighbor contacts (including H–P and P–P) is at most

`Ctotal <= 3n − Smin − (n−1)`.

If `Call_prefix` of these nonconsecutive contacts are already permanently realized in the partial walk and `KHH_prefix` are H–H, then

`remainingHH <= Ctotal − Call_prefix`.

The producer uses a dynamic-programming table of integer `Smin`; the independent verifier computes the minimum directly by independent nested integer loops. The **final safe upper bound** is `KHH_prefix` plus the minimum of the above three upper constraints and the projection-geometric remaining-contact bound. If this is less than the query threshold, the entire partial-walk subtree is correctly pruned.

For background on the discrete Loomis–Whitney inequality, see [Bounds on the Lattice Point Enumerator via Slices and Projections](https://doi.org/10.1007/s00454-021-00310-7), equation (1.2). The method here applies that *existing* result; it does not prove the original inequality anew.

## 2. Exact 16-mer contact-energy gap without full conformational enumeration

Sequence: **`HHHPPHHHHHHHPPHH`** (16 residues, 12 H). The independent verifier accepted three distinct byte-for-byte cut-tree certificates:

| Question | Cut-tree nodes visited | Verified pruned prefixes | Exact number of completions |
|---|---:|---:|---:|
| K≥10 | 25,998,290 | 20,300,541 | **2,032** |
| K≥11 | 10,851,705 | 8,444,857 | **2,032** |
| K≥12 | 4,071,184 | 3,163,781 | **0** |

A separately checked explicit 16-site, self-avoiding Hamiltonian coordinate witness has **K=9**; another has **K=11**. Therefore:

- The exact maximum H–H contact count is **M=11**.
- There are exactly **g=2,032** conformations at that score, with the same rooted first edge normalization.
- There are **zero** score-10 conformations; score 9 does occur. The first excitation gap is **exactly two model contact energies**.
- These 2,032 conformations are *not* claimed to be a single fold or one spatial symmetry orbit.

**No exhaustive n=16 conformational oracle was used.** As an external scale reference, [OEIS A001412](https://oeis.org/A001412) lists **21,175,146,054** unrestricted-start-direction 15-step cubic self-avoiding walks. Fixing one of the six equivalent first-edge directions leaves **3,529,191,009** full anchored walks. The 10,851,705-node certificate is far smaller than explicitly constructing those complete conformations, although its checker still explores millions of partial prefixes and is not polynomial-time in chain length.

The compact *binary* certificates have one preorder bit per visited node and a small header: approximately **3.1 MB, 1.3 MB, and 0.5 MB** for thresholds 10, 11 and 12. The checker does **not** trust cryptographic hashes in place of its own traversal; SHA-256 manifests are provided only for byte-level reproduction and data integrity.

## 3. All-q **upper Boltzmann partition envelope** from checked pruned subtrees

The independent checker returns more than a yes/no decision for each cut: it calculates the contact-score upper bound `U(prefix)<threshold` and the **maximum number of full nonbacktracking completions** beneath that prefix, `5^(n−depth)`. Because the cut prefixes are disjoint, this gives a **separately checked energy-resolved upper envelope**, counting every possible true completion at most via a positive upper majorant.

For a checked threshold-10 proof, write `C[k]` for the sum of `5^(n−depth)` over cut prefixes with `U(prefix)=k`, plus any directly visited non-ground terminal leaves of score k. The independently verified coefficients are:

| Contact-score cap k | Checked maximum completion count C[k] |
|---:|---:|
| 1 | 4,760 |
| 2 | 44,500 |
| 3 | 953,160 |
| 4 | 6,145,140 |
| 5 | 34,661,428 |
| 6 | 235,667,776 |
| 7 | 528,467,884 |
| 8 | 1,598,203,797 |
| 9 | 1,788,702,372 |

All 2,032 non-cut leaves have K=11, so **for every real q>=1**,

`Z_non_ground(q) <= sum_{k=1}^9 C[k] * q^k`.

This is an **upper bound** on Boltzmann weight, not a reconstructed exact density of states. It can greatly overcount conformations, but cannot miss any if the pruning inequalities and independent checker are correct. The separately checked threshold-11 cut tree yields another valid upper polynomial, and the bound `Z_non_ground(q) <= (5^(n−2)−g)*q^(M−2)` holds without examining any prefixes. Taking the *minimum* of their three independently valid upper envelopes gives the strongest checked bound in this module.

The all-contact fixed-first-edge population satisfies `N <= 5^(n−2)` because every later bead has at most five nonbacktracking directions. Consequently the **ground-population probability lower bound** is

`P_ground(q) = g*q^M / (g*q^M + Z_non_ground(q))`

`>= g*q^M / (g*q^M + min(U_10(q),U_11(q),(5^14-g)*q^9))`.

This is a **continuum** inequality, not a finite-temperature-grid experiment. Only exact integer counts and rational arithmetic enter the verifier and the Python consumer.

| Dimensionless contact weight q | Certified lower bound P_ground(q) |
|---:|---:|| 100 | >1.11% |
| 1,000 | >53.16% |
| 2,000 | >81.95% |
| **2,816** | **>90.005%** |
| 5,000 | >96.59% |
| 10,000 | >99.12% |

The stronger score-resolved envelope guarantees **90% ground occupancy for every real q>=2816**. The coarse total-walk-only bound would require integer q>=5200. Neither value is Kelvin, a physical native-state melting point, nor a measured folding result; these enormous toy contact weights are a warning about transfer to biology.

## 4. Run, tamper, and reproducibility

The producer and checker are **separate C++17 programs**. The consumer is Python 3.11+. They do not trust upstream data and do not modify the main FoldBounds verifier:

```bash
python -m pytest tests/test_cubic_cut_certificates.py -q
python experiments/cubic_cut_demo.py --output-dir ./cubic_cut_proof
# Or invoke directly:
g++ -std=c++17 -O3 experiments/cubic_cut_producer.cpp -o cubic_cut_producer
g++ -std=c++17 -O3 experiments/cubic_cut_checker.cpp -o cubic_cut_checker
./cubic_cut_producer HHHPPHHHHHHHPPHH 11 ground.cut geometry
./cubic_cut_checker ground.cut
```

The tests compare many independently checked C++ threshold counts against the existing, fully audited 3D exact atlas for short chains, verify both `capacity` and `geometry` producer modes agree, and reject false root cuts, forged counts, truncated/extended payloads, mismatched thresholds, and invalid geometric witnesses. One optional heavy test executes all three n16 cut trees and checks the exact two-contact gap and native-occupancy bound.

The producer is untrusted. Its certificate is accepted only after the separate checker reconstructs the legal successor search tree, checks **every** cut, scores every surviving terminal configuration, checks cardinalities and bitstream termination, and computes score-resolved partition majorants. The checker itself is ordinary C++ and thus remains in the trusted computing base; there is **no Lean-kernel executable refinement** of this step. An independent Lean formalization of cut-tree coverage and the projection inequalities is a priority before presenting stronger formal-verification claims.

## 5. Scientific significance and limitations

This contribution demonstrates a checkable **3D, not 2D**, global high-contact tail and continuum-native-mass certificate that does not enumerate billions of completed conformations. It is a meaningful exact-reference capability and a stronger target for FoldBounds' independent certifier research. It does not generalize to realistic protein potentials or experimentally supported folding behavior without substantial new modeling, testing, and independent biological evidence. The certificate checker still has exponential worst-case scaling and is **not** a general 16-residue uniqueness or protein-design solution.
