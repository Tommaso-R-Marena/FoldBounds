# FoldBounds: exact joint-covariance GNM ranking certificates on real ubiquitin 1D3Z

**Status:** independent computational methods result, not protein-folding prediction or mathematical-priority claim. No results from a sealed/held-out target were accessed.

## Question

Can interval guarantees for per-residue Gaussian-network-model (GNM) flexibility be tightened when residue coordinates and hence the underlying contact graph are uncertain?

The previous conservative strategy certifies the ordering `C_ii > C_jj` by comparing a lower bound on `C_ii` to an upper bound on `C_jj` independently. This loses correlation between the two diagonal terms.

## Exact rank-two interval theorem

Let \(L_-\preceq L\preceq L_+\) be graph Laplacians of an *arbitrary* connected optional-contact network bounded by its certainly-on and possibly-on contacts (unit springs including mandatory covalent chain edges). Work on the zero-mean subspace:

\[
C_+ = L_+^+ \preceq C = L^+ \preceq C_- = L_-^+,\qquad D=C_- - C_+\succeq 0 .
\]

For two distinct residues \(i,j\), put
\(a=D_{ii}, b=D_{jj}, c=D_{ij}\) and \(f=C_{+,ii}-C_{+,jj}\).
Over *every* symmetric matrix \(C\) in the PSD interval \([C_+,C_-]\), the attainable contrast range is exactly

\[
f+\frac{a-b-\sqrt{(a+b)^2-4c^2}}2
\;\le C_{ii}-C_{jj}\le\;
f+\frac{a-b+\sqrt{(a+b)^2-4c^2}}2 .
\]

This is the optimal solution of a rank-two semidefinite interval problem (it is **not** necessarily attainable by a physical contact graph, so it remains conservative for the graph family).

**Proof sketch.** Let \(X=C-C_+\), so \(0\preceq X\preceq D\). For \(H=e_ie_i^\top-e_je_j^\top\), optimizing \(\mathrm{tr}(HX)\) over that interval reduces to the positive/negative spectral parts of \(D^{1/2}HD^{1/2}\). This rank-two symmetric matrix has trace \(a-b\), product of its two nonzero eigenvalues \(c^2-ab\le0\), hence eigenvalues \((a-b\pm\sqrt{(a+b)^2-4c^2})/2\). Their sums of positive and negative parts give the displayed sharp interval. For an actual intermediate contact graph, the PSD sandwich follows from inversion reversing Loewner order on the connected zero-mean subspace.

The conventional separate diagonal interval is
\(f-b\le C_{ii}-C_{jj}\le f+a\), and is never tighter.
A strict signed certificate uses **no floating-point square root**:
set \(s=2f+a-b\), \(\Delta=(a+b)^2-4c^2\).
The difference is definitely positive when \(s>0\) and \(s^2>\Delta\);
definitely negative when \(s<0\) and \(s^2>\Delta\).
All strict tests are exact-integer inequalities after multiplication by positive common inverse denominators.

This rank-two PSD optimization and Laplacian inverse monotonicity are standard linear algebra; we claim a **new independently audited implementation and measured model-certificate improvement inside FoldBounds**, not a newly invented theorem.

## Producer untrusted, checker exact

For a connected graph on \(n\) sites, form the strictly SPD integer matrix

\[
B=nL+ \mathbf{1}\mathbf{1}^{\top}.
\]

An untrusted SymPy `DomainMatrix.inv_den()` producer emits integer \(N,d>0\). The **separate stdlib-only checker** independently builds \(B\) from the declared graph and multiplies all \(n^2\) entries to verify \(BN=dI\), also checking symmetry and dimensions. Then \(C=nN/d-\mathbf{1}\mathbf{1}^{\top}/n\) is the **exact** GNM pseudoinverse, without numerical condition-number or residual assumptions. Both endpoint graph witnesses are verified before any count is reported.

For optional contacts arising from independent ±epsilon Cartesian boxes about deposited Cα coordinates, the pairwise Euclidean distance shift is at most \(2\sqrt3\,\epsilon < (7/2)\epsilon\). The checker classifies guaranteed/possible contacts using **exact rational distance-squared tests** and the slightly conservative \(7\epsilon/2\) slack, not floating-point cutoff comparisons. Contact definition: \(|i-j|\ge 3\), Cα distance strictly <8 Å, plus mandatory backbone edges.

## Actual public structural evidence, no target tuning

- Protein: human ubiquitin, deposited NMR ensemble **PDB 1D3Z** (ten submitted conformers); not a validated equilibrium dynamical trajectory.
- Test: **model 1**, Cα residues **2..70** selected in public ENCORE example `filtered_1D3Z.pdb`.
- Archived source: `encore-similarity/encore`, commit `2c99cb3fdf345c7679801378ec63634e1d007a5f`, file `examples/experimental_ubq_ensembles/filtered_1D3Z.pdb`, Git blob SHA-1 `a438d36816e15647177e32d2eafa6aa9a0782747`. The checked-in 69-position fixture retains exact deposited **0.001 Å** precision. Opt-in provenance script hashes raw source and verifies every coordinate and all ten model identities.
- Hypothetical independent per-Cartesian-coordinate structural uncertainty `epsilon`. This is **not** a calibrated experimental measurement.
- Baseline and improved **conditional GNM rank certificates**, exactly verified by integer graph-inverse identities:

| ε (Å) | Guaranteed contact edges (nonlocal) | Possible nonlocal edges | Old separate intervals /2346 | New joint intervals /2346 | Strictly added |
|---:|---:|---:|---:|---:|---:|
| 0.025 | 154 | 161 | 2189 | **2207** | 18 |
| 0.05 | 152 | 167 | 1884 | **1919** | **35** |
| 0.10 | 145 | 188 | 1065 | **1131** | **66** |
| 0.20 | 137 | 207 | 629 | **699** | **70** |
| 0.30 | 124 | 245 | 198 | **246** | **48** |

The largest *relative* gain of these five predeclared uncertainty settings is **+24.24%** (198→246, ε=0.30 Å), but only 246/2346 residue pairs remain provably ordered there. At the 0.05 Å condition, 35 additional pairs are gained (1884→1919, +1.86% relative). For each row the result is valid over **all combinations of optional edges**, not just the two endpoint networks. The proof is conditional on the exact source coordinates, hypothetical uncertainty model, unit-spring GNM, correct graph specification, and exact verifier code.

## Independent tests

- All 64 seeded small connected graph fixtures (n=5..8) were compared to a **separate Fraction Gaussian-elimination oracle** over *every* combination of their optional edges. Every signed pair certificate was confirmed against every exact graph state.
- Mathematically invalid graphs, non-PD input, malformed/negative/tampered integer inverse witnesses, and bypassing exact verification are rejected.
- Five 69-residue real-structure cases reproduce all graph-edge and pair-count numbers from full integer inverses and a second direct proof verifier.
- Python 3.11+ tests need SymPy **only** for generating untrusted integer proposals; the checker needs standard library only.

Reproduce after installing the repo:

```sh
python -m pip install -e '.[test]' sympy
python -m pytest -q tests/test_gnm_joint_rank.py tests/test_gnm_joint_rank_real.py
python experiments/real_1d3z_joint_rank.py --eps 0.025 0.05 0.1 0.2 0.3
python experiments/verify_1d3z_source.py  # networked public-source provenance check
```

## What this emphatically does not show

No new successful blind protein-folding experiment, improvement over GNM in *biological motion prediction*, native stability, kinetic prediction, or generalization across proteins. The NMR conformer spread is an underdetermined **structural-ensemble** quantity, not necessarily a thermodynamic fluctuation measurement. A preliminary descriptive comparison on the 10-model 1D3Z ensemble gave standard GNM Spearman ~0.27 against aligned Cα variance, similar to a simple inverse-contact-degree baseline (~0.27). This result should **not** be framed as predictive success; no model was fit to these labels.

This is a **mathematical certification and efficiency** gain on actual deposited molecular geometry. The next scientific gate is benchmarking the new certificate's utility on proteins with independent experimental motions/uncertainties and comparing robust coverage, false claims and runtime with other bounded model families. It would be premature to call it a physical protein-folding breakthrough.

The Python checker is not yet Lean-refinement-proved and is not a claim of original spectral mathematics. It offers exact, replayable arithmetic obligations that could be formalized later.
