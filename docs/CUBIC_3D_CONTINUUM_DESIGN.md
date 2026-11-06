# First independently audited 3D cubic continuum design reference in FoldBounds

**Scope:** An entirely additive, stacked research branch based on #57 → #54 → #52 → #48 → #47 → #44. No modification to existing 2D geometry, Lean, Python certificate checkers, CI, FoldBounds protein-motion code, RNase A blind experimental data, or other agents' branches. This is a finite 3D **lattice HP** model, not a realistic protein force field.

## 1. Cubic geometry with a distinct provenance type

The first polymer edge is anchored from (0,0,0) to (1,0,0); every subsequent step moves to an unoccupied lattice neighbor in **six** possible coordinate directions. Only nonconsecutive lattice-adjacent pairs of H residues are stabilizing contacts. Exact contact weight is q>=1, dimensionless; the Boltzmann partition is the weighted sum over **every anchored walk**.

The symmetry group fixing the oriented first bond contains the eight symmetries of the transverse yz-plane (D4). An orbit representative is uniquely selected by the first nonaxial turn toward +y, then the first departure from the xy plane toward +z. Its orbit weight is:

- **1** for the unique completely straight path;
- **4** for a nonstraight path lying in the xy plane;
- **8** when some vertex has nonzero z after the first transverse turn.

Each labelled path's symmetry orbit preserves its contact bitmask because lattice isometries preserve all inter-residue Manhattan distances. Therefore this reduction enumerates the complete fixed-first-edge 3D conformational universe, including mirror images. The Python builder incrementally records exact contact-bitmasks, and its immutable `CubicContactAtlas` class rejects malformed masks and multiplicities. These congruences are **not** a completeness certificate by themselves.

For independent validation, `audit_cubic_atlas` separately enumerates **every** 3D walk without symmetry, using a different direction order and recomputing all contact bits from full leaf paths. It must match the *entire contact-mask multiplicity table*, not only the total walk count. A separate language implementation `experiments/independent_cubic_audit.cpp` similarly enumerates all walks without symmetry and scores contacts from leaf paths, compiling via a standard C++17 compiler. It emits an exact two-column integer histogram sorted by bitmask.

Known exact anchored cubic walk totals, by number of vertices:

| n | Anchored 3D walks |
|---:|---:|
| 4 | 25 |
| 5 | 121 |
| 6 | 589 |
| 7 | 2,821 |
| 8 | 13,565 |
| 9 | 64,661 |
| 10 | 308,981 |
| 11 | **1,468,313** |
| 12 | **6,989,025** |

For **n=11**, the Python independent *unsymmetrized* auditor confirmed all 1,468,313 walks and **2,180** contact masks. `audit_cubic_landscape` independently reconstructed the complete exact `t=q-1` partition coefficient vector for **all 2,048 H/P sequences**, directly from each contact map's score and binomial coefficients, rather than using the builder's edge-subset expansion or subset-zeta transform.

For **n=12**, the independent C++17 unsymmetrized enumerator verified all **6,989,025** walks and all **8,738** contact-mask multiplicities, exactly matching the Python D4-reduced enumeration. The sequence-independent positive polynomial expansion covers **4,096 sequences** (12 residues) from that independently checked 3D geometry. The hosted test includes a C++ compiler-conditional full-histogram recheck; it skips only that cross-language test if no compiler is present.

## 2. Universal exact polynomial extension from 2D to 3D

The all-sequence graph-subset expansion from #47 is **dimension independent**: the contribution of a contact edge (i,j) is `1+t*x_i*x_j`, and all exact integer polynomial coefficients result from expanding products of contact edges. `build_cubic_landscape` uses the same abstract `UniversalLandscape` type but requires the distinct `CubicContactAtlas` input type and its **3D-specific independent audit by default**. It cannot call the existing 2D atlas audit by mistake. The downstream exact rational Sturm comparison works on the resulting polynomial table without requiring lattice dimension to enter its algebra.

`full_skyline_phase.py` is a deliberately conservative independent reference implementation of complete fixed-composition global optimization. Unlike the optional additional pruning in #54/#52, this checker **retains every coefficientwise-undominated polynomial class** and explicitly isolates all positive roots of all pairwise polynomial differences. It audits coverage of *every original sequence* against final coefficientwise survivors, orders root-free intervals using exact rational evaluations, and marks root events as real global switches, global tangencies, or irrelevant pair crossings. Shared/unseparated root brackets and resource overruns cause **fail-closed abstention**, not guessed results.

### All 3D compositions, across the entire real q>1 continuum

For n=11, exact continuum certificates cover **every one of 2,048 sequences** across 12 fixed-H composition classes and establish **11** global-optimum changes in total.

For n=12, exact continuum certificates cover **every one of 4,096 sequences** across 13 fixed-H composition classes, yielding **14** global-optimum changes. The per-composition count is:

| H residues k | Global switches for n=12 |
|---:|---:|
| 0 | 0 |
| 1 | 0 |
| 2 | 0 |
| 3 | 1 |
| 4 | 2 |
| 5 | 1 |
| 6 | 2 |
| 7 | 3 |
| 8 | 1 |
| 9 | 3 |
| 10 | 1 |
| 11 | 0 |
| 12 | 0 |
| **Total** | **14** |

These are mathematically exact *finite-model partition-function* claims conditional on source enumeration and exact Python computation. The rational root intervals prove the ranking across the full real q continuum; finite q witness checks are an independent regression, not an interpolation argument.

## 3. Dimension-dependent **reentrant** ranking — exact counterexample

At n=11, k=4 H residues, consider

    A = HPPHPPHPPHP
    B = HPPHPHPPHPP.

The complete, independently audited **2D** polynomial difference in t=q-1 is

    Z_A^(2D) - Z_B^(2D) = 2*t*(t^3 - 2*t^2 - 21*t + 269).

Its cubic factor has **zero positive real roots**, so A has strictly greater partition function than B for *every real q>1*.

For the independently audited cubic **3D** geometry the exact difference is

    Z_A^(3D) - Z_B^(3D) = 4*t*(51*t^3 - 314*t^2 - 3187*t + 9755).

The cubic factor has **exactly two distinct positive roots**. Thus the 3D ranking is

- `A>B` for `1<q<q1`,
- `B>A` for `q1<q<q2`,
- `A>B` for `q>q2`,

with q1≈**3.664164**, q2≈**11.397653**. The exact checker returns rational isolating brackets for each root and proves no further positive crossings.

More strongly, these two sequence classes actually exchange the **global fixed-composition maximum** in the n=11 cubic model: the 3D k=4 ensemble has exactly two global switches B↔A, and independent exhaustive evaluation at q=10 identifies **B and its reversal** as the only 3D global maximizers. In contrast, the **2D** global maximizers at q=10 are **A and its reversal**, with no overlapping sequence. This is a model-dependence counterexample, *not* proof that either toy model accurately describes experimental folding.

## 4. Run and trust boundaries

After checkout of the stacked branch, run:

```bash
python -m pytest tests/test_cubic_contact_atlas.py
python experiments/cubic_3d_phase_demo.py --length 11
python experiments/cubic_3d_phase_demo.py --length 12
g++ -O3 -std=c++17 experiments/independent_cubic_audit.cpp -o cubic_audit
./cubic_audit 12 > cubic12_histogram.tsv
```

Tests cover known 3D SAW counts; full independent Python geometry through n=11; all 2,048 exact n=11 sequence spectra and polynomials; independent C++ full n=12 histogram; exact n=11/n=12 composition global phase envelopes and exhaustive rational witness maxima; strictly model-dimension-dependent 2D/3D ordering; malformed inputs and resource-exhaustion failure modes.

**Model and originality limits:** Cubic HP is still a simplified coarse-grained polymer with only H/P labels, uniform contact energies, no solvent, rotamer states, hydrogen bonding, real temperatures, unfolded-state calibration, or kinetics. It is closer in dimensionality to a protein than a square lattice but **not a physics-validated folding model**. Exact SAW counts, HP designability, lattice symmetry, polynomial transforms, and Sturm root isolation are existing mathematical techniques and literature topics; no original theorem or scientific discovery is asserted. Both Python and C++ are independently compiled/executed research programs, **not Lean-kernel verified**. Existing Lean CI does not prove these new code paths. The meaningful next frontier is verified non-enumerative bounds in richer physical models and blind experimental validation.
