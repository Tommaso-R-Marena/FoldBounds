# Publication and prior-art gate: exact lattice-polymer thermodynamics

**Status:** Source-grounded research plan, not a completed external benchmark.
The references below concern related but sometimes **different** lattice
Hamiltonians. Do not describe their results as direct head-to-head FoldBounds
comparisons without the explicit model alignment described here.

## Directly relevant prior work

1. **Giaquinta and Pozzi (2013).** *An Effective Exact Algorithm and a New
   Upper Bound for the Number of Contacts in the Hydrophobic-Polar
   Two-Dimensional Lattice Model.* Journal of Computational Biology 20(8),
   593–609. https://doi.org/10.1089/cmb.2012.0266.
   Proposes improved upper bounds for attainable HP contacts and the
   minwalk exact optimization method, reporting benchmark optima on
   sequences of up to 80 residues. This is essential prior art for any
   claimed *new HP contact-score cap*. Their optimization objective is
   the ground state, not a full temperature-dependent partition function.

2. **Lee (2018).** *Transfer matrix algorithm for computing the exact
   partition function of a square lattice polymer.* Computer Physics
   Communications 228, 11–21.
   https://doi.org/10.1016/j.cpc.2018.03.022.
   Develops an exact transfer method for nearest-neighbor **interacting
   self-avoiding walks (ISAW)**, reaching up to 42 monomers, with a
   reported asymptotic computation-time base near 1.6 versus ~2.7 for
   straightforward enumeration. This is a powerful methodological
   baseline, not automatically an HP-sequence baseline.

3. **Hsieh, Chen and Hu (2016).** *Efficient algorithm for computing exact
   partition functions of lattice polymer models.* Computer Physics
   Communications 209, 27–33.
   https://doi.org/10.1016/j.cpc.2016.08.006.
   Reports exact ISAW partition functions up to 40 monomers on the
   square lattice and discusses possible extensions to the HP model.

4. **Jensen (2013).** *A new transfer-matrix algorithm for exact
   enumerations: self-avoiding walks on the square lattice.*
   https://arxiv.org/abs/1309.6709.
   Counts lattice self-avoiding walks up to length 79 using improved
   transfer-matrix enumeration. This is an important counting baseline,
   especially at q=1, not an interacting HP thermodynamics baseline.

5. **Bachmann and Janke (2003; arXiv posting 2007).** *Density of States for
   HP Lattice Proteins.* Acta Physica Polonica B 34.
   https://arxiv.org/abs/0710.4074.
   Provides exact density-of-states results for short HP sequences and
   multicanonical methods for longer chains. **Their displayed 14-residue
   example is three-dimensional**, as evidenced by the 13-step total
   943,974,510, which agrees with the cubic-lattice sequence OEIS
   A001412 rather than the square-lattice A001411. Their numerical
   histogram is therefore NOT directly comparable to our 2D HP counts;
   keep this paper as methodological HP thermodynamics prior art.

6. **Chen and Lin (2005).** *Partition function zeros of the
   two-dimensional HP model for protein folding.* Physica A 350, 45–51.
   https://doi.org/10.1016/j.physa.2004.11.035.
   An existing exact-enumeration study of 2D HP partition functions.
   This establishes that exact finite-HP thermodynamic calculations
   are established scientific work.

7. **OEIS A001411 / Iwan Jensen (2004).** Exact sequence of the
   number of n-step square-lattice self-avoiding walks.
   https://oeis.org/A001411. For FoldBounds' anchored first step,
   the expected count of N-monomer conformations is exactly
   A001411(N-1)/4. In particular, N=2 through N=10 should count
   1, 3, 9, 25, 71, 195, 543, 1479, 4067 walks respectively.
   This normalization still counts both mirror-image orientations.

## Critical apples-to-apples controls

Before quoting a numerical performance improvement, record in a machine-
readable manifest for **every** run:

- Number of monomers N versus number of bonds N-1; do not conflate.
- Square lattice in two dimensions versus simple cubic in three.
- Root anchored at (0,0) and first step fixed to (1,0).
- Whether translation/rotation/reflection/reversal are factored out.
  FoldBounds fixes the first bond and does **not** quotient reflections.
- Whether contacts include only nonconsecutive H-H adjacency (HP),
  or **all** nonconsecutive monomer adjacency (homogeneous ISAW).
- The sequence, H density, contact threshold, and any constraints on
  end points or compactness.
- The interaction factor q; report contact polynomial coefficients
  when possible rather than relying on rounded exp(beta*epsilon).
- Interval width, exact benchmark oracle where tractable, elapsed CPU
  time, peak memory, certificate size, verification time and producer
  time separately.
- Host architecture and number of threads; one branch expansion is not
  equivalent to one unit of wall-clock time.

**Valid special case:** for an all-H sequence, FoldBounds' H-H contact
count agrees in physical definition with ISAW all-contact counting on
the same grid and the same rooted walk space. The path-counting and
normalization conventions must still match.

**Invalid general comparison:** comparing FoldBounds' arbitrary mixed
H/P sequence partition function to Lee/Hsieh homogeneous ISAW results
and claiming an accuracy or runtime improvement.

A fast *optimizer* for the maximum contact count does not automatically
compute or tightly bound the weighted partition sum over all folds.
Keep HP ground-state optimization and finite-temperature thermodynamics
as separate benchmark categories.

## Required baseline experiments

**Gate A: exact small-model validation.** Independent oracle and Lean/Python
differential tests for 2 <= N <= 10 with varied H/P patterns; compare
the exact contact histogram and Z(q), q in {1, 3/2, 2}. Test incorrect
symmetry normalization as a negative control.

**Gate B: finite-temperature 2D HP baselines.** Prefer Chen and Lin
(2005) and other published **square-lattice** HP densities, after
checking exact sequence, symmetry and orientation conventions.
Compare contact polynomials rather than only a floating-temperature
probability. The Bachmann/Janke 14-mer 3D table is an explicit
negative control and must not be presented as a 2D baseline.

**Gate C: all-H ISAW baseline.** Compare contact-count generating
functions against Lee/Hsieh where table/data/code can be obtained.
Normalize counts by the exact anchoring convention; do not guess or
assume a symmetry multiplicity factor.

**Gate D: useful non-enumerative certificates.** For N in {10, 12, 14,
16, 18} and varied sequences/events, report sound interval widths at
controlled CPU times, together with certificate construction time,
verified checker time and certificate bytes. Explicitly record
failure cases with near-vacuous intervals.

**Gate E: formal novelty.** Before claiming publishability, separately
establish (i) a theorem not already implied by established results
under a trivial renaming, (ii) an algorithmic advantage under a matched
Hamiltonian and computational budget, or (iii) a real proof-carrying
assurance contribution that validates compact scientific computations
not previously addressed by existing methods.

No benchmark or theorem novelty is claimed from this document alone.

## Relevance to current Aristotle work

The highest-value independent assurance advance is a small executable
Lean checker of a **compact subtree or generating-function certificate**
that proves event mass and probability bounds without scanning every
conformation. A polynomial-time verifier is a distinct scientific
contribution from solving protein folding in polynomial time.

The literature makes a clear warning: matching tiny exact enumeration
results cannot establish a major algorithmic breakthrough. The research
case should rest on proof-carrying guarantees, their verification cost,
and independently measured bound quality.
