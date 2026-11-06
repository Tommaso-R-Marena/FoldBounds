# Exact global minimax repair brackets for two real processed NOE restraint sets

**Status:** additive experimental FoldBounds work, based on PR #80. Exact
arbitrary-integer certificates for a *selected* set of atom-level constraints.
Not a Lean-refined algorithm, chemically constrained structure refinement,
new biological finding, new mathematical priority claim, or a blind protein
prediction. Do not merge without source and chemical-assignment review.

## Objective and assumption boundary

For a processed NMR coordinate model with original atom positions `X_a` and
one-row, unambiguously assigned interval restraints `(a,b,L,U)`, define

    rho* := inf { max_a ||Y_a - X_a|| :
                  L_ab <= ||Y_a - Y_b|| <= U_ab for EVERY selected record }.

This is a **joint** nonconvex feasibility optimization, unlike PR #80's
individual-pair violation flags. The `Y_a` are unconstrained 3D locations:
**no covalent bond/angle, torsion, stereochemistry, clash or excluded NOE
constraint is imposed**. Atom identifiers refer to the PDBStat-converted
NEF+AMBER-PDB files, not unprocessed spectra.

The exact integer source data are inherited from PR #80:
`data/pinned_unambiguous_noe_2png_6nbn.json` (Git blob
`96daf4e88743df4e623cdf5495f28b8e2113fc32`). That frozen parser
selects only 304 of 1,732 original 2PNG restraint IDs and 553 of 1,694
6NBN IDs. It excludes mixed assignments, pseudoatom/wildcard groups,
missing hydrogens and unsupported combination semantics. The selected graph
has respectively 170 and 381 mapped atoms; the actual protein comprises
more atoms. *A repaired subset is not a repaired NMR structure.*

## Mathematical construction: sharp geometric dual + full primal witness

For any one original observed pair `(a,b)` with distance `d`, lower L and
upper U, displacement at most `rho` per endpoint implies

    | ||Y_a-Y_b|| - d | <= ||Y_a-X_a||+||Y_b-X_b|| <= 2*rho.

Therefore if the original pair violates U, `rho* >= (d-U)/2`. If it
violates L, `rho* >= (L-d)/2`. A source restraint whose violation is
strictly too large for a proposed radius gives a valid, global lower
certificate regardless of other constraints. This inequality is standard,
not claimed novel.

The producer applies alternating Euclidean pair projections at floating-point
precision, then rounds **every change** to integer ticks of 0.000001 Å.
It stores only changed-atom displacement vectors. The independently written
verifier never trusts the optimization procedure: it recomputes all original
atom coordinates, applies every displacement and checks **all** constraints
using arbitrary-precision integer squared distances. It also checks every
atom's movement norm lies within the claimed upper radius. This furnishes a
feasible joint 3D arrangement and proves `rho* <= upper`.

For the lower certificate, let `S=1000` fine ticks per source milliångström.
If original squared distance is the integer `s`, then violating the upper
constraint even after radius `r` is exactly

    s*S*S > (U*S + 2*r)^2.

For a lower-bound violation the corresponding exact inequality is

    L*S > 2*r and s*S*S < (L*S - 2*r)^2.

A verified strict inequality proves `rho* > r`. **No floating-point
roundoff is used to accept the certificate.** Combined with the complete
integer candidate, the proof is a near-optimal global minimax interval.

## Reproducible results on real, processed NMR data

| Model | Selected restraints | Mapped atoms | Moved atoms | Lower witness ID | Certified optimal minimum max atom displacement |
|---|---:|---:|---:|---:|---|
| 2PNG | 304 | 170 | 10 | 1701 | **0.027283 Å < rho* <= 0.027284 Å** |
| 6NBN | 553 | 381 | 9 | 470 | **0.368720 Å < rho* <= 0.368722 Å** |

The numerical bracket widths are 0.000001 Å and 0.000002 Å, respectively,
*as exact mathematical statements about frozen 0.001 Å rounded model
coordinates and chosen interval endpoints*. They are **not** statements of
experimental accuracy to a millionth of an ångström. The large 6NBN repair
is dominated by the selected converted Val69 H–HB upper-distance conflict.
These precise comparisons may change under alternate atom assignments,
original raw-NOE interpretation, stereochemical restrictions or different
coordinate models.

## Negative result: frozen five-fold source-only restraint holdout

We also test whether unambiguous converted NMR **upper** restraints alone
can imply withheld restraint intervals. A fixed FNV-1a hash of
`<sample>:<restraint ID>`, modulo five, assigns each original restraint to
exactly one test fold. No PDB coordinates or withheld limit enters the
shortest-path producer. A separately replayed path certifies each inferred
upper bound (triangle inequality) on withheld atom pairs.

| Sample | Withheld exactly once | Pair connected by training graph | Conditional upper <=8 Å | Training path upper strictly below withheld listed upper | Complete withheld interval guaranteed |
|---|---:|---:|---:|---:|---:|
| 2PNG | 304 | 257 | 93 | 3 | **0** |
| 6NBN | 553 | 348 | 87 | 2 | **0** |

The **180 coarse <=8 Å upper bounds** (93+87) are logical inferences from
training upper constraints, *not newly measured Cα contacts* and not
validated on an independent protein. The 857 holdouts came from the same
processed experiment. The **zero** exact-interval result is a meaningful
negative control: the chosen upper-only graph cannot recover the strong
positive lower limits and typically cannot match measured upper limits.
The complete lower+upper NOE system, chemical bonds and ambiguous NOE
semantics have NOT been tested in that held-out protocol.

A source-only upper graph admits a collapsed one-dimensional embedding of
all atoms and hence cannot prove a strictly positive lower distance.
The previous PR #80 primal-dual shortest-path construction also establishes
that shortest-path upper bounds cannot be tightened from **upper-only**
constraints. This makes the modest holdout performance unsurprising, but
quantifies it without selecting examples after examining the scores.

## Reproduction and verification

    python -m pytest tests/test_noe_minimax_repair.py -q
    python experiments/noe_minimax_repair_demo.py

The tests recompute the Git blob SHA of the inherited original selected
fixture, reject forged coordinates, source IDs, movement bounds, absent
witness edges, and mutated primal displacement vectors. All live GitHub
workflows must pass before the PR could be considered for merge.

## Next scientific gate

The next strong test must incorporate raw/ambiguity-aware NOE semantics,
peptide chemistry and a genuinely held-out experimental assay. In those
conditions the current near-optimal unconstrained repair certificate is
not a prediction; it is an honest numerical baseline and an exact conflict
witness for a carefully delimited mathematical problem.
