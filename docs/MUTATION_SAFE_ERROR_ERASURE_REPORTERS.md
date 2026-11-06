# Mutation-site invariant, erasure-aware GA/GB geometric codes

## Scientific question

Can the **original 2009** GA95/GB95 structural distance-reporter
fingerprints be redesigned so every reporter site has exactly
**the same amino acid** throughout the full 2012 one-residue
fold-switch cascade, while also handling missing readouts and
arbitrary channel faults in a deterministic, fail-closed manner?

The existing PR #81 correctly recognized all 70 deposited external
NMR structural conformers using the previously selected code,
but several reference codes used **site 25**, which mutates T→I
in 2LHG/2LHE. Thus those codes were mathematically useful,
but not mutation-invariant biochemical labeling suggestions.

This new branch is **strictly stacked** on PR #81 and adds
independent modules. It modifies no previously reviewed code or
other agents' branches.

## Data source, model domain and non-leakage

Only experimental **2KDL GA95 and 2KDM GB95**, twenty NMR
models each, train the code. The four 2012 experimental entries
2LHC GA98 (20), 2LHD GB98 (20),
2LHG GB98-T25I (10), 2LHE GB98-T25I,L20A (20)
are **test only**.

The *mutated-site list* [20,25,30,45] is known from the published
mutation cascade BEFORE fitting. Any reporter that uses these
sites is disallowed, including site 25 that affected prior tests.
All allowed reporter sites have identical amino-acid identities
across all six protein variants, confirmed in the code. This
is site-sequence conservation, **not** validation of label chemistry
or side-chain accessibility.

All six raw PDB files must exactly match previously authenticated
SHA-256 bytes. Every experimental sequence identity, NMR conformer
count and complete CA trace are revalidated.

## Primal/cover optimality certificates

Training-derived candidate reporter edges are restricted to
positions 1..14 paired with 15..56, residue separation at least
4, model interval distances 12..40 Angstroms, strictly disjoint
GA/GB coordinate distance bands, and no mutated site.

For each odd code length K, the checker maximizes minimum
distance-band gap g across **K residue-disjoint reporter edges**.
A witness comprises:

- K explicit chosen disjoint edges, all with gap≥g.
- A set of fewer than K residue vertices meeting every edge
  with strictly larger gap; this gives an independent finite
  optimality obstruction.
- Exact integer milli-Angstrom lower/upper bands for every
  chosen edge, reproducible only from the official reference PDBs.

The receiver independently regenerates every admissible edge,
checks each reported interval and computes the strict
larger-margin cover. It rejects forged, missing, duplicated
or mutation-site-bearing evidence.

As in PR #74, the certificate proves optimality **within this
explicit bipartite candidate family**, not among all possible
protein measurements or labeling sites.

## New fault model: erasures and erroneous channels jointly

For K pairwise residue-disjoint reporter channels, E unavailable
channels (perhaps a failed residue label) and T arbitrary **wrong**
remaining reporter votes, binary majority has at least
K−E−T correct votes and at most T wrong votes. Therefore the
true fold wins with a strict majority if

    2*T+E < K.

This is the classical sharp repetition-code error/erasure bound.
For example, a 5-channel code tolerates one wrong binary decision
plus TWO missing channels (2*1+2=4<5). The 7-channel
code tolerates two wrong decisions and two missing channels
(2*2+2=6<7).

The residue-disjoint matching ensures each independently missing
residue can disable **no more than one channel**. Thus up to
E missing sites cause no more than E missing reporter channels.
Arbitrary chemical or geometric correlation in the resulting
incorrect votes must still be counted inside the T budget.

The decoder does **not** invent a class when the surviving
vote counts tie or measurements land exactly on a boundary;
it returns an explicit abstention. It rejects untrusted channels
and non-exact readouts. The bound is **conditional** on every
other channel correctly representing the relevant source-defined
class; no experimental noise calibration is assumed.

## Geometry-based exact experimental verification

On every one of 70 heldout deposited NMR conformers, the
runner computes precise enclosing integer distance brackets
for every selected reporter. It evaluates all erasure subsets
of up to 2 channels and, for each, every wrong-vote subset
obeying 2*T+E<K. Every classification success and failure,
per-model, is retained.

It **separately reports**:
- Clean structural fold recognition and abstentions.
- Correct combined error/erasure patterns.
- The smallest truth-labeled signed distance to a fixed threshold.
- Every channel lying *outside* the 2009 training conformer
  distance envelope. An unseen conformer may classify correctly
  despite not belonging to the training energy/geometry envelope,
  so this is a critical out-of-distribution warning.

No heldout conformer can change a reporter pair, optimize a
threshold or be omitted from the audit.

The new Lean \`FoldBounds/MutationSafeCodes.lean\` kernel-checks
six precise counting/margin theorems, including the
sharp error/erasure bound, tie counterexample at equality,
lost-site to channel-erasure budget implication, and combined
pre-existing error budget. The complete 3D atom geometry,
sequence-to-structure mapping and Python-to-Lean input
decoder remain separate audited Python obligations.

## Limits and potential importance

A CA–CA geometric distance is not automatically the distance
between fluorophores in a real FRET measurement.
"Mutation-invariant residue" does **not** guarantee
label feasibility, solvent exposure, zero fold perturbation,
or measured calibration. The six deposited structures are
not six independently sampled real biological experiments;
the 70 conformers are not 70 independent observations.

This is a rigorous, experimentally sourced **geometry-based
measurement-design benchmark** and a potentially useful
proof-carrying experimental planning component—not
a validated protein folding physics breakthrough,
a prospective blind sequence-to-fold predictor,
a thermodynamic population model, or a kinetics predictor.

The true test outcomes are deliberately **not predeclared**
here. They must emerge from successful immutable-source
GitHub Actions on the new branch, and all failures remain
visible in the committed reproducible JSON.

Reproduce:

    python -m pytest -q tests/test_ga_gb_mutation_safe.py
    python -m experiments.mutation_safe_real_audit \
      --2kdl <2KDL.pdb> --2kdm <2KDM.pdb> \
      --2lhc <2LHC.pdb> --2lhd <2LHD.pdb> \
      --2lhg <2LHG.pdb> --2lhe <2LHE.pdb> \
      --out new_safe_results.json

The dedicated
\`.github/workflows/mutation-safe-codes.yml\`
performs all steps using the previously verified exact SHA-256
coordinates and preserves the unfiltered JSON output.
