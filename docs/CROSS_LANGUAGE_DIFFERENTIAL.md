# Python/Lean HP lattice differential testing

This test is an independent, *finite regression comparison*, not a proof
that an arbitrary Python computation refines a Lean specification.

The test driver \`scripts/cross_language_diff.py\` independently enumerates
canonical 2D self-avoiding walks with first edge (0,0) to (1,0), **without
using the production Python walk generator**. It computes contact scores
from direct nonconsecutive H-H pair adjacency, **without using production
Python's parity-filtered contact-pair algorithm**. It independently checks
production Python's path validator, contact scorer, and, through length 9,
the complete exact contact-count spectrum.

It then generates Lean fixtures for a deterministic sample of complete
walks at lengths 2 through 10, with all-H, all-P, alternating and mixed
H/P sequences. The pinned Lean executable evaluates:

* \`walkFrom\` path decoding from four absolute directions;
* \`contactCount\` scores on the same explicit coordinates;
* exact rational weights \`(3/2)^score\`.

The Python command **fails** unless every group outputs its expected
pass marker, and rejects Lean evaluation errors and mismatch markers.

## Reproduction

After installing Python 3.11+ and pinned Lean 4/Mathlib:

    PYTHONPATH=src python3 scripts/cross_language_diff.py

For only the independent Python checks:

    PYTHONPATH=src python3 scripts/cross_language_diff.py --python-only

The full test runs as a required step in the existing Lean CI workflow.

## Limits

1. These tests **do not establish universal Python/Lean equivalence**.
   The full bridge requires a mathematical proof of the canonical
   representation correspondence and a verified parser or byte-level
   schema, preferably with a Lean-native witness checker.
2. Most Lean tests are runtime evaluations of explicit finite fixtures,
   not theorem declarations proved by kernel reduction.
3. The finite test does not independently prove Python's maximum-flow,
   isoperimetric, witness search or weighted transfer relaxations.
4. The rational q model remains distinct from physical energy models.
5. Fixture sizes and samples are fixed for reproducibility; do not
   claim statistical inference or large-system performance from this run.

The goal is to quickly identify discrepancies while Aristotle develops
the full Lean representation and executable-refinement theorems.

## External published-count control

In addition to cross-language comparisons, the independent enumerator
and production Python spectrum must match OEIS A001411's published
square-lattice self-avoiding-walk counts. For an N-monomer rooted walk
with the first direction fixed, the count is `c_(N-1)/4`, where c_k is
the published number of unrestricted k-step walks. This division removes
only four initial rotations, **not reflection symmetry**.

Reference: https://oeis.org/A001411 (Iwan Jensen and preceding work).
The test rejects any disagreement with the published counts, so shared
implementation mistakes in both custom enumerators have an additional
external numerical check. This remains finite-sample validation.
