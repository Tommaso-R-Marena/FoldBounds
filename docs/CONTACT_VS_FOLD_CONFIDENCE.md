# Same-model dual certificate: dominant contact event, negligible single-fold confidence

## Why this matters

FoldBounds' 3D HP contact-mass checker can establish that an enormous
set of compact, high-contact conformations accounts for most of the
equilibrium ensemble. That does **NOT** show any one conformation is
preferred, unique, or even likely.

The new \`contact_vs_fold.py\` checker requires the *exact same*
sequence, complete cubic target, contact fugacity, correction interval,
and unrestricted rooted 3D SAW ensemble for both independent statements:

1. \`P(K_HH >= t) >= L\`, by the positive event-mass lower family and
   universal upper bound on non-event mass in \`cubic_tiling.py\`.
2. \`P(one exact target) <= U\`, by explicitly verified spatially
   independent local rewires and a rigorous relative mass denominator
   in \`local_rewire_confidence.py\`.

It exposes both proven endpoints together, not a point estimate.

## Fail-closed assumptions

The coupled checker refuses to combine incompatible models. It requires:
- One sequence consisting only of H and P.
- A complete, canonical, validated 3D lattice path that occupies
  **all sites** of the claimed \(2n_x\times2n_y\times2n_z\) solid.
- Positive exact rational HH fugacity q >= 1, HP=PP=1, and all
  requisite pairwise weights explicitly supplied.
- The same exact positive bounds on arbitrary per-state correction
  multipliers in both certificates.
- An integer contact threshold t>=1, with the specified target
  independently verified to have at least t nonconsecutive HH contacts.

It **does not** silently compare a real-protein interaction matrix
with a simple-cubic HP Hamiltonian; the general structural checker
can use arbitrary residue pairs, but this combined certificate
explicitly excludes those different Hamiltonians.

The exact target belongs to the event, yet its probability can
still be extremely small because the event includes exponentially
many other different conformations.

## Independently reproducible examples

All statewise correction factors lie anywhere in [1/2,2] with no
independence assumptions.

| Synthetic model | Event | Contact fugacity | Event probability lower | Exact target probability upper |
|---|---|---:|---:|---:|
| 216 all-H residues, full 6³ cube | K_HH>=108 | 4 | >99% | < 10⁻⁸ |
| 216 residues, 20 P at i=0,11,... | K_HH>=54 | 5 | >99% | < 10⁻⁸ |
| 512 all-H residues, full 8³ cube | K_HH>=256 | 4 | >99% | < 10⁻¹⁸ |

The 216-all-H target also saturates the classical sharp cubic
edge-isoperimetric bound and is a global energy minimum for q>1;
nevertheless it has a tiny *individual* equilibrium probability.
This is an explicit mathematical example where **a guaranteed
majority high-contact event coexists with exclusion of a single
majority fold**.

None of these numbers are estimates for an experimental protein,
native fold, atomistic force field, or real folding temperature.

## Independent validation

    python -m pytest -q tests/test_contact_vs_fold.py
    python experiments/contact_vs_fold_demo.py

The tests independently compute contact counts directly by checking
all pairs of indexed residues, compare with the classical cubic
ground-state calculation, and confirm both quantitative bounds
under the identical declared HP model. Negative controls reject
conflicting energy matrices, coordinate boxes, noninteger
thresholds, missing/incorrect target membership, and non-HP
residue labels.

## Assurance boundary

The contact-mass inference and target-confidence mass inference have
generic Lean theorems elsewhere in this branch/its parent PR. Their
actual 3D coordinate scoring, tile enumeration, geometric injection,
and model identity checks are currently implemented in Python and
tested against independent finite exact oracles, **not end-to-end
Lean-kernel verified**.

This model-level diagnostic is useful for scientific claim hygiene:
**certified collapse/compactness is not certified structure
identification**. The distinction is familiar in statistical
mechanics, so we do not claim a new discovery in protein folding.
