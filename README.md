# FoldBounds

**Certified Protein-Folding Thermodynamics — mathematical research prototype**

FoldBounds studies rigorous algorithms for partition-function and
conformational-population bounds in explicitly specified statistical-mechanical
models. The long-term research challenge is to compute informative guarantees
without exhaustively enumerating every conformation.

**Current status:** a Lean 4 / Mathlib formal thermodynamics library
with an executable finite rational certificate kernel, alongside a
tested, fully rational Python research prototype for the 2D square-lattice
hydrophobic/polar (HP) model. Lean proves abstract thermodynamic bounds
and the correctness of the new finite-universe checker. Python enumerates
small contact spectra and checks model-relative frontier certificates.
**The Python implementation is not yet formally refined to the Lean kernel.**

**Trust boundary:** Guarantees are relative to this simplified finite model.
The Python implementation is **not yet formally verified by Lean** and does not
predict experimental structures of real proteins. No novelty or publication
claim is made at this stage.

## Formal verification

Pinned to Lean 4.28.0 and Mathlib 4.28.0:

    lake build
    lake env lean FoldBounds/AxiomAudit.lean
    lake env lean FoldBounds/ExecutableTests.lean

The [Lean anchored HP model](docs/LEAN_LATTICE_MODEL.md) defines
self-avoiding lattice conformations, contact scores and rational weights
inside Lean itself. Its finite-universe certificate checker does NOT yet
prove the separate Python search algorithms correct.

The [Lean executable certificate checker](docs/EXECUTABLE_CERTIFICATE_KERNEL.md)
proves that accepted finite rational witness/cap certificates enclose
event masses and their coupled probability ratio. The rational finite
universe is currently enumerated during upper-cap verification; this is
**not** an efficient, proved lattice-folding algorithm or a Lean verification
of the Python code.

## Reproduce

Requires Python 3.11+:

    python -m pip install -e '.[test]'
    python -m pytest
    python -m foldbounds HHHH --q 2 --depth 4 --event-at-least 1 --exact-spectrum
    python -m foldbounds HHHHHHHH --q 2 --depth 2 --event-at-least 1 --enhanced-lookahead 4
    python experiments/compare_bounds.py
    python experiments/compare_isoperimetric.py
    python experiments/compare_adaptive.py > adaptive_comparison.csv
    python experiments/compare_transfer.py > transfer_comparison.csv
    python experiments/compare_witnesses.py > witness_comparison.csv
    python experiments/scaling_witnesses.py > scaling_witnesses.csv
    python -m foldbounds HHHHHHHH --q 2 --enhanced-lookahead 0 --matching --transfer
    python -m foldbounds HHHHHHHH --q 2 --adaptive-steps 12 --adaptive-policy probability --matching
    python -m foldbounds HHHHHHHH --q 2 --depth 2 --enhanced-lookahead 2 --matching
    python -m foldbounds HHHHHHHH --q 3/2 --depth 4 --event-at-least 1

The optional enhanced lookahead checker uses geometric upper bounds, explicit
full-length conformation witnesses, and nested one-step subtree refinements.
See [proof sketch and limitations](docs/BOUND_HIERARCHY.md) and the
[maximum-flow matching relaxation](docs/MATCHING_BOUND.md). This is an
internally tested research method, not a formal Lean proof or claimed novelty.

An optional [isoperimetric upper bound](docs/ISOPERIMETRIC_BOUND.md)
limits all possible H-H contacts using the perimeter and bounding box
of the current partial conformation. Its mathematical assumptions are
explicitly restricted to the finite 2D square-lattice model.

The [adaptive search algorithm](docs/ADAPTIVE_CERTIFICATES.md) optionally
refines promising branches while independently rechecking the
resulting interval bounds. Its auditable history retains previous checked
frontiers to ensure reported intervals never widen. It also provides checked zero-completion terminal
certificates for trapped conformations. These are software experiments, not
claimed discoveries.

An optional [weighted prefix-obstacle transfer relaxation](docs/TRANSFER_RELAXATION.md)
upper-bounds the Boltzmann mass with an exact-integer dynamic program
that includes every legal completion within a relaxed larger set.
It remains a tested 2D lattice model, not a Lean-certified predictor.

A [checked complete-conformation witness mechanism](docs/CHECKED_WITNESSES.md)
lets an untrusted bounded DFS propose full self-avoiding folds. An
independent verifier rejects malformed/duplicate paths and derives
rigorous event/complement *lower* masses by exact rational summation.
Use `--witness-nodes 256` with `--enhanced-lookahead 0` to enable it.
This is computational model-relative assurance, not a Lean theorem.
See [the controlled witness experiment log](docs/WITNESS_EXPERIMENT_LOG.md)
for small-system successes, failures, and explicit limits.

The parameter q is an exact rational Boltzmann factor exp(beta * epsilon) >= 1
for the model energy E = -epsilon * (number of nonconsecutive H-H lattice contacts).
A reported probability interval bounds the event that contact score is at least
the chosen event threshold. A depth equal to sequence length produces an exact
rational result; smaller depths return conservative bounds without exploring
all conformations. Broad [0,1] intervals are expected for shallow trees.

## Research plan

1. Establish foundational partition, probability, and free-energy theorems in Lean 4.
2. Build an independently checkable finite-model executable specification.
3. Improve admissible bounds and search strategies with rigorous comparisons
   against exact enumeration and published baselines.
4. Extend to richer, explicitly qualified coarse-grained energy models.

See [soundness and trust boundaries](docs/SOUNDNESS.md) for the mathematical
model, certificate argument, limitations, and pending Lean bridge.

## Collaboration and reproducibility

Clearly distinguish mathematical theorems, software tests, physical assumptions,
and speculative research claims. Benchmarks should record exact parameters,
published baselines, computational budgets, and failed experiments.

Licensed under Apache 2.0.


## Formal thermodynamics — Lean 4 / Mathlib

The [Aristotle formal library](docs/ARISTOTLE_FORMAL_REPORT.md) studies finite
Boltzmann partition functions, event probabilities, coupled optimal enclosures,
convergence criteria, symmetry, and free-energy bounds.

- Pinned Lean 4 / Mathlib 4.28.0, with separate formal CI and a 34-theorem
  axiom audit. See the original [run summary](docs/ARISTOTLE_RUN_SUMMARY.md).
- **Trust boundary:** Lean proves abstract real-number theorems with explicit
  model/energy certificates as hypotheses. It does **not** currently verify
  the Python lattice algorithms, certificate parser, numerical exponentials,
  molecular physical model, or practical runtime performance.
- The imported original Lean sources receive no algorithmic edits; only
  the obsolete RequestProject Lake template entry is removed.
