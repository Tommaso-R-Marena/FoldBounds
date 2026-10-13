# FoldBounds

**Certified Protein-Folding Thermodynamics — mathematical research prototype**

FoldBounds studies rigorous algorithms for partition-function and
conformational-population bounds in explicitly specified statistical-mechanical
models. The long-term research challenge is to compute informative guarantees
without exhaustively enumerating every conformation.

**Current status (v0.1):** a conservative, fully rational Python reference
implementation for the 2D square-lattice hydrophobic/polar (HP) model.
It enumerates exact contact spectra for small sequences and checks the
completeness and geometric validity of partial-enumeration certificates.

**Trust boundary:** Guarantees are relative to this simplified finite model.
The Python implementation is **not yet formally verified by Lean** and does not
predict experimental structures of real proteins. No novelty or publication
claim is made at this stage.

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
