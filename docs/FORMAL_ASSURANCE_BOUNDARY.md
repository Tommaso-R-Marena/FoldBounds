# Formal assurance boundary and integration map

## Status

Aristotle run `8becb5ac-4678-49a8-8c98-84d3faa3f7b7` supplies a
Lean 4 / Mathlib **abstract finite-conformation thermodynamics** library.
The sources are imported verbatim except for removing an unused
RequestProject entry in `lakefile.toml`. No computational-lattice
implementation is contained in these theorems. This document is an
engineering interface specification, not a theorem.

## Mathematical meaning

Lean models a finite type `C`, energy `E : C → ℝ`, inverse temperature
`β : ℝ`, weight `exp (-β * E c)`, partition function `Z`, and event
probability `ZOn/Z`. Abstract proofs handle exact real numbers and
supply interval inequalities under stated assumptions.

Python models **anchored 2D square-lattice H/P self-avoiding walks** with
energy `E(c) = -ε k(c)`, where `k(c)` is the count of nonconsecutive
H-H contacts. The reference implementation accepts an exact rational
`q >= 1` and sums `q^k`. For the physical Boltzmann interpretation,
one would need `q = exp (β ε)`, and that equality generally cannot be
assumed for arbitrary rational `q`. The rational Python model is a
well-defined combinatorial model in its own right.

## What Lean currently guarantees

- Positivity and normalization of abstract finite-space Boltzmann
  partition functions and probabilities.
- Validity and optimality of certain abstract probability enclosures
  **given** sound tail/event/complement and block energy bounds.
- Abstract monotonic refinement, convergence-width, free-energy,
  compositional, symmetry, and Gibbs–Bogoliubov theorems under explicit
  hypotheses. See [original theorem index](THEOREMS.md) and
  [Aristotle's technical report](ARISTOTLE_FORMAL_REPORT.md).
- A 34-declaration axiom audit uses only `propext`,
  `Classical.choice`, and `Quot.sound`, checked by CI.

## What is NOT proved yet

- That the Python prefix-frontier generator enumerates precisely the
  required canonical conformation space, without missing or duplicating
  folds. Reflection convention must be explicit.
- That the Python parity, reachability, maximum-flow, isoperimetric,
  transfer-DP, and witness constructions produce valid certified
  **energy or weighted-mass bounds** for the exact Lean definitions.
- That serialized input is parsed without semantic ambiguity or that
  Python's implementation refines a proved Lean decision function.
- That rational `q^k` equals Lean's transcendental exponential
  weights for an arbitrary modeled physical temperature.
- That any algorithm improves asymptotic complexity or accurately
  models real proteins or experimental thermodynamics.

The package is therefore **formally proved mathematics alongside a
tested computational reference**, not a formally verified Python
protein-folding predictor.

## Next proof obligations

1. Define canonical self-avoiding lattice conformations as an explicit
   `Fintype` and formally define the H-H contact energy.
2. Prove frontier partition/coverage, no duplicates and independently
   checked exact-witness lower mass for that type.
3. Prove the isoperimetric, graph b-matching, reachability and transfer
   relaxations each bound the true conformational weighted mass.
4. Establish a sound `q`-parameterized rational formulation, or a
   verified rational interval enclosure for `exp (β ε)`, and relate
   either to the real-energy theorems.
5. Build an independently executable certificate format whose accepted
   results refine the formal definitions, explicitly scoping the TCB.
6. Benchmark against published exact/combinatorial methods with matched
   Hamiltonians, symmetry conventions, input sets and compute budgets.

## Reproducibility

```bash
python -m pip install -e '.[test]'
python -m pytest
lake build
lake env lean FoldBounds/AxiomAudit.lean
python3 scripts/check_lean_sources.py
```

Python checks and Lean proofs run in separate CI jobs. A successful job
does not establish an implication between the two implementations.
