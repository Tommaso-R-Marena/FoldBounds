# Prior work, comparison requirements and novelty gate

**Status:** FoldBounds is currently an exact-rational reference implementation
for a simplified 2D HP lattice model. We do not currently claim a new
scientific theorem, superior method versus published research, or physical
predictive accuracy for real proteins.

## Established prior work

- Iwan Jensen (2004), *Enumeration of self-avoiding walks on the square
  lattice*, Journal of Physics A 37, 5503.
  https://doi.org/10.1088/0305-4470/37/21/002
  Establishes advanced exact square-lattice walk enumeration.

- Iwan Jensen (2013), *A new transfer-matrix algorithm for exact
  enumerations: self-avoiding walks on the square lattice*.
  https://arxiv.org/abs/1309.6709
  Demonstrates transfer-matrix methods outperform elementary enumeration
  for suitable sizes.

- Emanuele Giaquinta and Laura Pozzi (2013), *An effective exact algorithm
  and a new upper bound for the number of contacts in the hydrophobic-polar
  two-dimensional lattice model*, Journal of Computational Biology 20,
  593-609. https://doi.org/10.1089/cmb.2012.0266
  Important direct prior art on H-H contact bounds and exact optimization;
  parity and H coordination alone are not research novelty.

- Julian Lee (2018), *Transfer matrix algorithm for computing the exact
  partition function of a square lattice polymer*, Computer Physics
  Communications 228, 11-21.
  https://doi.org/10.1016/j.cpc.2018.03.022
  Exact partition-function algorithms already exist. Its interaction
  Hamiltonian is not necessarily identical to our H/P sequence model;
  comparison requires precisely aligned energy definitions.

- *Exact enumeration of conformations for two and three dimensional
  lattice proteins* (2016), Computer Physics Communications 199, 8-11.
  https://doi.org/10.1016/j.cpc.2015.09.017
  Provides additional exact-lattice-enumeration prior art.

## Conditions for a research claim

1. Identify a contribution missing from published results: a new theorem,
   provably stronger bound, materially faster method, or a genuinely
   verified executable checker and its trusted computational boundary.
2. Match state-space and physical conventions: anchoring, rotations,
   reflection multiplicity, contact definitions and Boltzmann weights.
3. Compare with relevant published alternatives under matched parameters,
   work budgets and evaluation hardware.
4. Record full benchmark distributions, adversarial cases, failure cases,
   precision guarantees and scalability; not just selected success stories.
5. Distinguish rigor in the toy model from validation of real molecular
   physics, which requires a separate evidence/approximation bridge.

## Candidate next investigations (NOT accomplished discoveries)

- State-dependent transfer-matrix relaxations for rigorous weighted
  partition-function bounds.
- A small Lean-refined executable checker for frontier mass certificates.
- Adaptive allocation of exploration effort to maximize reduction in
  guaranteed probability interval width per unit of computational work.
- More expressive, explicitly specified coarse-grained molecular models.

The current internal speedups and exhaustive test results establish an
engineering baseline, not manuscript-ready novelty.
