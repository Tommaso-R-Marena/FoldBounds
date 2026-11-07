# Scalable local JEPA → exact rational Hastings proof gate

**Status:** experimental, unmerged Python research prototype. This is an
algorithmically defined Metropolis Monte Carlo process on the simplified
2D HP model, **not molecular dynamics**, an atomistic folding predictor,
a measured transition rate, or a Lean-refined implementation.

Stacks on draft PR #94 (dual kinetic flux certificates), #93
(exact pivot kinetics and stochastic JEPA teacher), and #92
(proof-gated JEPA search). Main and the 61-PR integration are untouched.

## Why this solves a real computational bottleneck

The exact optimal transition sampler in PR #94 has to enumerate all
conformations and calculate an exponential-size event-cut graph before
it can propose a new state. `jepa_local_hastings.py` instead computes
the **exact transition row for just the current conformation** by
enumerating available pivot actions, a polynomial-per-step operation.

A learned network or any other untrusted policy supplies a fixed
state-dependent pivot-class distribution `Q_x`, with positive,
normalized, rational weights. Geometrically invalid proposals count
toward staying at the current state; all forward and reverse proposals
are independently reconstructed.

For a pair of *different* valid full conformations `x,y`:

```text
w(x) = q ** contacts(x)
Q(x,y) = sum of exact class-weighted pivot proposals x -> y
alpha(x,y) = min(1, w(y)*Q(y,x) / (w(x)*Q(x,y)))
P(x,y) = Q(x,y) * alpha(x,y)
```

Because the policy may depend on `x`, the simpler Metropolis ratio
`w(y)/w(x)` in PR #93 **is generally insufficient**. Hastings' reverse
proposal factor must be included.

The local verifier checks:

```text
w(x)*P(x,y)
  = min(w(x)*Q(x,y), w(y)*Q(y,x))
  = w(y)*P(y,x)
```

using exact rational operations. Thus, **for a fixed deterministic
policy**, every transition satisfies the detailed-balance identity;
row normalization follows by retaining the rejected mass as a
self-loop. No full conformation enumeration is needed for this
*local* calculation.

The generic finite-state stationarity and reversible-mixture
identities are proved separately in Lean
(`FoldBounds/JepaKineticStationarity.lean`, in PR #94).
That theorem does not formally verify this Python policy parser or
its pivot geometry; the implementation remains software-checked.

## Frozen-policy requirement

The policy **must be deterministic and frozen throughout the Markov
run**. Updating neural parameters midway creates a sequence of
different transition kernels; the usual time-homogeneous stationary
distribution and convergence argument no longer applies without
additional adaptive-MCMC assumptions. This prototype does not claim
such an adaptive-MCMC convergence theorem.

`StabilityGuard` rejects observed changes to a proposal distribution
for previously encountered conformations. It does **not** prove a
network is stable on unseen states. The caller must independently freeze
weights, action ordering, feature computation, inference mode, and
preprocessing and preserve their checksum/provenance in serious tests.

`quantize_ranks` converts finite neural scores into strictly positive,
normalized exact `Fraction` weights. It is robust to scale and avoids
silently treating inexact floating-point softmaxes as exact certificates.
Invalid scores, inconsistent weights, bad paths and nonsymmetric edge
receipts fail closed. The exact seeded sampling function draws integer
random variates proportional to rational transition probabilities.

## Independent local receipts and falsification

`EdgeReceipt` contains the source and destination paths, exact
forward/reverse proposal probabilities, exact Boltzmann weights, both
corrected transition probabilities, and their common stationary
flow mass. `verify_edge_receipt` independently regenerates all those
values using the same **frozen policy** and refuses tampering.

`audit_finite_mh` additionally performs a *separate complete-universe
small-system test*: enumerates all legal conformations, rebuilds every
local row, checks all row sums and pairwise detailed balance, and
compares the exact partition function to the earlier HP oracle.

This is important: local arithmetic can check proposed edges without
enumeration, but the guarantee that the **entire deployed model** has
the intended stationary measure depends on the frozen policy and
correct geometric implementation for all states. Unvisited states
cannot be exhaustively software-audited for long polymers.

## Reproduce

```bash
python -m pip install -e '.[test]'
python -m pytest -q tests/test_jepa_local_hastings.py
python experiments/benchmark_jepa_local_mh.py --steps 8
python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
python experiments/benchmark_jepa_local_mh.py --steps 8 --train
```

The baseline equivalence tests reconstruct **exactly** the same
transition matrices as the classical uniform-pivot Metropolis
reference in #93. Separate tests verify nonuniform, conformation-
dependent rational policies at `q=5/3`, hostile receipts, NaN
score rejection, model drift, deterministic integer sampling, and a
small frozen JEPA actor over a completely enumerated held-out
eight-residue universe.

The default large-chain case uses a **32-residue HP sequence** with
locally checked pivot transitions and **does not enumerate the
conformation state space**. Its short sampled trajectory proves
only that the local implementation is runnable and finite; not
that the neural policy accelerates mixing or improves protein physics.

## Scientific significance and next gates

This is a viable integration of model-guided proposals with
Boltzmann-preserving **local accept/reject mathematics**. Metropolis-
Hastings and detailed balance are standard prior art; no new theorem
or biological breakthrough is claimed from using them.

Before a major scientific claim: freeze and compare independent
whole-sequence training/test splits, then show statistically reliable
improvements in integrated autocorrelation, effective sample size per
CPU/GPU second, rare-event recovery, and variance for several
observables. Compare against optimized classical proposals,
nonreversible samplers and active-learning methods with the
model inference time included. Finally, establish physically meaningful
dynamics and independent structural/binding validation in real
molecular simulation and experiments, neither of which has been done.

Key mathematical prior art: Metropolis–Hastings (1953/1970) and
Madras & Sokal's self-avoiding-walk pivot sampler (1988),
[DOI:10.1007/BF01022990](https://doi.org/10.1007/BF01022990).
