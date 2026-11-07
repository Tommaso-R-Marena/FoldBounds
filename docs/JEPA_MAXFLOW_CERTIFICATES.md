# Proof-carrying kinetic event flux via exact min-cut

**Status:** research prototype, Python exact-rational checker (not Lean
refined); benchmark concerns a finite 2D HP model and **algorithmic Monte
Carlo steps**, not real protein folding rates.

This module stacks on the reversible pivot-Markov oracle in PR #93 and the
proof-gated JEPA modules in PR #92. It does not modify them.

## Precise finite theorem being instantiated

Fix an explicitly enumerated finite conformation universe `S`, exact
Boltzmann weights `w(x)=q**contacts(x)`, normalization `Z=sum_x w(x)`,
a contact event `A subset S`, and an undirected move graph `G`. Edges
of `G` are conformations connected by a valid anchored pivot move.

Consider **any reversible row-stochastic** transition kernel `P` with
stationary probability `pi(x)=w(x)/Z` and all non-self transitions on
`G`. Reversibility implies symmetric nonnegative stationary edge
flows:

```text
f(x,y) = pi(x) * P(x,y) = pi(y) * P(y,x)
sum_{y != x} f(x,y) <= pi(x)
```

The one-step equilibrium flux out of event `A` is

```text
Phi(P,A) = sum_{x in A, y outside A} f(x,y).
```

Other edges cannot increase this specific objective, and every
nonnegative collection of **cross-event** flows satisfying the vertex
stationary-mass budgets gives a reversible kernel: set `P(x,y)=f(x,y)/pi(x)`
on these edges, and retain the remaining probability as a self-loop.

Consequently **the maximal possible event-crossing flux among such
reversible kernels** is exactly a capacitated bipartite max-flow:

- source `s -> x in A`, capacity `w(x)`;
- each allowed pivot edge `x in A -> y not in A`, capacity `Z+1`;
- each `y not in A -> t`, capacity `w(y)`.

The max-flow value divided by `Z` is the exact optimal stationary
event-crossing probability per step. The `Z+1` capacity ensures any
finite optimal cut can be chosen without cutting a cross edge. All
arithmetic uses `Fraction` and therefore no floating tolerances.

The theorem uses classical maximum-flow/minimum-cut duality and elementary
reversible-chain identities: **these are established mathematics**, not
novel proofs of a new physical theory.

## Why the certificate is independently checkable

The untrusted producer `propose_maxflow_certificate` runs a deterministic
exact Edmonds–Karp algorithm. Its serialized mathematical witness consists
of:

1. Positive rational flow amounts on named pivot edges crossing the
   event/complement partition;
2. A cut's source-side conformation indices;
3. Exact rational total flow and proposed cut capacity.

The verifier `verify_maxflow_certificate` separately reconstructs the
complete finite graph from the checked underlying pivot Markov kernel,
checks every used edge and positivity, every vertex capacity, uniqueness,
the candidate cut membership, and the exact equality

```text
sum(accepted_cross_flow) == claimed_flow == claimed_cut == recomputed_cut.
```

A feasible flow is a lower bound on the maximum, and *any* cut is an upper
bound. Their equality proves optimality of the finite math problem.
The verifier does **not trust the producer's search, residual BFS, or
numeric result**. It is still a software-checked mathematical guarantee,
not end-to-end Lean verification.

The witness `construct_reversible_kernel` blends the exact optimal
cut-flow kernel with the original uniform-pivot Metropolis kernel:

```text
P_blend = (1-alpha) * P_uniform + alpha * P_maxflow
0 < alpha < 1.
```

Both terms have the **same** exact Boltzmann stationary distribution
and are separately reversible. Their mixture is therefore exactly
reversible. All positive baseline edges remain positive for `alpha<1`,
preserving its communication graph support. The independent
`verify_reversible_kernel` reconstructs and checks the entire mixture.
Its event flux equals the same convex combination of the baseline and
max-flow flux, so cannot be smaller than baseline for the chosen event.

**Crucial:** This does not claim that the original pivot chain is
irreducible or that improved stationary cut flux improves global
mixing, asymptotic variance, MFPT, folding time, biological
accessibility or any measured protein rate. The optimized transition
kernel is global, state-dependent, and computationally expensive to
construct. It is not the simple state-independent proposal family of
PR #93.

## Reproducibility

```bash
python -m pip install -e '.[test]'
python -m pytest -q tests/test_jepa_flux_maxflow.py
python experiments/benchmark_jepa_maxflow.py
```

The benchmark runs finite root-anchored HP chains, `q=2`,
with a defined contact-event threshold. Every reported primal optimum is
also checked as a dual minimum cut; the 50/50 mixed witness is
independently rechecked for stochasticity, detailed balance and flux.
Input sizes are small (typically 71, 543 or 1,479 conformations).

In local exhaustive checks at `q=2`, threshold 1:

| Sequence | Conformations | Maxflux / uniform | 50% blend / uniform |
| --- | ---: | ---: | ---: |
| HHHHHH | 71 | 5.00 | 3.00 |
| HPHHPPHH | 543 | ~8.02 | ~4.51 |
| HHPPHHHH | 543 | ~7.96 | ~4.48 |
| HPHPHPHH | 543 | ~8.67 | ~4.84 |

These are **exact finite-model comparisons**, not experimental speedups,
and were known/inspected during development. They are not prospective
JEPA generalization evidence.

## Scientific next steps

1. Extend the formal assurance: encode a generic rational
   feasible-flow/min-cut dual-certificate theorem in Lean, and prove
   the finite kinetic decoder/refinement relation. That is not yet done.
2. Develop a learnable local, state-dependent **proposal producer** that
   approximates the certified max-flow witness without exhaustive state
   enumeration. A JEPA actor is one candidate; its outputs must still
   pass a separate reversibility/checker gate.
3. Quantify proposal construction time and sample-time mixing, autocorrelation,
   asymptotic variance and independent-event transfer; optimization on one
   event is likely to sacrifice other kinetic observables. Study the
   tradeoff as a multiobjective optimization, not merely cut flux.
4. Compare with classical pivot and adaptive/nonreversible MCMC baselines
   from the literature on matched ensembles.
5. Transfer only scientifically justified techniques to realistic
   coarse-grained / atomistic models, and validate using independent
   simulation trajectories and experimental data. **Do not identify
   Monte Carlo steps with physical time.**

References: Madras and Sokal (1988), *The pivot algorithm: a highly
efficient Monte Carlo method for the self-avoiding walk*,
[DOI:10.1007/BF01022990](https://doi.org/10.1007/BF01022990);
classical max-flow/min-cut duality and Metropolis–Hastings.


## Exact additional guarantee: lower stationary lag-one event correlation

Let `I_t=1_{X_t in A}` under a chain already in the common stationary
Boltzmann distribution, with `p=Pr(I_t=1)`. The joint probability that
both adjacent steps are in the event equals `p-Phi(A)`. Thus

```text
Cov(I_t, I_{t+1}) = p*(1-p) - Phi(A)
rho_1 = 1 - Phi(A)/(p*(1-p)),  for 0 < p < 1
```

A higher *checked* stationary event flux therefore gives a *lower*
**lag-one autocorrelation of that same binary event** at stationarity.
`stationary_event_lag1` and `stationary_event_lag1_reversible`
calculate this exactly from rational masses. Independent tests recompute
the correlation from `Pr(I_t=1,I_{t+1}=1)-p²` using the entire verified
transition matrix rather than reusing the flux identity.

This is a genuine additional model-relative performance identity,
but does **not** imply faster asymptotic decorrelation, lower integrated
autocorrelation, lower long-run statistical error, improved global
mixing, or faster experimental folding. The event's correlation can
become negative under an aggressive alternating policy, which is
not evidence of physically realistic kinetics.
