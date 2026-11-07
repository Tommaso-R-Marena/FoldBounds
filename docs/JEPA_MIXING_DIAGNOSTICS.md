# Kinetic mixing falsification: stationary lag-one is not enough

**Status:** exploratory finite 2D HP Markov-chain comparisons, not
physical protein folding. This work is additive to the exact-rational
flow certificates from PR #94 and the local Hastings interface in PR #95.

The max-flow optimizer in PR #94 has a rigorously checked improvement for
the selected *one-step* equilibrium event flux. That alone cannot prove
faster convergence or lower statistical error for a stationary time series.

## New numerical diagnostic

For each fully enumerated small HP world, the stationary distribution
and complete transition matrices are first checked with exact Fractions.
For a centered observable `f(x)` under an ergodic finite chain, the
integrated autocorrelation time is

```text
tau_int = 1 + 2*sum_{t>=1} Corr(f(X_0), f(X_t))
        = 2*<f,(I-P)^(-1)f>_pi/Var_pi(f) - 1.
```

The second formula is evaluated by a constrained sparse Poisson solve.
An **independent finite correlation series** is also summed to numerical
convergence and compared with the Poisson result. A large solve residual
or a mismatch causes an explicit error/abstention. These are numerical
diagnostics with finite precision, not Lean proofs or rigorous intervals.
An irreducibility/appropriate convergence assumption is needed; the
reversible mixture retains the baseline's positive edge support, but
this test is not a general theorem establishing baseline irreducibility.

Run:

```bash
python -m pip install scipy numpy
python experiments/benchmark_jepa_mixing.py
```

## Local experimental comparisons

All rows use `q=2`, contact-event threshold `>=1`, the original
uniform-pivot reference versus a 50/50 mixture with the max-flow
cross-event reversible kernel.

| Sequence | Event tau_int uniform | Event tau_int mixed | Ratio (higher=better) |
| --- | ---: | ---: | ---: |
| HHHHHH | 6.664 | 2.926 | 2.28x |
| HPHHPPHH | 12.430 | 5.340 | 2.33x |
| HHPPHHHH | 12.049 | 4.495 | 2.68x |
| HPHPHPHH | 12.377 | 3.139 | 3.94x |
| HPHHPPHHH | 14.117 | 5.615 | 2.51x |

**Negative controls matter:** for the *same* HPHHPPHH model, the
radius-of-gyration-squared observable worsens from an integrated
autocorrelation time of approximately 9.865 to 10.808. In HPHPHPHH it
worsens from approximately 8.594 to 10.248. The contact-count observable
can improve by much less than the targeted event. Therefore, neither a
higher event flux nor a lower event lag-one autocorrelation implies a
universal mixing-time improvement.

The values above were inspected during development and are NOT
prospectively held-out discoveries. Each diagnostic must be replicated
with frozen source revisions, multiple random seeds for sampled checks,
independent implementations, and larger sequences before publication.

## Appropriate next target

Optimize a multiobservable, resource-normalized objective that penalizes
degrading important independent observables, while retaining strict
reversibility/thermodynamic checks. A JEPA actor should be compared with
classical adaptive and nonreversible Monte Carlo, taking model inference
cost, training-oracle enumeration cost and correlation uncertainty into
account. Claim improved **scientific sampling efficiency** only after
prospective held-out sequences and realistic molecular models support it.

This file makes no claim about physical folding times, experimental
barriers, structural accuracy, or protein therapeutics.
