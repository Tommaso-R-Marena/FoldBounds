# JEPA-guided reversible lattice kinetics: exact-oracle pilot

**Research status (7 November 2026):** experimental Python exact-rational
checker, not Lean-refined; *no real-protein folding-rate claim* and *no
mathematical priority claim*. This branch stacks on draft PR #92 and is
independent of main and the larger draft integration PR #91.

## Why add kinetics to FoldBounds?

The existing static partition function does not determine kinetic rates:
a **transition mechanism must be specified separately**. We choose a
canonical, explicit, discrete-time **Metropolis pivot Monte Carlo** process
over all rooted self-avoiding walks for short sequences. A transition means
one algorithmic pivot proposal, not one unit of physical time. This is an
independently checkable training oracle for dynamics-based JEPA.

### Kinetic model and precise guarantee

For a full conformation `x`, let `c(x)` be its HP H–H nonbonded
contact score, `q>=1` an exact rational Boltzmann factor, and

```text
w(x) = q**c(x)     pi(x) = w(x) / sum_y w(y)
```

The initial two lattice sites are anchored. Pick an internal pivot index
`k=1,...,n-2` and one of two move families: a quarter turn chooses
+90 or -90 degrees with equal probability; a half turn chooses 180 degrees.
Collisions are rejected. Every proposed nonselfintersecting move is accepted
with exact rational probability `min(1, w(y)/w(x))`. Each quarter-turn
proposal has its reverse within the same move family, at the same class
weight, and 180-degree turns are self-inverse.

For any **state-independent** strictly positive rational probability
distribution over the move families, proposals are symmetric. Therefore

```text
pi(x) P(x,y) = pi(y) P(y,x)
```

by the elementary Metropolis detailed-balance identity. The independent
verifier `verify_kernel` reconstructs the entire sparse matrix, checks
every rational row sums to one, and explicitly checks both directions of
every nonzero edge. This guarantees a stationary distribution for that
finite software model, but does not establish unique stationary convergence
if the restricted move graph were disconnected.

For a contact-event subset `A={x: c(x)>=threshold}`, the exact stationary
per-proposal flux from A to its complement is

```text
Phi(A) = sum_{x in A,y not in A} pi(x) P(x,y)
```

and the conditional one-step event-exit probability is `Phi(A)/pi(A)`.
Neither quantity is a physical folding rate, a mean first-passage time or
a guarantee of improved overall Markov mixing.

### Exact flux-maximizing proposal weights

Let `F_c` be the exact equilibrium event-crossing flux if the entire
proposal mass were allocated to pivot/rotation class `c`. Because each
move class is individually reversible and Metropolis accept/reject is
independent of the class probability, flux is the linear function

```text
Phi(p) = sum_c p_c F_c
```

Constrain `p_c >= epsilon > 0` and `sum_c p_c=1`, where `epsilon`
is at most the minimum weight in the uniform-pivot baseline. Then the
**global optimum within this restricted class of proposals** assigns
`epsilon` to every class, and the remaining probability to the class
with maximal `F_c`. This is elementary linear optimization over a
simplex, not a newly discovered theorem.

The uniform pivot proposal lies in the feasible set; consequently the
optimum's equilibrium event-crossing flux cannot be smaller than
uniform's. The code independently reconstructs the chosen kernel and
replays detailed balance and the exact flux identity before returning a
result. This requires *full enumeration* and thus does not scale to
biomolecular conformational spaces.

### JEPA student and proof-gated release

`jepa_kinetic_teacher.py` creates stochastic training examples with exact
conditional future-state probabilities, including Metropolis rejections.
JEPA's target is the probability-weighted **mean of EMA-encoder successor
embeddings**. In a stochastic system this is not generally the embedding
of the mean successor state. Auxiliary targets predict actual conditional
event exit and expected contact-score change.

The student ranks pivot classes on entirely different **held-out HP
sequences**, without inspecting their crossing-flux labels while ranking.
The independently checked **proposal gate** quantizes that ranking into a
strictly positive exact rational class distribution. If the rebuilt
kernel's exact flux is worse than uniform, the gate discards the proposed
policy and uses uniform. Invalid lengths, NaNs, infinities, and nonfinite
models are refused. The gate additionally checks the proposal is no
better than the exact finite optimum.

Important limitation: the gate *itself* enumerates all conformations,
rebuilds transition matrices, and runs the exact oracle. Any apparent
JEPA performance guarantee therefore comes from the expensive independent
checker, **not the neural model**, and is not yet efficient for long proteins.

## Reproduce

```bash
python -m pip install -e '.[test]'
python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
python -m pytest -q tests/test_jepa_reversible_kinetics.py tests/test_jepa_kinetic_teacher.py
python experiments/benchmark_jepa_flux.py
python experiments/benchmark_jepa_kinetic_teacher.py --seed 11 --epochs 12 --train-states 24
```

The exact benchmark records all rational metrics, canonical state counts,
chosen actions, and measured execution time. The student benchmark keeps
`HPHHPPHH`, `HHPPHHHH`, and `HPHPHPHH` completely out of its
training-sequence set. It reports rejected/accepted proposals, oracle
regret and the **full-enumeration gate cost**; no published scientific
baseline superiority is claimed. Use multiple seeds and larger frozen
heldout suites before drawing a predictive conclusion.

## Initial algebraic finite-model finding (not a folding breakthrough)

At `q=2` and contact threshold 1, exact enumeration of the
eight-residue sequence `HPHPHPHH` gives a substantial *algorithmic*
flux increase with the optimal symmetric pivot allocation relative to
uniform pivot proposals (approximately 2.7x in the first smoke run).
This is a complete small finite toy model, **not experimental protein
data**. A larger physical system has different energetics, constraints,
and dynamics. Optimizing one event flux may worsen other observables,
mean first-passage times, asymptotic variance or global mixing.

## Prior art and scientific boundaries

Classical pivot algorithms, reversible MCMC, Metropolis-Hastings and
optimized/adaptive proposals are established prior art. See Madras and
Sokal, *The Pivot Algorithm* (1988), J. Stat. Phys. 50:109–186,
DOI: [10.1007/BF01022990](https://doi.org/10.1007/BF01022990).
Exact proposal-weight optimization here is a finite-model baseline
and **not** asserted to be novel.

For a potential consequential contribution, show, prospectively, that
JEPA predicts favorable proposals on unseen, larger *unverifiable*
worlds with useful confidence guarantees and lower *total computational
work*, while maintaining independently checked thermodynamics. Compare
against fixed heuristic, exact optimal small-world oracle, random/uniform,
traditional adaptive MCMC, and nonreversible methods with matched
transition/proposal budgets. Separately evaluate physical dynamics using
real OpenMM/GROMACS trajectories and experimental observables. None of
those claims has been established by these toy-model pilots.
