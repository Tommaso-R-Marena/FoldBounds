# Exact kinetic cut certificates: why equilibrium preference does not imply accessibility

**Scope:** an isolated research prototype for the finite, rooted, 2D square-lattice
HP toy model. These algorithms **do not establish real-protein folding rates,
novel physical laws, experimental predictive accuracy, or Lean verification**.
They complement equilibrium certificate work without editing it.

## Statement and mathematically specified dynamics

A conformation is a self-avoiding walk of n square-lattice sites beginning
`((0,0),(1,0))`; no reflection quotient. Contact score `k(x)` counts
nonconsecutive adjacent H-H pairs. For exact rational `q>=1`, the Gibbs
stationary weight is `pi(x)=q^k(x)/Z` on the *entire* rooted state space.

At each discrete update, choose uniformly one of the `n-2` movable bead
indices 2..n-1 and one of the four neighboring positions of its immediately
preceding bead. Reject if the proposed position overlaps a fixed bead,
breaks an adjacent chain bond, or equals the original position; otherwise
accept with `min(1,q^(k(y)-k(x)))`. Invalid and rejected proposals are
self-loops. Thus every legal distinct edge has symmetric proposal probability
`1/M`, with `M=4(n-2)`, and exact detailed balance

    pi(x) P(x,y) = q^min(k(x),k(y)) / (M Z).

**The move graph need not be irreducible.** The global Gibbs distribution is
stationary, not necessarily the unique limiting distribution; this is a
central limitation rather than an omitted assumption.

## Two independent finite-time no-go guarantees

Choose a set `S` containing an explicitly identified start `x` and disjoint
from the target set `A={y : k(y)=max k}`. Define its directed boundary flow

    F_q(S) = sum_{u in S, v outside S, u~v} q^min(k(u),k(v)).

Each adjacent pair is counted once in the sum oriented *out* of S. Let
`Z_S=sum_{u in S}q^k(u)` and `T` be a nonnegative integer number of updates.
Then both inequalities hold **without requiring S to be metastable,
reversible mixing within S, or an irreducible global state graph**:

    P_x(hit A by T) <= P_x(exit S by T)
                     <= min(1, T F_q(S)/(M q^k(x))).       [fixed start]

    P_{pi(.|S)}(hit A by T)
                     <= min(1, T F_q(S)/(M Z_S)).         [basin-Gibbs start]

Proof: `delta_x <= pi/pi(x)` componentwise and `pi(.|S)<=pi/pi(S)`.
Multiplication by any Markov transition kernel preserves this pointwise
inequality, while `pi P^t=pi`. Thus at every update, the probability of
crossing the cut is at most the corresponding equilibrium edge flow divided
by `pi(x)` or `pi(S)`; union-bound all T steps. This remains valid even when
paths leave and reenter S, and hence bounds **first** exit. If `F=0`, S is
closed and neither start can reach any target outside it, ever.

Let `h=max_{S}k`, `a=# {u in S : k(u)=h}`, `B=# boundary edges`,
`c=max_boundary min(k(u),k(v))`. For all real `q>=1` (not just sampled q),

    F_q(S)/(M Z_S) <= (B/(M*a)) q^(c-h),
    F_q(S)/(M q^k(x)) <= (B/M) q^(c-k(x)).

These are global in q, elementary monomial dominance certificates. In the
closed case both upper bounds are identically zero.

## Proof-carrying exact minimum-cut improvement

The untrusted producer first proposes a score-superlevel component, then uses
an exact-integer augmenting-path algorithm to find a *different*, potentially
stronger boundary `S` minimizing `F_q(S)` at a declared rational q among all
cuts separating the **same fixed starting state** from *every* maximum-score
state. It clears rational denominators by `denominator(q)^max_score`.

The independently written verifier does **not** call that optimizer or trust
its ranking. Instead it independently enumerates all canonical SAWs and legal
moves, validates S and the target exclusion, recomputes all scores and exact
cut capacity, then checks the supplied sparse integral **max-flow dual**:
edge existence, capacity constraints, antisymmetric net edge flow, flow
conservation at all intermediate states, net inflow into target sinks, and
exact equality of feasible flow value and cut capacity. By max-flow/min-cut
duality, this proves the proposed cut is **globally minimum for that fixed
start, target class, move graph, and q**. The same geometric cut gives sound
(non-optimal) bounds at other q, but its flow proof is only accepted at the
q recorded in the certificate.

This verifier is tested Python, **not Lean-kernel-checked**. Both finite-state
enumerators run exponentially and enforce a strict all-state budget; no
partial enumeration is ever silently accepted.

## Exact reproducible examples (`q=10000`, `T=10^9` proposals)

| HP sequence | All rooted walks | Equilibrium maximum-contact probability | Cut | Finite-time no-go |
|---|---:|---:|---:|---|
| H x 9 | 1,479 | 0.998681473 | 5-state basin, 2 edges, barrier exponent 3 | fixed-start hit <= 1/14,000 |
| H x 11 | 11,025 | 0.999379432 | 445-state optimal cut, 15 edges, barrier exponent 4 | fixed-start hit <= 1/24,000,000 |
| H x 12 | 30,073 | 0.998569532 | one *completely frozen* nonnative conformation | reaching a maximum-contact fold is impossible for every T |

For `n=11`, the originally constructed 321-state score-superlevel basin had
195 cut edges of contact score zero. The exact mincut improves this to 15
cut edges **13x** at the same q and fixed seed. The independently checked
conditional-basin bound for the optimized cut is even tighter than the
fixed-start bound, but makes a specific initial-distribution assumption.
This is not an experimental folding-time measurement: one step is exactly one
local bead/compass proposal. The 12-mer is frozen only under this restricted
local move set; an explicit legal nonlocal pivot invalidates a general
"permanently trapped polymer" claim.

## Reproduce and try to falsify

    python -m pytest tests/test_kinetic_barrier_certificates.py -q
    python experiments/kinetic_cut_benchmark.py

The 25 tests include independently re-enumerated counts n=4..12, two
independent local move definitions and inverse/detailed-balance checks,
exact one- and three-step dynamics versus the union bound, q=1 and rational
q=3/2 max-flow/min-cut duality, nonnative closed-class detection, every
short H/P composition, adversarial member deletion/addition, spurious paths,
corrupted flow balances/capacities, q substitution, and fail-closed budgets.

## Prior art and scientific interpretation

Kinetic traps, local moves, high-barrier folding pathways, and the failure of
simple Metropolis dynamics to reproduce equilibrium on reasonable simulation
times are **established** research topics. See H. S. Chan & K. A. Dill,
*J. Chem. Phys.* **100**, 9238-9257 (1994), DOI 10.1063/1.466677, on exact
2D HP dynamics and move-set dependence; O. Collet, *Phys. Rev. E* **67**,
061912 (2003), DOI 10.1103/PhysRevE.67.061912, on kinetic rigidity and
correct dynamics. Maximum-flow/minimum-cut and Markov-chain inequalities
are established tools as well. **No priority or new physics is claimed.**

The useful FoldBounds extension is a completely explicit, independently
replayable, adversarially checked certificate that equilibrium confidence
can be misleading about kinetics. A future research step is to use stronger
move sets, measure actual protein data, and formalize the general Markov
cut lemma in Lean; none has been accomplished in this branch.
