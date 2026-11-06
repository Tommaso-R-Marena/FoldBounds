# Sharp proof-carrying thermal response and compressed spectral witnesses

**Scope:** additive Python research extension stacked on PR #46; no change to
the existing square-lattice checker, Lean sources, workflows, or other agents'
branches. These are exact-rational, **model-relative** certificates, not
Lean-verified Python or experimentally calibrated protein predictions.

## The exact rational variance relaxation

After validating the existing HP frontier, let known conformations have
scores `k_j` and exact weighted masses `w_j > 0`. Each unexpanded
subtree has a known contact-score enclosure `[a_i,b_i]` and a rational
mass cap `U_i`. *Relax* the unknown ensemble: allow any nonnegative
rational mass assigned to any integer score in the enclosure, with
mass in each subtree summing to at most `U_i`. This contains the actual
unknown HP Boltzmann score distribution.

For rational center `c` and candidate variance upper bound `V`,
the following **dual inequality** proves variance <= V for **every**
admissible score distribution:

```text
sum_j w_j * ((k_j-c)^2 - V)
  + sum_i U_i * max(0, max((a_i-c)^2, (b_i-c)^2) - V) <= 0.
```

Indeed, every possible score contributes at most its worst endpoint
squared deviation; the displayed inequality implies
`E[(K-c)^2] <= V`, hence `Var(K) <= V`.

To prove *sharpness for the relaxation*, the returned certificate also
contains a **feasible rational primal mass assignment** that attains
`E[K] = c` and `Var(K) = V`. The independent verifier recomputes
the inequality, validates score support, checks all mass caps, and
recomputes both moment equalities. This verifies exact optimality without
trusting the search procedure. It does **not** claim that the attaining
relaxed distribution corresponds to real self-avoiding conformations.

The optimizer uses exact fractional programming to choose cap/no-cap
mass assignments at the score interval endpoints for each trial center.
An adaptive envelope of the resulting rational quadratics is optimized
at its means, support endpoints, and pairwise crossings. Once the full
oracle agrees with the restricted envelope, at most two active models
are mixed to create the primal attainment. A resource budget fails closed.

The result is the sharp upper variance **for the stated continuous mass
relaxation**. If all mass is unknown, it reduces to the classical range
variance bound. It can be far from the true physical variance when the
subtree mass caps are loose.

## Exponentially large, independently counted score witnesses

The optional `use_monotone_spectra=True` derives exact HP contact-score
histograms for the union of four axis-monotone continuation languages
(E/N, E/S, W/N, W/S), using compressed dynamic programming.

Each cone has strictly increasing signed coordinate potential. Hence:
- No suffix self-intersection can occur.
- Nonconsecutive suffix sites cannot be lattice neighbors.
- The exact contact score consists of prefix contacts plus explicitly
  counted prefix-to-suffix contacts.
- The DP needs states `(endpoint, cross-contact score)`, not an
  explicit list of exponential many walks.

Each straight cardinal ray belongs to two cones. Subtracting the four
ray histograms from the sum of four cone histograms therefore counts each
valid family member **once**. Different frontier prefixes define
disjoint subtrees. Any separately supplied full-walk witness lying in a
cone is removed once by score before cone masses and moments are credited,
preventing overlap/double counting.

These genuine, nonenumeratively counted paths contribute known weighted
first/second moments and event/complement moments. Their total mass is
subtracted from the remaining subtree upper cap, allowing strictly
positive variance and event-sensitivity lower bounds even at length 30.

**Practical limitation:** for length 30 at q=1 and frontier depth 2, the
new lower bound is strictly positive but extremely small relative to the
loose upper bound. Counting exponentially many actual states does not,
by itself, solve the general exponential SAW problem or establish
predictive biophysics.

## Reproduce

```bash
python -m pip install -e '.[test]'
python -m pytest tests/test_variance_dual.py tests/test_monotone_spectral.py tests/test_thermal_response.py -q
python experiments/benchmark_dual_thermal.py
```

Independent checks include 160 randomized rational cases compared with an
exhaustive extreme-model and pairwise-mixture oracle; counterfeited dual and
primal certificates; independently enumerated 2D SAW suffix membership and
cone score histograms through length 8; all binary HP sequences through
length 7 against true Gibbs moments; and overlapping/invalid witness tests.

## Novelty and trust boundary

Variance/convex duality, Bhatia–Davis bounds, and monotone walk languages are
well-established. The implementation combines these tools with **externally
checkable primal and dual rational certificates** and disjoint,
nonenumerative *spectral* subset credit in FoldBounds' existing frontier
semantics. Do not claim discovery of a new inequality or publication
priority. See Bhatia & Davis (2000), Lim & McCann (2022), and the existing
`docs/MONOTONE_FAMILY_LOWER_BOUNDS.md`.

The base Python frontier checker and its interpretation remain outside
Lean's trusted kernel; no new Lean refinement has been proved. This work
does not validate 3D HP, empirical NMR, real temperature calibration,
protein kinetics, or experimentally observed conformational switching.


## Additional refinement: submultiplicative exact rooted-SAW block caps

The optional `use_saw_count_cap=True` tightens the generic unexpanded
subtree bound `3**r * q**kmax` without enumerating full-length walks.
For `k <= block_max` (default 8), count **all** rooted anchored 2D
self-avoiding extensions of `k` sites by short direct enumeration,
obtaining the integer `A(k)`. Each continuation of an arbitrary valid
prefix with `r = k1 + ... + kj` steps injects into its sequence of
rigidly oriented, independently legal rooted step blocks; forgetting
collisions between blocks only enlarges the set. Thus

```text
number_of_full_extensions <= A(k1) * ... * A(kj).
```

An integer dynamic program picks the smallest block-count product over all
partitions of `r` with allowed block sizes. The certificate stores the
partition and exact product. The validator re-enumerates all short base
blocks and checks that the partition sums to `r`, the product is
exact, and it does not exceed `3**r`. Multiplying by the existing
parity-aware `q**kmax` cap then gives a strictly stronger **model-relative**
upper bound in many cases. This is a standard submultiplicativity argument,
not a new asymptotic growth constant or proof that the full SAW model can
be solved efficiently.

In local all-H, `q=1`, depth-two tests, enabling count caps in addition
to cone spectra reduced the original *partition* upper bound by factors
1.61 (n=10), 4.80 (n=30), 311 (n=100), and 130,848 (n=200), relative
to cone spectra with the original `3**r` cap. At n=100 the positive
variance lower bound increased by ~96,539x. **Absolute lower bounds
remain extremely small.** No experimental protein conclusion follows.

Run `python -m pytest tests/test_saw_count_cap.py` to reproduce the
independent rooted-walk oracle, forged-product rejection, randomized
exact-Gibbs enclosure tests, and long-chain smoke tests.
