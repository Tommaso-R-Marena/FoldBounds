# Thermodynamic response enclosures (Python research prototype)

This is a strictly additive module for the anchored two-dimensional HP model. It consumes
an existing reference Certificate, checks it with verify_certificate, and optionally accepts
a set of extra full-conformation witnesses. No Lean sources, core search implementation,
CI configuration, or existing solver interfaces are changed.

**Trust boundary:** Guarantees are conditional on the explicit finite 2D HP model and
the mathematical validity of the Python reference verifier. The Python implementation is
not proven equivalent to Lean, and these outputs are not experimental protein predictions
or physical heat-capacity measurements.

## Run

    python -m pytest tests/test_thermal_response.py -q
    python experiments/thermal_response_demo.py

Python 3.11+ with the existing package installed; pytest is only needed for tests.

## Precisely what is bounded

For q >= 1 rational and an integer HP contact count k(c), the finite Gibbs ensemble
has Z(q) = sum_c q**k(c). The new consumer returns exact-rational intervals for:

1. Mean contact score E_q[k].
2. Second moment E_q[k**2].
3. Variance Var_q(k) = d E_q[k] / d log(q).
4. Threshold-event susceptibility d P_q(k>=t)/d log(q) = Cov_q(1_{k>=t}, k).
5. A witness-improved partition lower bound and event probability interval.

The derivatives are finite-sum identities, not numerical finite differences.
Both susceptibilities are nonnegative because threshold membership is monotone in k.

## Proof sketch of the finite-model enclosure

The existing verifier checks certificate coverage, validity of every SAW prefix,
nonoverlapping frontier nodes, dead-end claims, and the exact rational parameter q.
For each terminal node the new module recomputes the full contact score, weight, and
first two score moments.

For every unexpanded prefix with r remaining residues, every completion has a
contact count in [a,b], with a already verified contact score and b obtained by
conservatively counting all parity-compatible future H-H pairs. A universal
mass cap is U = 3**r * q**b. The unknown actual subtree mass x obeys 0 <= x <= U.

For a known total mass W and known first moment M, the unknown-subtree mean
score lies in [a_i,b_i]. Lower/upper ensemble mean bounds follow from exactly
optimizing the relaxed fractional problems

    min (M + sum a_i*x_i) / (W + sum x_i)
    max (M + sum b_i*x_i) / (W + sum x_i)

over 0 <= x_i <= U_i. With W > 0 the exact minimum is attained by filling a
prefix of groups sorted by increasing a_i, and the maximum by filling a prefix
sorted by decreasing b_i. This follows from the sign of the partial derivative
of a weighted average with respect to each new mass. Check every sorted prefix;
no floating point or heuristic nonlinear optimizer is used. When W=0, take the
minimum/maximum possible group score, conservatively. Apply the same method
to the nonnegative squared score k**2.

For variance, the pair identity gives

    Z**2 * Var(k) = sum_{i<j} w_i*w_j*(k_i-k_j)**2.

The pairwise contributions from known states cannot be cancelled by unknown
states; dividing their sum by a verified upper bound on Z**2 provides a valid
lower bound. Another lower bound is max(0, M2_lower - mean_upper**2).
The upper bound is the smaller of M2_upper - mean_lower**2 and the
Bhatia-Davis inequality (b-mu)*(mu-a), maximized over the certified mean
and score-support intervals.

For a threshold event, every event score exceeds every complement score.
Known event/complement cross-pairs therefore provide a nonnegative lower
bound on Cov(1_event,k) even if further states remain unexpanded. The upper
bound uses p*(1-p)*(b-a) maximized over the certified event probability interval.

## Independently checked witnesses

Each optional witness must be a tuple of exactly n integer-coordinate lattice
sites, anchored by (0,0),(1,0), with adjacent steps and no self-intersection.
Duplicate witnesses, witnesses already represented by a terminal frontier
node, and witnesses not belonging to exactly one bounded subtree are rejected.
Their exact Boltzmann mass is credited to known moment sums and subtracted from
the corresponding upper mass cap. This prevents double counting and makes
the residual uncertainty smaller without modifying the original certificate.
If all residual caps reach zero, every response interval becomes an equality.

## Test coverage and limitations

The regression suite checks every H/P sequence of length 2 through 6,
using an independently coded self-avoiding-walk enumerator, exact fractional
weights, four event thresholds, partial and complete frontiers, and multiple
frontier depths. Further length 7-9 randomized and deterministic cases test
the interval enclosures. Fractional optimization is separately compared
against all corner assignments for small artificial instances. Adversarial
cases test invalid frontiers and witness rejection; the last tests verify
strictly improved lower response bounds and complete witness exhaustion.

Shallow frontiers without verified exact states may still return broad
intervals. These are safe model-relative enclosures, not automatically
informative physical predictions. Incorporation of stronger locally verified
subtree caps and a Lean refinement theorem remains future research. No
claim of scientific novelty or real-protein validation is made.
