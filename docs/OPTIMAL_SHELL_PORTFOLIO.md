# Optimal certificate-budgeted antichains for FoldBounds shell lower masses

**Status:** Exact rational Python research implementation, tested against an
independent exact geometry oracle and brute-force enumeration of every feasible
small prefix-tree portfolio. **Not Lean verified**, and not a claim of
algorithmic novelty or real-protein predictive utility. Stacked on shell
construction of PR #28, independent of Aristotle's formal PR #26.

## Mathematical construction

Let p be an anchored self-avoiding prefix, and let F(p) denote the finite
set of its valid completed extensions. For every p, the extremal-quadrant
algorithm of PR #28 defines a computably scored subset S(p) ⊆ F(p), with
exact rational event and complementary masses from a finite formula. For
q≥1 its mass is nonnegative. More generally positivity of rational q>0
is sufficient for the subset lower bound.

A valid **portfolio** is the leaf set L of a finite complete rooted prefix
tree. A node may either STOP, obtaining the restricted family S(p), or
EXPAND all its geometrically legal one-step extensions. No node can be both
an accepted leaf and have descendants, and every live child must be covered.

The full completion sets F(p) of distinct prefix-free leaves are disjoint:
each completed SAW has only one prefix of any given length, and comparable
prefixes cannot both be leaves. Thus the union of selected families
⋃_{p∈L} S(p) is disjoint. Its mass is **exactly** the sum of independently
re-scored shell masses at its leaves. This is a lower mass for the full model,
not an assertion of exhaustive coverage by the selected shell families.

## Exact tree knapsack

For a fixed maximum core depth D and integer certificate leaf budget B,
let DP_p[b] be the maximum restricted event mass available from any
complete prefix-tree antichain rooted at p using **exactly b leaves**.
For complement optimization substitute the complement mass. The recurrence
is

    DP_p[1] >= shell_event_mass(p)               (STOP)
    DP_p[b] >= max_{b_1+...+b_s=b} Σ_i DP_child_i[b_i]  (EXPAND)

for all possible allocations among all s legal children. The code computes
an exact maximum in every reachable budget cell (not merely a bound), stores
the corresponding leaf paths, and finally selects the highest mass for any
feasible b≤B. The proof is structural induction: every legal antichain either
stops at p or expands all its children, and conversely every selected child
portfolio union is legal. Base case is D or dead ends. All DP coefficients
are Fraction, with no floating-point objective approximations.

The optimizer's portfolio is **untrusted**. The separate
`verify_shell_portfolio` validates source parameters and geometry, rejects
duplicates, shadowing, omitted branches, malformed leaves and unreachable
nodes, then recomputes all shell masses from source conformation geometry.
No proposer-supplied totals or contact scores are accepted as evidence.

Because event and complement are disjoint event classes, independently
optimized event and complement portfolios can provide separate lower bounds;
**within each class** overlapping portfolios must be combined using max,
not sum. Their event and complement lower masses can be added to lower-bound
Z since the target events are disjoint. Neither portfolio is a Lean proof.

## Complexity

The number of candidate prefixes to depth D grows exponentially with D.
Per candidate, the DP considers up to B leaf allocations per child, with
convolutions costing O(B²) in the simple implementation. The optimizer can
therefore be more expensive than constructing a uniform-depth shell
certificate, even if its resulting leaf certificate is much smaller.
The verifier of a supplied portfolio recomputes exact local shells and
traverses its finite prefix tree; for fixed depth and bounded coordinates,
its work is polynomial in the listed certificate and geometric dimensions.
No polynomial-time claim is made for unconstrained n or portfolio search.

The optimizer aims to maximize **restricted event mass** or **restricted
complement mass**, not necessarily to minimize the final probability interval,
since independent upper masses and cross-class coupling matter.

## Matched certificate-budget evidence

Reproduce with `python experiments/compare_shell_portfolio.py` and
`python -m pytest tests/test_shell_portfolio.py`.

For mixed HP length 12, q=3/2 and event threshold 2:

| Method | Core depth (fixed or maximum) | Leaf count | Event lower |
|---|---:|---:|---:|
| Uniform shells | 5 | 25 | 576 |
| Optimal variable-depth portfolio, budget 25 | 9 | 25 | 4455/4 (1113.75) |
| Uniform shells | 6 | 71 | 576 |
| Optimal variable-depth portfolio, budget 50 | 9 | 50 | 7587/4 (1896.75) |

For all-H length 12 at q=2, threshold 2, uniform depth 5 has event
lower 1,536 from 25 leaves, while optimized budget 25 depth 9 has
2,848; for all-H length 16 the values are 24,576 versus 46,048.

The comparisons match certificate **leaf counts**, not proposer wall-clock
work. Optimal planning to maximum depth 9 explored 2,326 candidate prefixes;
fixed depth 5 examined 25. Any performance claim must disclose that extra
search cost. The independent tests include brute-force optimality for small
prefix trees and full-completion membership/oracle contact scoring for n≤7.

## Trust, interoperability and limitations

A complete validated portfolio is *necessary for these exact selected-family
sums* but not necessary for arbitrary valid lower bounds. Rejecting an
incomplete portfolio is conservative. No generated tuple is assumed sound
merely because the optimizer returned it.

No theorem yet proves equivalence of this Python checker to Aristotle's Lean
frontier checker or correctness of the parser. Future work: Lean 4 theorem
for shell-family cardinal-ray inclusion/exclusion and sound antichain
aggregation, then certificate parsing and byte/parameter correspondence.