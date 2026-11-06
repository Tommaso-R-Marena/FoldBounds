# Exact global design of contact-event probabilities, for all q>1

**Additive research contribution:** This draft stacks on #54, #52, #48, #47 and #44. It creates only new files. It leaves existing Lean sources, physical-motion research, RNase A blind tests, main model verifier, CI configuration, and every other agent's branch untouched.

## Why a new objective matters

A large partition sum `Z(q)` measures total Boltzmann weight and is **not**, by itself, the probability of a specific conformational feature. A more targeted finite-model observable is the event

    A_r = {complete conformations with at least r nonconsecutive H-H contacts}.

Define its exact positive partition denominator `Z(q)` and event numerator `E_r(q)`. The model-relative contact-event probability is `P_r(q)=E_r(q)/Z(q)`. For two sequences A and B, because `Z_A(q),Z_B(q)>0` for all q>=1,

    sign(P_A(q)-P_B(q))
      = sign(E_A(q)*Z_B(q) - E_B(q)*Z_A(q)).

**The right-hand side is a finite integer polynomial in t=q-1.** In contrast to the partition function, probability optima generally have their own crossing patterns, and positivity of partition-polynomial coefficients does **not** imply event-probability dominance.

## Complete all-parameter reference experiment

The source is the independently exhaustive length-12 anchored 2D HP walk atlas:

- **30,073** complete conformations, with rotations normalized by the first edge and reflections counted separately.
- **4,096** H/P sequences; their partition polynomials are all independently audited against direct per-map contact scoring, not only a sampling of sequences.
- **13** fixed-H composition classes, k=0 through k=12.
- **Five** contact-event thresholds r=1,2,3,4,5.

The checker completed **all 65 composition-by-event optimization problems** and found, across these separate optimization problems:

- **74** genuine global-optimum transitions as q varies over all real q>1;
- **136** globally optimal probability-function classes, counted separately by optimization problem;
- no unresolved coincident/overlapping pair-root event in this tested dataset.

These are exact results in the specified finite toy model. The statements do **not** imply a phase transition in experimental proteins, native structure prediction, a new physical law, or mathematical originality.

### A five-family event design case

For n=12, fixed k=7 hydrophobic residues and event `K>=4`, there are `binom(12,7)=792` candidate sequences. Coefficientwise exact event/partition cross-polynomial dominance leaves 14 sequences in seven rational probability classes. All seven survive the strict all-q Sturm dominance check; their **25 pairwise positive-root events** contain only **four actual global-maximizer changes**.

The winning sequence representatives (each with its reversal-equivalent partner) in successive q>1 regions are:

| q range | Optimal representative |
|---|---|
| `1<q<q1` | `HPHHPHHHPHPP` |
| `q1<q<q2` | `HPHHPPHHPHPH` |
| `q2<q<q3` | `HPHPHHHPPHPH` |
| `q3<q<q4` | `HPHHHPHPPHPH` |
| `q>q4` | `HPHHPPHHHPPH` |

All four crossover roots have exact rational isolating brackets printed by the reproduction script, and lie approximately at

    q1 ≈ 1.289018116
    q2 ≈ 4.545797055
    q3 ≈ 9.189805457
    q4 ≈ 27.577772154

At each crossover the two neighboring maximum probability classes tie; outside the isolated roots the winner is fixed across its full real parameter interval. The 21 other pairwise crossing/tangency events cannot change the optimum. This conclusion follows from an exhaustive exact pairwise root inventory, **not** from sampling the listed parameter values.

A separate k=8,r=4 event problem also yields five globally optimal classes and four switches. Both outcomes illustrate why partition optimization is not a substitute for event-probability optimization.

## Exact certificate algorithm

The new `src/foldbounds/event_probability_envelope.py` provides the following integrity-preserving computational structure:

1. For each sequence, reconstruct exact integer coefficient polynomials for `E_r(1+t)` and `Z(1+t)`. Reject negative/nonintegral coefficients, a nonpositive partition count, inconsistent q=1 conformation counts, or event coefficients exceeding partition coefficients.
2. Define the cross-difference polynomial `D_AB(t)=E_A(1+t)Z_B(1+t)-E_B(1+t)Z_A(1+t)`. Since q>=1 implies t>=0, a **nonnegative coefficient vector** for D is a sufficient, checkable guarantee of `P_A(q)>=P_B(q)` for the entire continuum.
3. Apply an online sufficient-condition antichain filter, then **recheck every original sequence** against the final antichain.
4. Group exactly equal **rational probability functions**, even if their numerator/denominator representations differ.
5. Use exact Sturm sequences and rational arithmetic to identify pairwise strictly dominated probability classes (no positive root and fixed strict sign), then **recheck direct all-q dominance** of every eliminated class by a final survivor.
6. Construct a complete pairwise positive-root inventory for the final classes. Isolate every distinct root into disjoint rational brackets; bounded refinement attempts to separate roots. **Fail closed** on ambiguous root coincidence/overlap.
7. Evaluate exact rational probabilities once per root-free open interval. The absence of any pairwise equality root in that interval proves its winner cannot change there. Classify each single-pair root bracket as a genuine global switch, a global tangent tie, or an irrelevant loser/loser crossing.

This approach works even when event probabilities individually approach 1 as q becomes large; comparing `E_A/Z_A` directly with machine floats in that regime would obscure small but exact differences.

## Reproduction and assurance

```bash
python -m pip install -e '.[test]'
python -m pytest tests/test_event_probability_envelope.py
python experiments/event_probability_design_12.py
```

The dedicated tests:

- independently re-enumerate **all 30,073 full walks**;
- audit **all 4,096** contact spectra and exact partition polynomials against the separate geometry oracle;
- certify all **65** problem instances and assert **74** genuine switches and **136** winning rational probability classes;
- compare every candidate in every composition at independent exact rational q values 2, 10, and 50 against direct contact-spectrum event/partition ratios;
- verify the five k=7,r=4 optimal representatives by direct contact-spectrum scoring;
- reject malicious spectra with impossible event/partition counts and invalid budgets, and verify correct behavior for events impossible for the entire sequence class.

The exact q witness checks are software regressions, **not** the basis of continuum completeness. Complete real positive-root counting plus direct dominance/coverage postconditions are the mathematical grounds for the continuum statement, conditional on correct Python execution.

### Scientific and formal limits

- This remains an anchored **2D HP lattice** experiment, not a physically calibrated protein force field, kinetic model, or proof of native-fold probability.
- q is the model's dimensionless exact Boltzmann contact weight `exp(beta*epsilon)`, not a Kelvin temperature.
- Even an exactly optimal contact-rich ensemble need not be a uniquely folded or experimentally stable protein.
- Classical HP enumeration, polynomial comparison and Sturm root isolation are prior art. No novelty/publication claim is made without deeper literature work and a meaningful experimental benchmark.
- The code is **not Lean-kernel-executed**, and its Python implementation has not been formally refined to the Lean logic. Existing Lean CI passing only establishes what its Lean jobs actually check.

The next meaningful research step is to replace full geometry enumeration by separately certified frontier/partition enclosures in richer coarse-grained physical models and prospectively validate held-out experimental structures.
