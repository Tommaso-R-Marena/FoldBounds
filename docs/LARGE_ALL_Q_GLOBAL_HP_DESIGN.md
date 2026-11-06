# Exact global thermodynamic design of 12,870 model sequences

**Parallel scope:** This draft is stacked on #48, then #47, then #44. It only adds new files; the existing numerical algorithms, Lean proofs, RNase blind experiment, motion code, CI and other agents' branches are unchanged.

## Complete n=16, eight-H sequence phase diagram

For the 16-residue anchored 2D lattice HP model, the independent geometry enumerator gives **1,604,149 conformations** represented by **31,940 distinct contact masks**. There are precisely `binom(16,8)=12,870` H/P sequences with eight H residues.

The universal nonnegative polynomial expansion from #47 produces an exact integer polynomial `Z_S(1+t)` for each sequence in `t=q-1`. Because `t>=0`, if the coefficient vector of A dominates the vector of B in every degree, then `Z_A(q)>=Z_B(q)` for **every real q>=1**. An online antichain (skyline) retains sequences whose coefficient vectors are not strictly dominated. A separate final replay verifies that every one of the original 12,870 sequences is coefficientwise dominated by at least one final survivor.

### First elimination: coefficientwise skyline

12,870 sequences are reduced to **nine** coefficientwise-undominated sequences corresponding to **five** distinct partition polynomials (the others arise from sequence-reversal equivalence).

### Second elimination: whole-continuum Sturm dominance

For each pair of distinct candidate polynomials, form their exact integer difference in `t=q-1`. The standalone Sturm checker from #48 counts **all** positive real roots. If the difference has zero positive roots and a fixed positive sign, one polynomial is strictly better than the other for all real q>1. Three of the five candidate polynomial classes are eliminated by this stronger all-q comparison.

The checker verifies again that each removed class is directly globally dominated by a **final** surviving class. This does not depend on trusting a heuristic ranking or a numerical grid.

**Only two polynomial classes remain** for the entire full-continuum optimum, each representing a pair of reversal-related sequences:

Low-to-moderate q maximizers:

    HPPHPHHPHHPPHPPH
    HPPHPPHHPHHPHPPH

High-q maximizers:

    HPHHPHPPHHPPHPPH
    HPPHPPHHPPHPHHPH

The difference `Z_high - Z_low` is exactly

    t*(-158788 - 95134*t - 23630*t^2 - 3224*t^3 - 160*t^4 + 18*t^5 + 6*t^6).

The Sturm checker proves this polynomial has **exactly one distinct positive root** after the already-known `t=0` factor. The unique global switch is isolated as:

    24049038191/2147483648 <= q* <= 751532857/67108864

That is approximately **q*=11.198709**, with exact rational enclosure width less than 1/100000. For every real **1<q<q***, the low-q pair are the only global partition maximizers among all 12,870 sequences; for **q>q***, the high-q pair are the only global partition maximizers. At q*, the four sequences tie for the maximum. At q=1, all 12,870 sequences have exactly Z=1,604,149 and therefore tie.

This classification is proven from finite polynomial identities and Sturm's theorem applied to independently re-enumerated geometry, rather than guessed from sampled q. It is **not** a result on experimentally measured native folding, or a claim of physically meaningful temperature.

## Algorithm, limits and assurance

`src/foldbounds/extended_global_phase.py` implements the finite exact algorithm:

1. Enumerate each fixed-composition sequence exactly once.
2. Maintain the coefficientwise maximal antichain, retaining identical-polynomial ties.
3. **Recheck every sequence** against the final antichain.
4. Group identical polynomials into equivalence classes.
5. Compute exact Sturm root counts/signs for each candidate-class pair; eliminate only globally strictly dominated classes.
6. **Recheck** that every removed class is dominated by a FINAL class.
7. If exactly one or two polynomial classes remain, return a global phase certificate for all real q>1; otherwise **abstain**.

The new tests check all 12,870 sequences at exact rational witnesses q=2,10,20, reconstruct each finalist polynomial directly from the independently scored contact spectrum, and assert both full-walk atlas replay and exact rational crossover/root counting. These sampled witnesses are regression controls; they are not the basis of the full-continuum proof.

Reproduce (Python 3.11+, no external scientific libraries):

```bash
python -m pip install -e '.[test]'
python -m pytest tests/test_extended_global_phase.py
python experiments/extended_global_phase_16.py
```

The large model requires explicitly raised but bounded resource budgets (`max_length=16`, `max_edge_subsets=3000000`, `max_sequences=20000`). Work-limit failures return no partial certificate.

**Trust boundaries:** The result assumes the anchored two-dimensional lattice HP model, a correct input atlas and the correctly executed Python logic. It is not Lean-kernel verified nor an experimental protein prediction. Known published HP-design literature already contains enumeration and designing-sequence work (Irbäck & Troein 2002; Schiemann et al. 2005). Polynomial root counting by Sturm is classical. No scientific originality or protein-physics breakthrough is claimed. The contribution is a model-specific exact global parameter-optimization and audited test-oracle capability.
