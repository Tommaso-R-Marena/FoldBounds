# Exact multiway phase certification of all 65,536 length-16 HP sequences

**Research scope:** This entirely additive branch builds on the exact geometry atlas (#44), all-sequence positive polynomial transform (#47), exact positive-root Sturm checker (#48), and two-finalist skyline (#52). No existing Lean, Python model-verifier, CI, motion-model, blind-test, or other agent file is modified. Keep the PR draft until hosted tests and review pass.

## Main result

An independently full-walk-audited length-16 anchored 2D HP geometry model has **1,604,149 conformations** and **31,940 distinct contact masks**. The total H/P design space is `2^16 = 65,536` sequences. For **every composition k=0,...,16**, the new checker certifies which sequence equivalence classes maximize the exact partition function `Z(q)` for **all real q>1**, including all isolated global switches and global tangencies. At q=1, every sequence has the same Z and therefore ties within every composition class.

Across 17 fixed-H composition classes:

- **65,536** total sequences audited against their final coefficientwise skyline
- **25** global-optimum switches over q>1
- **42** distinct globally optimal partition-polynomial classes, counted separately by composition
- All composition classes fully classified; no unresolved overlapping-root group in this particular dataset

The full phase diagram is reproducible with `python experiments/multiway_all_compositions_16.py`.

### Per-composition summary

| H residues k | Sequences | Global switches | Winning polynomial classes |
|---:|---:|---:|---:|
| 0 | 1 | 0 | 1 |
| 1 | 16 | 0 | 1 |
| 2 | 120 | 0 | 1 |
| 3 | 560 | 1 | 2 |
| 4 | 1820 | 0 | 1 |
| 5 | 4368 | 3 | 4 |
| 6 | 8008 | 1 | 2 |
| 7 | 11440 | 2 | 3 |
| 8 | 12870 | 1 | 2 |
| 9 | 11440 | **4** | **5** |
| 10 | 8008 | 3 | 4 |
| 11 | 4368 | 3 | 4 |
| 12 | 1820 | 1 | 2 |
| 13 | 560 | 2 | 3 |
| 14 | 120 | 2 | 3 |
| 15 | 16 | 2 | 3 |
| 16 | 1 | 0 | 1 |
| **Total** | **65536** | **25** | **42** |

### Five global winners for k=9

Among the 11,440 nine-H / seven-P sequences at length 16, the coefficientwise skyline retains 20 sequences in 10 distinct polynomial classes. Exact positive-root dominance eliminates 3 classes globally, leaving 7 possible finalists. The complete pairwise root atlas contains 21 distinct positive-root events, only **four** of which change the global maximum. The remaining 17 events are crossings among nonwinning polynomial classes.

The winners in successive open q-regions (each also includes the exact reversal-equivalent sequence) are:

| Region | One optimal representative |
|---|---|
| `1<q<q1` | `HPPHHPHHPHHPHPPH` |
| `q1<q<q2` | `HPPHPHHHPHHPHPPH` |
| `q2<q<q3` | `HPPHPHHHHPHPHPPH` |
| `q3<q<q4` | `HPHPHHHHPHPPHPPH` |
| `q>q4` | `HHHHPHPHPPHPHPPH` |

The exact rational isolating intervals printed by the script enclose four distinct irrational/algebraic crossover weights approximately:

- q1 ≈ **3.0899613**
- q2 ≈ **11.3236405**
- q3 ≈ **12.7490368**
- q4 ≈ **18.7960383**

At each crossover, the two adjacent optimal polynomial classes tie. No other model sequence is optimal in any open region between them. All statements are relative to the dimensionless lattice contact weight `q>=1`; **q is not temperature in Kelvin**.

## Why this constitutes a complete continuum certificate

The core algorithm in `src/foldbounds/multiway_envelope.py` is exact rational/integer:

1. Enumerate all `binom(n,k)` sequences at one fixed hydrophobic composition. Compute the coefficient polynomial `Z_S(1+t)` in `t=q-1>=0` from the previously audited exact contact geometry.
2. Maintain a coefficientwise antichain of nondominated partition polynomials; then **independently recheck that every original sequence is coefficientwise dominated by a member of the FINAL antichain**. Nonnegative powers of t mean such dominance holds for all real q>=1.
3. Group identical polynomials (including reflection/reversal-related sequences); retain every tied sequence name, not just one representative.
4. For each distinct class pair, compute the exact integer polynomial difference, remove factors of t, and use complete **Sturm root counting** and exact rational sign isolation. Prune a class only if another has an **everywhere strictly greater** partition polynomial for all q>1.
5. **Recheck that every pruned class has a direct everywhere-strict dominator among FINAL surviving classes.** An intermediate eliminated chain is never trusted as the sole evidence.
6. Collect **every positive real root of every final pair difference** into exact rational isolating intervals. Refine brackets until roots are pairwise separated. If distinct polynomials have coincident or unresolved overlapping roots, **abstain** rather than claim an ordering.
7. Between the globally sorted root intervals there are **no equalities between any distinct finalists**, hence pairwise ordering cannot change there. Exact rational evaluation at one point in each root-free interval determines its unique maximal polynomial class over that *entire* open interval.
8. Within an isolated bracket only one pair can change relative order. If the global winner changes and both adjacent winners are that pair, the event is a **true unique global switch**. If a winner touches another class at an even-multiplicity root and remains maximal, it is an **isolated global tie**. If neither crossing class can win, mark the pair crossing as irrelevant.

Thus the answer applies to the entire continuum, not just a temperature grid. The exact algebraic crossover root itself can be irrational; a rational isolating bracket plus the polynomial root identity specifies it without floating point. At every isolated global switch the two adjacent winner classes tie at the unique root. At q=1 all compositions have the same anchored-walk count.

**Fail-closed limitations:** If any root brackets overlap and remain ambiguous after bounded refinement (including a shared root of two distinct pair-difference polynomials), the routine raises ValueError rather than issuing a full phase claim. Resource budgets bound enumerated compositions, skyline classes, root events, and refinement attempts. This checker is not a general algebraic-number root equality solver.

## Independent validation and reproduction

```bash
python -m pip install -e '.[test]'
python -m pytest tests/test_multiway_envelope.py
python experiments/multiway_all_compositions_16.py
```

The dedicated test module does all of the following:
- Re-enumerates **all 1,604,149 anchored self-avoiding walks** independently of the symmetry-reduced contact atlas.
- Checks all 17 compositions, every one of the **65,536** sequence polynomials against the skyline coverage postcondition, every discarded class against a final strict all-q dominator, and the **25** global switches / **42** winning polynomial classes.
- Re-evaluates all 65,536 exact partition polynomials at independent rational witness weights q=2, 10 and 30 and matches the winning sequences for every composition. These checks are software regressions, **not** the source of continuum completeness.
- Checks all five k=9 winning classes by independently scoring complete contact maps and summing exact rational Boltzmann weights.
- Includes synthetic adversarial polynomial suites (explicitly NOT physical HP geometries) for three-way switching, an irrelevant pair crossover, an even-multiplicity global tangency, coincident roots that must cause abstention, and invalid budget parameters.

The local prototype also reproduced all 17 composition envelopes using an **independently written NumPy-based** integer edge-subset/zeta transform over separately enumerated full anchored walks. The committed implementation uses only the standard Python library.

## Research-integrity boundary

This is a **substantial exact combinatorial optimization capability** for a specified small lattice HP model, not an experimentally validated result on native folds or sequence fitness. The use of contact-map enumeration, graph-subset expansions, coefficientwise polynomial dominance, and Sturm root isolation relies on established methods; this document makes **no claim to have discovered a new mathematical theorem or physical protein-folding mechanism**.

The core polynomial and Sturm arguments are mathematical. Their *Python execution*, parsing, data production, and the 65,536-sequence computation are **not kernel-verified Lean theorems**. The Lean workflow may pass because no Lean code changed, but that does not itself make these new Python certificates Lean-certified. Any step toward real protein claims requires richer physical modeling, independent prospective benchmarks, and scientific prior-art review.
