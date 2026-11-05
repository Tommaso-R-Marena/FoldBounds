# Universal HP sequence landscapes: exact nonnegative interaction polynomials

**Scope:** This additive extension consumes the isolated contact-map atlas from PR #44. No changes to the existing Lean proofs, original Python checker, CI configuration, protein-motion experiments, or other research branches.

**Prior-art boundary:** Contact-set reuse for HP sequence enumeration is not new: Schiemann, Bachmann and Janke, *J. Chem. Phys.* 122, 114705 (2005), [doi:10.1063/1.1814941](https://doi.org/10.1063/1.1814941). Graph-subset expansions, subset-zeta transforms, and Descartes' rule of signs are likewise established. This contribution is a carefully tested FoldBounds implementation and a useful research oracle, **not a claim of a new scientific theorem, biological prediction, or original protein-folding breakthrough**.

## Positive all-sequence partition identity

For a geometry with nonconsecutive lattice contact edges E, let S be the H-residue subset, and let q>=1 be the rational HP Boltzmann contact weight. Put t=q-1>=0. For binary H indicators x_i,

    q^(x_i*x_j) = 1 + t*x_i*x_j.

Expanding over all subsets F of the contact edges of every geometry (with its exact conformation multiplicity) yields

    Z_S(1+t) = sum_{T subset S} sum_{k>=0} C[T,k] * t^k,

    C[T,k] = sum_G multiplicity(G)
                    * #{F subset E(G) : |F|=k, vertices(F)=T}.

**All C[T,k] are nonnegative integers** and C[empty,0] is exactly the anchored-walk count. The program accumulates edge-subset coefficients once, then computes every one of the 2^n sequence polynomials using the Boolean subset-zeta transform. No numerical interpolation, independence assumption, floating-point arithmetic, or repeated walk enumeration is needed for subsequent sequence queries.

This identity holds for every real q>=1, even though the Python evaluator accepts exact rational q only. Exact density-of-states coefficients N_k are recovered from the t-basis via binomial inversion:

    N_k = sum_{j>=k} (-1)^(j-k) * binom(j,k) * C_total[j].

Moments follow from exact derivatives of Z(q), with mean K=q*Z'/Z and variance K=q*Z'/Z+q^2*Z''/Z-(q*Z'/Z)^2.

## Global guarantees for higher-order mutation interactions

For initial H-set S and disjoint new H-residue set U, define the finite mixed difference

    Delta_U Z(S) = sum_{V subset U} (-1)^(|U|-|V|) Z_(S union V).

The polynomial identity makes this equal to a nonnegative-coefficient sum,

    Delta_U Z(S) = sum_{W subset S} sum_k C[U union W,k] * (q-1)^k >= 0.

Thus **every order of hydrophobic-addition mixed difference of the partition sum** is nonnegative for all q>=1, not merely at a sampled grid of contact weights. It is strictly positive for q>1 whenever a relevant positive-order coefficient is present. This theorem does **not** imply the same sign for event probabilities, log Z, physical fitness, or experimentally observed protein mutations.

Example, length 7 with baseline PPPPPPP and additions at zero-based indices 0,2,5:

    Delta_{0,2,5} Z = 2*(q-1)^2.

The coefficient describes an emergent three-residue *partition-sum interaction from pair contacts*, not a literal three-body physical force.

## Exact continuum sequence-order certificates

Every sequence of a fixed length has the same Z at q=1. The separate `sequence_crossing.py` module subtracts two exact integer polynomials in t=q-1 and removes the common power of t. It then uses the number of coefficient sign variations:

- Zero variations (nonzero difference): fixed partition ordering **for every q>1**.
- One variation: exactly one positive real root by Descartes' rule of signs and continuity; the checker returns an exact rational isolating interval or the exact rational root itself.
- More than one variation: **abstain**, without guessing uniqueness.

A concrete same-composition example has ten residues and three H residues each:

    A = HHPHPPPPPP: Z_A = 4067 + 692*t.
    B = PHPHPPHPPP: Z_B = 4067 + 542*t + 24*t^2.
    Z_A - Z_B = t*(150 - 24*t).

Therefore the **unique crossover** for q>1 occurs at t=25/4, or **q=29/4**. A has larger Z for 1<q<29/4; B has larger Z for q>29/4. Here q is a dimensionless contact weight, **not Kelvin temperature**. This is a partition-sum ordering within a toy lattice model, not an experimental sequence-fitness comparison.

## Independent checking and resource limits

`audit_landscape` checks the entire coefficient table by a **different route**: it first independently enumerates the original walk geometry, then recomputes each sequence's t-basis coefficient directly as `sum_G multiplicity(G)*binom(K_G,k)`, and finally applies an inverse Boolean zeta transform to verify all interaction coefficients. It checks *every* sequence at allowed small lengths and fails closed on resource limits. This is valuable independent Python cross-validation, **not** a Lean-kernel proof.

`build_landscape` defaults to independently auditing its input atlas for n<=12; longer runs require explicit opt-out together with an external independent audit. Both the subgraph-expansion and audit routines reject over-budget runs without returning partial tables. Manually constructed dataclass instances should not be trusted without an audit.

Tests compare every H/P sequence through n=8 against separate exact contact scoring at multiple q, verify signed inclusion-exclusion identities for every baseline and up to three simultaneous added residues at n=7, check moment/spectrum identities, reject forged polynomial tables, and test rational and irrational sequence-crossover cases.

Local reproducibility:

```bash
python -m pytest tests/test_contact_map_atlas.py tests/test_universal_landscape.py tests/test_sequence_crossing.py
python experiments/sequence_crossover_demo.py
python experiments/universal_landscape_benchmark.py
python experiments/universal_landscape_benchmark.py --lengths 14
```

In the local Python 3.13 benchmark, at n=14 the independently audited geometry contained **220,375 anchored walks** represented by **5,953 contact maps**. All **16,384** sequence partition sums at q=2 agreed with the direct per-sequence contact-map method, and the polynomial construction plus all queries took about **51x less elapsed time** than that method. This comparison excludes the common initial geometry-enumeration cost, is dependent on the machine and query batch, and is **not** an improvement for arbitrarily long proteins. Benchmark scripts verify equality before reporting any speedup.

**Not certified:** Python executable semantics, input parsing, and this transform are not formally refined to the Lean kernel. No claims of novel mathematical discovery, publishability, experimental biological accuracy, or scalability to full-size real proteins are made.
