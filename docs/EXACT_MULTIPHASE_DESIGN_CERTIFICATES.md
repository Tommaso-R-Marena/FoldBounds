# Exact multi-crossover and composition-optimality certificates

**Parallel scope:** This research extension is stacked on PR #47 (universal sequence polynomials), which is stacked on #44 (contact-map atlas). Every change in this PR is a *new file*. The original Python algorithms, Lean, CI, molecular-motion research, and other agents' branches are unchanged.

## 1. Whole-continuum phase comparisons

The earlier `sequence_crossing.py` in #47 proves fixed ordering or a unique crossover when there are zero/one polynomial sign variations. Its conservative abstention for multiple variations remains intact. The separate `exact_phase.py` lifts that restriction using **Sturm's theorem and exact rational arithmetic**, counting and isolating all distinct positive roots, including even-multiplicity tangencies.

For two fixed-length sequence partition functions, compute the integer polynomial `D(t)=Z_A(1+t)-Z_B(1+t)`, with `t=q-1`. Remove factors of t representing equality at q=1; square-free-reduce by exact rational polynomial gcd; construct the signed Euclidean remainder (Sturm) sequence. A strict rational Cauchy bound encloses all positive roots. Differences of Sturm sign variations count roots between non-root rational endpoints. Recursive subdivision yields exact rational root brackets and definite signs in every intervening root-free cell. No temperature sampling or floating-point root solver is involved.

### Event-probability comparison

If `E_A(q)` and `E_B(q)` are unnormalized event masses for `K>=threshold`, their probabilities are `P_A=E_A/Z_A` and `P_B=E_B/Z_B`. Because both fixed-universe partition functions are strictly positive for q>=1, the sign of the probability difference is exactly the sign of the integer polynomial

    D(q) = E_A(q)*Z_B(q) - E_B(q)*Z_A(q).

Thus `certify_sequence_comparison(...,event_at_least=threshold)` gives exact root-count and sign-cell certificates for event probabilities as well as partition sums.

### Genuine two-crossover case missed by the one-variation test

For n=12, with four H residues in each sequence:

    A = PPPPPHHHPPHP
    B = PHPPHHPPPPHP

Independently checked exact contact spectra from the 30,073-walk anchored-SAW universe:

    A: N_0=26063, N_1=3800, N_2=210
    B: N_0=25639, N_1=4322, N_2=108, N_3=4

For event K>=1, the exact probability-comparison numerator in t=q-1 is

    D(t) = -12750952 - 10390070*t + 2256630*t^2 - 104252*t^3.

**Exactly two distinct positive roots** were isolated with the exact Sturm count:

    q_1 in (9.8938717842, 9.8938791752)
    q_2 in (14.7520341873, 14.7520415783)

These decimal values describe exact rational enclosure endpoints, not numerical-only evidence. For all q>=1: `P_B>P_A` below q1; `P_A>P_B` between q1 and q2; `P_B>P_A` above q2. Checks at q=1,2,10,12,16,20 separately reproduce signs by exact ratio-of-spectra arithmetic. No biological inference or Kelvin temperature interpretation follows from this toy model.

## 2. All-q fixed-composition partition maximizers

The positive partition polynomial from PR #47 is expressed in powers of `t=q-1>=0`. Coordinatewise coefficient dominance of Z_A over Z_B guarantees `Z_A(q)>=Z_B(q)` for every real q>=1; at q>1 the inequality is strict when any positive-degree coefficient strictly exceeds the other.

`fixed_composition_dominance.py` exhaustively compares the entire `binom(n,k)` sequence class with precisely k hydrophobic residues. `universally_dominant` contains guaranteed all-q partition maximizers. `coefficientwise_nondominated` is a conservative remaining set; lack of a global coefficientwise winner **does not rule out** a true global partition champion.

| H residues | Length-12 sequences | Coefficientwise survivors | All-q maximizers certified |
|---:|---:|---:|---|
| 3 | 220 | 4 | none from sufficient condition |
| 4 | 495 | 12 | none from sufficient condition |
| 5 | 792 | 2 | `HPPHPHPPHPPH` and `HPPHPPHPHPPH` (same polynomial) |
| 6 | 924 | **1** | **`HPPHPHHPHPPH`**, strictly better for all q>1 |

At n=14, k=6, an independent full-walk geometry audit verified **220,375 walks**, compressed into **5,953 contact maps**. All **3,003 six-H sequences** reduce by coefficientwise dominance to just two globally relevant candidates:

    A = HPPHPHPPHPHPPH: Z_A(t) = 220375 + 140328*t + 41454*t^2 + 6812*t^3 + 620*t^4 + 36*t^5 + 2*t^6
    B = HPPHPPHHPPHPPH: Z_B(t) = 220375 + 142046*t + 35966*t^2 + 5524*t^3 + 640*t^4 + 56*t^5 + 6*t^6

Every eliminated sequence is dominated by **at least one** survivor through transitivity of finite coefficientwise dominance. Thus the global maximum over all 3,003 sequences must be attained by A or B at every q>=1. The exact difference

    Z_A-Z_B = t*(-1718 + 5488*t + 1288*t^2 - 20*t^3 - 20*t^4 - 4*t^5)

has **exactly two distinct positive roots**, proved with the full Sturm count. This completes the **global parameter phase diagram**, not merely a pruning list:

- B is the unique global maximizer for 1<q<q1;
- A is the unique global maximizer for q1<q<q2;
- B is the unique global maximizer for q>q2;
- A and B tie at the two crossover roots; all 3,003 sequences tie at q=1.

The two exact roots are isolated in rational intervals corresponding approximately to **q1≈1.293014** and **q2≈7.520025**, with the exact brackets produced by the demo and regression test. This is complete across **all real q>1**, not a finite temperature grid. Its statements remain model-relative and conditional on the independently audited geometry and correctly executed Python.

These results concern the full partition sum, not a biological native-structure design objective, event probability, or measured fitness.

## 3. Reproduce and audit

Run after installing the repository in editable mode with test dependencies:

```bash
python -m pytest tests/test_exact_phase_and_dominance.py
python experiments/exact_phase_design_demo.py
```

Tests use independent exact spectrum ratios at witness q values; randomized products with known rational roots and multiplicities; root-count and budget negative controls; and exact all-sequence composition checks. Local additional diagnostics compared the same Sturm code against SymPy's root counting for **400 randomized integer polynomials**, with agreement, and ran 200 separately factored random-polynomial checks.

The input contact atlas is separately audited by full independent SAW enumeration (#44), and the universal sequence polynomial is separately audited by a direct-score inverse transform (#47). This implementation and its root-checking algorithms are in **Python**, not verified by the Lean kernel or refinement-proven. Exact HP sequence enumeration and designability are well established (Irbäck & Troein 2002; Schiemann et al. 2005), as is Sturm root counting. This work does **not** establish a new physical folding law, originality of Sturm's theorem, experimental predictive success, or a research publication.
