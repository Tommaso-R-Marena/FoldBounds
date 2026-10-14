# Summary of changes for run 8becb5ac-4678-49a8-8c98-84d3faa3f7b7
I built a Lean 4 + Mathlib library for certified thermodynamic bounds on finite conformation spaces, set up as the `FoldBounds` package. It builds cleanly with `lake build`, with no `sorry`, `admit`, new axioms or `@[implemented_by]`. `FoldBounds/AxiomAudit.lean` checks 34 headline theorems, and each uses only `propext`, `Classical.choice` and `Quot.sound`.

**Deliverables**
- `README.md`: the technical report.
- `docs/THEOREMS.md`: the exact statements of all 99 public declarations, extracted automatically from the sources.
- Lean files under `FoldBounds/`.

**The seven requirements**
1. **Positivity of Z** (`Z_pos`): holds for every real β when the space is nonempty.
2. **0 ≤ P(A) ≤ 1** (`prob_mem_Icc`): also proved are the complement rule, additivity over disjoint sets, and `prob_eq_ratio`, which writes P(A) as \(Z_A/(Z_A+Z_{A^c})\).
3. **Sound bounds on Z from partial enumeration**, for β ≥ 0. There are three forms: a generic tail form, a uniform energy-bound form, and a per-block energy-interval form (`Z_mem_Icc_of_blocks`).
4. **Bounds on P(A)**:
   - A naive bound from separate intervals.
   - A coupled bound that uses the numerator appearing in the denominator: \(a_L/(a_L+b_U) \le P(A) \le a_U/(a_U+b_L)\).
   - An end-to-end version built directly from enumeration (`prob_mem_of_blocks`).
   - The coupled bound is never worse than the naive one, strictly better on the upper side when the numerator interval has positive width, and always within [0, 1].
   - Its endpoints are exactly the minimum and maximum of \(a/(a+b)\) over the interval box (`coupled_lower_isLeast`, `coupled_upper_isGreatest`).
5. **Refinement**:
   - Intersecting valid intervals never invalidates the enclosure (`prob_enclosure_refine`), and tighter intervals give tighter enclosures.
   - Finer blocks or more enumerated conformations never loosen the bounds on Z.
   - Width bounds show the enclosures collapse as block energy gaps shrink (`prob_width_of_blocks`).
6. **Free energy**, assuming β > 0: certified Helmholtz bounds computed from enumeration. Positivity of the lower Z bound is proved, not assumed. Also proved: \(F \le E(c)\) for every conformation, \(F \ge E_{\min} - \log|C|/\beta\), \(P(A) = e^{-\beta(F_A-F)}\), and certified bounds on free-energy differences.
7. **Strengthening**:
   - A counting (density-of-states) formula, and symmetry reduction for an energy-preserving group action (`Z_eq_sum_orbits`), with bounds certified at the level of reduced states.
   - Adaptive block partitions.
   - Two-domain compositional bounds for Z and F.
   - The Gibbs–Bogoliubov inequality, including its two-domain form.

**Research results** (each is formally proved)
- **Data-optimality** (`Optimality.lean`): both endpoints of the enumeration-based P(A) bound are reached by energy landscapes that match the same certificate data. So no sound method that sees only that data can do better.
- **Convergence certificate**: a guaranteed probability resolution follows from the block energy gaps.

The report classifies results as established, newly derived, conjectural, or needing a literature check before any novelty claim. I found no asymptotic speed-up over known algorithms. The proved advantages are tighter bounds for the same data, the optimality result, and the stopping rule. Four conjectures are listed as unproved: greedy block splitting, multi-domain composition, symmetry-aware optimality, and a hybrid bound.

**Real-number proofs vs. computation**
All main theorems are about real numbers. `Certificate.lean` shows that exact rational (ℚ) arithmetic gives exactly the real formulas. As a worked example, `toy_prob_mem` proves \(0.6652 \le P \le 0.6654\) for energies (0, 1, 2) at β = 1. Not covered: a verified `exp` for general arguments, floating-point arithmetic, and any runnable or benchmarked checker.

**Assumptions**
- Energy certificates (lower and upper energy bounds) are hypotheses; producing them is outside the formalization.
- β ≥ 0 is assumed only where energy bounds become weight bounds.
- β > 0 is assumed for all free-energy statements.

**Project changes**: I renamed the Lake package to `FoldBounds`. The original template `RequestProject/Main.lean` is still there but is not part of the library.