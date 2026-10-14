This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# FoldBounds

Machine-checked bounds on Boltzmann partition functions, equilibrium probabilities and
Helmholtz free energies over **finite conformation spaces**, written in Lean 4 + Mathlib.

* Lean `v4.28.0`, Mathlib `v4.28.0` (see `lean-toolchain`, `lake-manifest.json`).
* Build: `lake build` builds the library and `FoldBounds/AxiomAudit.lean`, which prints the
  axioms used by each headline theorem.
* There is no `sorry` or `admit`, no new `axiom`, and no `@[implemented_by]`. Every audited
  theorem uses only the standard axioms `propext`, `Classical.choice` and `Quot.sound`.
* **Exact statements** of all 99 public declarations: [`docs/THEOREMS.md`](THEOREMS.md)
  (extracted automatically from the sources).

This README is the technical report. Formulas are written in plain Unicode, with
`Z_A = ∑_{c∈A} e^{-βE(c)}` and so on.

---

## 1. Setting and conventions

| Lean | Mathematics |
|---|---|
| `C : Type*`, `[Fintype C]` | finite conformation space |
| `E : C → ℝ` | energy |
| `β : ℝ` | inverse temperature `1/(k_B T)` |
| `weight β E c = exp (-β * E c)` | Boltzmann weight |
| `Z β E = ∑ c, weight β E c` | partition function |
| `ZOn β E A = ∑ c ∈ A, weight β E c` (`A : Finset C`) | partial partition function `Z_A` |
| `prob β E A = ZOn β E A / Z β E` | equilibrium probability `P(A)` |
| `freeEnergy β E = -log (Z β E) / β` | Helmholtz free energy `F` |
| `freeEnergyOn β E A = -log (ZOn β E A) / β` | macrostate free energy `F_A` |
| `blockBound β blk e T = ∑ c ∈ T, exp (-β * e (blk c))` | block bound (see §3) |

**Assumptions.** Every assumption appears as an explicit hypothesis or instance argument:

* `[Nonempty C]` wherever positivity of `Z` is needed.
* `0 ≤ β` only when an energy bound is turned into a weight bound. `Z > 0` and
  `0 ≤ P(A) ≤ 1` hold for **every** real `β`.
* `0 < β` (positive, finite temperature) for every free-energy statement. Lean defines
  `x / 0 = 0`, but no theorem relies on that convention.
* Certificates such as "`ℓ (blk c) ≤ E c` for every unexplored `c`" are hypotheses. They are
  the interface to external energy-bounding procedures (for example pairwise lower bounds or
  dead-end elimination), which this library does not formalize.

No theorem takes its own conclusion as a hypothesis. For example, the generic tail theorem
`Z_mem_Icc_of_tail` assumes a bound on `Z_{Sᶜ}` and proves a bound on `Z`. The
block-certificate theorems assume only pointwise energy bounds and derive the bounds on
`Z_{Sᶜ}` themselves.

## 2. Module map and dependencies

```
Basic ─► PartitionBounds ─► Probability ─► FreeEnergy ─┬─► Symmetry
                                   │                   ├─► Compositional ─► Variational
                                   └─► Optimality      └─► Certificate
FoldBounds.lean imports all of them; AxiomAudit.lean prints axioms.
```

All modules import `Mathlib`. The Mathlib facts used are: finite sums (`Finset.sum_*`,
`Finset.sum_fiberwise`, `Finset.sum_mul_sum`, `Fintype.sum_prod_type`), `Real.exp`/`Real.log`
monotonicity and algebra, `Real.add_one_le_exp`, `Real.exp_one_gt_d9`/`Real.exp_one_lt_d9`,
and `MulAction.orbitRel`.

## 3. Results, by requirement

Each item gives the Lean name, followed by the statement in words. Exact Lean text is in
`docs/THEOREMS.md`.

### (1) Positivity of `Z` — `Basic.lean`
* `Z_pos [Nonempty C] : 0 < Z β E`. This holds for every real `β`.
* Also proved: `ZOn_pos`, `ZOn_nonneg`, `ZOn_mono`, `ZOn_le_Z`, and
  `ZOn_add_ZOn_compl : Z_A + Z_{Aᶜ} = Z`.

### (2) `0 ≤ P(A) ≤ 1` — `Basic.lean`
* `prob_mem_Icc : prob β E A ∈ Set.Icc 0 1`. No hypotheses are needed: if `C` is empty then
  `Z = 0` and `P = 0`.
* Also proved: `prob_univ = 1`, `prob_empty = 0`, `prob_compl : P(Aᶜ) = 1 − P(A)`,
  `prob_union_of_disjoint` (additivity), `prob_pos`, and
  `prob_eq_ratio : P(A) = Z_A / (Z_A + Z_{Aᶜ})`.

### (3) Sound bounds on `Z` from partial enumeration — `PartitionBounds.lean`
Let `S` be the explored set. The unexplored conformations are labelled by blocks
`blk : C → ι`, and each block `i` carries energy bounds `ℓ i ≤ E c ≤ u i`.
* `Z_mem_Icc_of_tail`: if `L ≤ Z_{Sᶜ} ≤ U`, then `Z ∈ [Z_S + L, Z_S + U]`.
* `Z_mem_Icc_of_blocks` (`0 ≤ β`): `Z ∈ [Z_S + blockBound u Sᶜ, Z_S + blockBound ℓ Sᶜ]`.
* `Z_mem_Icc_of_uniform` (`0 ≤ β`):
  `Z ∈ [Z_S + |Sᶜ|·e^{-βEmax}, Z_S + |Sᶜ|·e^{-βEmin}]`.
* `blockBound_eq_sum_card` gives the computable form
  `blockBound = ∑_i #{c ∈ T | blk c = i} · e^{-β e_i}`. It needs only block sizes and one
  number per block.

### (4) Certified bounds on `P(A)` — `Probability.lean`
* `prob_mem_naive`: from `aL ≤ Z_A ≤ aU` and `0 < Zlo ≤ Z ≤ Zhi`, we get
  `aL/Zhi ≤ P(A) ≤ aU/Zlo`.
* `prob_mem_coupled` (the numerator also appears in the denominator): from
  `Z_A ∈ [aL,aU]`, `Z_{Aᶜ} ∈ [bL,bU]`, `aL, bL ≥ 0`, we get
  **`aL/(aL+bU) ≤ P(A) ≤ aU/(aU+bL)`**.
* `prob_mem_of_blocks`: the coupled enclosure built directly from enumeration plus block
  certificates. The bounds on `Z_A` and `Z_{Aᶜ}` are split into explored and unexplored parts.
* Comparison with the naive enclosure, using the same data with `Zlo = aL+bL` and
  `Zhi = aU+bU`:
  * `naive_lower_le_coupled_lower` and `coupled_upper_le_naive_upper`: the coupled bounds are
    never worse.
  * `coupled_upper_lt_naive_upper`: the coupled upper bound is *strictly* better when
    `aL < aU` and `aU > 0`.
  * `coupled_bounds_mem_unit`: the coupled enclosure always lies in `[0,1]`. The naive upper
    bound can exceed 1.
* **Optimality** (`coupled_lower_isLeast`, `coupled_upper_isGreatest`): the coupled endpoints
  are exactly the minimum and maximum of `a/(a+b)` over the box `[aL,aU]×[bL,bU]`.

### (5) Refinement — `PartitionBounds.lean`, `Probability.lean`
* `mem_Icc_max_min`: intersecting valid enclosures gives a valid enclosure inside both.
* `prob_enclosure_refine`: given two independently certified interval pairs, the coupled
  enclosure built from their intersections contains `P(A)` and lies inside the original
  coupled enclosure. In short, refinement never invalidates the enclosure.
* `coupled_enclosure_mono`: tighter input intervals give a tighter probability enclosure.
* `blockBound_antitone_refine` and `blockBound_refine_partition`: finer blocks with child
  bounds at least as large as the parent bound can only lower the upper bound on `Z`.
  `blockBound_monotone_refine` is the dual statement for the lower bound.
* `upper_antitone_enumeration` and `lower_monotone_enumeration`: enumerating more
  conformations (`S ⊆ S'`) never loosens either end of the `Z` enclosure.
* Convergence (the width shrinks to zero as information increases):
  * `blockBound_sub_le`: `width ≤ β ∑ e^{-βℓ}(u−ℓ)`. This holds for every real `β`.
  * `blockBound_sub_le_mul`: `width ≤ βδ · blockBound ℓ` when every gap `u − ℓ` is at most δ.
  * `coupled_width_le`: `width(P) ≤ (width(Z_A)+width(Z_{Aᶜ}))/(aL+bL)`.
  * `prob_width_of_blocks`: combining the above,
    `width(P) ≤ βδ · (unexplored upper mass) / (certified lower Z)`.

### (6) Free energy — `FreeEnergy.lean` (requires `0 < β`)
* `freeEnergy_mem_Icc`: from `0 < Zlo ≤ Z ≤ Zhi`, we get `F ∈ [−log Zhi/β, −log Zlo/β]`.
* `freeEnergy_mem_Icc_of_blocks`: the same enclosure computed from enumeration plus block
  certificates. Positivity of the lower end is *proved* (`certifiedZLower_pos`), not assumed.
* `freeEnergy_le_energy : F ≤ E c` for every `c`.
  `le_freeEnergy_of_energy_ge : Emin − log|C|/β ≤ F`.
* `freeEnergy_le_freeEnergyOn : F ≤ F_S` (variational bound from any nonempty explored set).
* `prob_eq_exp_freeEnergy : P(A) = exp(−β(F_A − F))`.
  `freeEnergy_diff_mem_Icc` turns a certified `P(A)`-enclosure `[pL,pU]` with `pL > 0` into
  a certified enclosure of `ΔF = F_A − F`.

### (7) Strengthening: symmetry, counting, adaptive partitioning, composition
* **Counting / density of states** (`Symmetry.lean`): `Z_eq_sum_multiplicity` states that
  if `E = e ∘ f`, then `Z = ∑_d m_d e^{-β e_d}`.
  `Z_mem_Icc_of_reduced` gives a certified enclosure that only uses the multiplicities, the
  exact energies of explored *reduced* states, and energy intervals for unexplored reduced
  states.
* **Symmetry** (`Symmetry.lean`): `Z_eq_sum_orbits` states that for an energy-preserving
  group action, `Z = ∑_{orbits ω} |ω| · e^{-βE(rep ω)}`. `multiplicity_orbit_eq` identifies
  each multiplicity with the orbit size.
* **Adaptive partitioning**: `blockBound` with refinement monotonicity and the width bounds
  from (5). `blockBound_eq_ZOn_of_exact` shows that blocks on which `E` is constant (for
  example orbits) give *exact* block bounds.
* **Composition** (`Compositional.lean`): take `E(x,y) = E₁x + E₂y + J(x,y)` with
  `J ∈ [Jlo, Jhi]`.
  * `Z_coupledEnergy_mem_Icc`: `Z₁Z₂e^{−βJhi} ≤ Z ≤ Z₁Z₂e^{−βJlo}`.
  * `freeEnergy_coupledEnergy_mem_Icc`: `F₁+F₂+Jlo ≤ F ≤ F₁+F₂+Jhi`.
  * `Z_coupledEnergy_zero`: the exact factorisation `Z = Z₁Z₂` when `J = 0`.
* **Variational / Gibbs–Bogoliubov** (`Variational.lean`):
  * `Z_ge_gibbs_bogoliubov`: `Z(E₀+V) ≥ Z(E₀) e^{−β⟨V⟩₀}` for every real `β`.
  * `freeEnergy_le_gibbs_bogoliubov`: `F ≤ F₀ + ⟨V⟩₀`.
  * `refAverage_le`: the bound is never weaker than the worst-case bound.
  * `Z_coupledEnergy_ge_gibbs`: the two-domain case of the bound.

### Data-optimality of the enumeration enclosure — `Optimality.lean`
`prob_lower_endpoint_attained` and `prob_upper_endpoint_attained`: for any valid block
certificate there is an energy landscape `E'` with exactly the same certificate data. It
agrees with `E` on `S` and satisfies `ℓ ≤ E' ≤ u` blockwise off `S`. Its probability
`P_{E'}(A)` **equals** the certified lower (respectively upper) endpoint of
`prob_mem_of_blocks`. As a result, no sound procedure that sees only
`(S, E|_S, blk, ℓ, u)` can return a strictly narrower interval.

## 4. Real-number proofs vs. finite-precision certificates

All theorems in §3 are statements about **real numbers**. `Certificate.lean` connects them
to computation:

* `coupledLowerQ` and `coupledUpperQ` are computable functions on `ℚ`. They use exact
  rational arithmetic with no rounding.
* `ZOn_mem_Icc_of_rat`: rational enclosures of individual weights sum, exactly in `ℚ`, to an
  enclosure of `Z_T`.
* `prob_mem_Icc_of_rat`: rational interval data certified as real inequalities makes the
  *rationally computed* endpoints enclose `P(A)`.
* `toy_prob_mem` is a fully checked instance. It has energies `(0,1,2)`, `β = 1` and
  `A = {0}`, and proves `0.6652 ≤ P(A) ≤ 0.6654`. The true value is about `0.66524`. The
  bounds on `e^{−1}` and `e^{−2}` are derived from Mathlib's 9-digit bounds on `e`.

What is **not** covered:

* There is no verified enclosure of `exp` for general rational arguments. A production
  checker would need a verified interval-`exp` library.
* There is no IEEE floating-point model. Floating-point code would need directed rounding
  and a separate proof.
* No executable checker has been extracted or benchmarked.

The rational layer shows only that the abstract theorems transfer exactly to exact-rational
computations once the weight enclosures are certified.

## 5. Research investigation: classification of results

| Result | Status |
|---|---|
| `Z>0`, `0≤P≤1`, `Z = Z_S + Z_{Sᶜ}`, uniform and block tail bounds, free-energy monotone transfer, `F ≤ E_min`, `F ≥ E_min − log|C|/β` | **Established mathematics** (textbook). The contribution here is the formalization. |
| Density-of-states and orbit formulas for `Z` | **Established** (orbit counting). Formalized. |
| Product factorisation and the `e^{±βJ}` compositional bound | **Established**, elementary. Formalized. |
| Gibbs–Bogoliubov inequality | **Established** (classical statistical mechanics). Formalized here with a tangent-line proof. |
| Coupled ratio enclosure `aL/(aL+bU) ≤ P ≤ aU/(aU+bL)`, with its optimality over the box and strict improvement over the naive ratio | **Elementary consequences derived here and formally proved**. Ratio bounds like this are probably used, at least implicitly, in existing algorithms. **Novelty would need a literature check.** |
| Data-optimality of the enumeration enclosure (endpoints attained by consistent landscapes) | **Derived here and formally proved**. It is simple, but I have not seen it stated as such. **Novelty would need a literature check.** |
| Width and convergence certificate `width(P) ≤ βδ·tail/Z_lo` | **Derived here and formally proved**. Related to the ε-stopping rules of provable partition-function algorithms in computational protein design, such as the K*/BBK* family. **Novelty would need a literature check.** |
| Refinement and enumeration monotonicity | **Elementary, proved here.** They match the anytime-bound behaviour of A*-based methods. Not claimed as new. |

**Does any of this give a nontrivial advantage over known algorithms?** No asymptotic
speed-up was found. Computing `Z` exactly is #P-hard in general for structured energy
models, and none of the results above changes the worst case. The proved advantages are:

1. **Tightness at no extra cost.** The coupled enclosure is never worse than the
   independent-ratio enclosure. It is strictly better whenever the numerator interval has
   positive width, and it is optimal for the given interval data (§3(4), `Optimality.lean`).
2. **A proved stopping and convergence rule.** `prob_width_of_blocks` says how small the
   block energy gaps must be to guarantee a given probability resolution.
3. **Symmetry.** Energy evaluations drop from `|C|` to the number of orbits, with *exact*
   block bounds on orbits. This is classical, but the proved interface lets it compose with
   the other bounds.
4. **Composition.** A two-domain system needs `|C₁|+|C₂|` evaluations instead of
   `|C₁|·|C₂|`. The price is a factor `e^{β(Jhi−Jlo)}`, which Gibbs–Bogoliubov reduces on the
   lower side.

**Conjectures and open directions.** None of these is proved.

* **(C1) Greedy splitting.** Splitting the block with the largest width contribution
  `|B|(e^{−βℓ_B} − e^{−βu_B})` is a natural heuristic. Whether it is optimal at each step, or
  globally within a constant factor, is open. It is plausibly false in general, because the
  achievable bounds after a split depend on the bounding oracle.
* **(C2) Multi-domain composition.** The `k`-domain extension
  `Z ≥ ∏Z_i · exp(−β ∑_{i<j}⟨J_ij⟩)` should follow from `Z_ge_gibbs_bogoliubov` by induction
  over products. It has not been formalized.
* **(C3) Symmetry-aware data-optimality.** The analogue of `Optimality.lean` for orbit-level
  certificates should hold by the same construction. It has not been formalized.
* **(C4) Hybrid enclosure.** Combining a Gibbs–Bogoliubov lower bound with block upper
  bounds inside the coupled enclosure gives valid bounds by the theorems above. How much this
  helps on realistic protein landscapes is an empirical question that has not been tested.

## 6. Limitations

* Energy certificates (`ℓ`, `u`) are hypotheses. Producing them, for example from
  pairwise-decomposable force fields, is outside this formalization.
* Conformation spaces are finite types. Continuous or very large implicit spaces enter only
  through `Fintype` and block labels.
* The numerical example is a toy. No benchmarks on real proteins were run.

## 7. Repository layout

```
FoldBounds.lean                root import
FoldBounds/Basic.lean          definitions; Theorems 1–2
FoldBounds/PartitionBounds.lean  Theorem 3; Z-refinement and width bounds
FoldBounds/Probability.lean    Theorems 4–5; coupled/naive/optimal/width
FoldBounds/FreeEnergy.lean     Theorem 6
FoldBounds/Symmetry.lean       counting and orbit reduction
FoldBounds/Compositional.lean  domain-decomposition bounds
FoldBounds/Variational.lean    Gibbs–Bogoliubov
FoldBounds/Optimality.lean     data-optimality of the enclosure
FoldBounds/Certificate.lean    exact-rational certificates; checked numerical example
FoldBounds/AxiomAudit.lean     #print axioms for headline theorems
docs/THEOREMS.md               exact statements (auto-extracted)
```

(`RequestProject/Main.lean` is the original project template and is not part of the
library.)