# CEC: Lean 4 formalization of `CEC_THEORY.md` §1–§3

**Status (2026-11-13).** Machine-checked in Lean 4.28.0 / Mathlib v4.28.0 with no proof holes.
Every audited declaration depends only on `propext`, `Classical.choice` and `Quot.sound`; none
depends on `Lean.ofReduceBool`. This document maps the specification to Lean, lists every
deviation from the spec with a counterexample, and states what is **not** established.
Independent review is still pending, as the spec requires.

| File | Content |
|---|---|
| `FoldBounds/CoupledRatio.lean` | §1: exact tests, extremes, dominance, refinement, T1′, rational checker |
| `FoldBounds/CoupledSeparator.lean` | §2: Lemma 2.0, Theorems 2.1–2.6, Corollary 2.2′, separator checker |
| `FoldBounds/CoupledBinding.lean` | §3: Theorems 4.1–4.3 |
| `FoldBounds/CoupledTests.lean` | kernel-evaluated accept/reject tests (not imported by the library) |
| `scripts/check_coupled_lean.py` | source and axiom audit for `FoldBounds/Coupled*.lean`. *Review M1 (2026-11-13): the original parser accepted a forged axiom printed as `«Classical.choice»`. Axiom names are now parsed strictly and the forged log is rejected. The source scan still covers only `Coupled*.lean`; axioms introduced elsewhere are caught only through the strict log parse.* |
| `.github/workflows/coupled-lean.yml` | CI: source audit, build, tests, both axiom checkers |

## Conventions

All three namespaces are `FoldBounds.CoupledRatio`, `FoldBounds.CoupledSeparator` and
`FoldBounds.CoupledBinding`.

**§1.** Regions are `k ∈ s : Finset ι`, with data `L U r R : ι → ℝ`.
- `WellFormed s L U r R` is the §1.1 data condition `∀ k ∈ s, 0 ≤ L k ≤ U k ∧ 0 ≤ r k ≤ R k`.
- `IsAdm s L U r R α β` is the §1.1 admissible set, including `α, β ≥ 0` and `Σα > 0`.
- `Q s α β = Σβ/Σα`, and `QSet` is its image on the admissible set.
- `gLo`/`gHi` are g₋/g₊ in their `min`/`max` form.
- `lamLo = sInf QSet` and `lamHi = sSup QSet`.
- `sepLo = Σ r L / Σ U` and `sepHi = Σ R U / Σ L`.

**§2.** `T`, `IA`, `IB` and `O` are `Fintype`s. Weights are curried: `a : IA → T → ℝ`,
`b : IB → T → ℝ` and `e : T → O → ℝ`.
- `Zsys x e = Σ_{p : I × T × O} x p.1 p.2.1 * e p.2.1 p.2.2`.
- `Phi x τ = Σ_i x i τ` and `Psi e τ = Σ_o e τ o`.
- `rho φA φB τ = φB τ / φA τ` and `Tplus φA = {τ | 0 < φA τ}`.
- `rhoMin`/`rhoMax` are `inf'`/`sup'` of `rho` over `Tplus`.

Lean's `x / 0 = 0` is a junk value. Wherever the spec's value would be `+∞`, an explicit
positivity hypothesis is added. Each such case is flagged below.

## Spec item → Lean declaration

Status: **P** = proved as stated, or with weaker hypotheses (see Deviations);
**Part** = partially formalised; **N** = not attempted.

### §1 Coupled linear-fractional enclosures

| Spec | Lean name(s) | Exact hypotheses | Status |
|---|---|---|---|
| g₋, g₊ and their `ite` forms | `gLo`, `gHi`, `gLo_eq_ite`, `gHi_eq_ite` | `∀ k ∈ s, L k ≤ U k` (`ite` forms only) | P |
| Thm 1.1 (⇐) soundness | `le_Q_of_gLo_nonneg` | `0 ≤ gLo s L U r t`, `IsAdm …` (no data condition) | P |
| Thm 1.1 (⇒) completeness | `gLo_nonneg_of_forall_le_Q` | `WellFormed` | P |
| Thm 1.1 iff | `forall_le_Q_iff` (`forall_le_Q_iff_of_nonempty` = spec-literal) | `WellFormed` (non-emptiness not needed) | P |
| Thm 1.2 (both directions, iff) | `Q_le_of_gHi_nonpos`, `gHi_nonpos_of_forall_Q_le`, `forall_Q_le_iff` (`…_of_nonempty`) | as for 1.1 | P |
| Meaning of `Adm ≠ ∅` | `adm_nonempty_iff` | `WellFormed`; `(∃ α β, IsAdm …) ↔ ∃ k ∈ s, 0 < U k` | P |
| Thm 1.3: existence and vertex attainment | `exists_isLeast_vertex`, `exists_isGreatest_vertex`, `isLeast_lamLo`, `isGreatest_lamHi`, `lamLo_attained_at_vertex`, `lamHi_attained_at_vertex` | `WellFormed`, `∃ k ∈ s, 0 < U k` | P |
| Thm 1.3: `g₋(λ₋) = 0`, `g₊(λ₊) = 0`; test ⇔ `t ≤ λ₋` | `gLo_lamLo_eq_zero`, `gHi_lamHi_eq_zero`, `le_lamLo_iff`, `lamHi_le_iff` | same | P |
| Thm 1.3: sorted-threshold description | — | — | N |
| Thm 1.4 lower: `λ₋ ≥ max(S₋, min{r_k : U_k>0})` | `coupled_dominates_lower` (parts: `sepLo_le_lamLo`, `minRatio_le_lamLo`) | `WellFormed`, `(s.filter (0 < U ·)).Nonempty` | P |
| Thm 1.4 upper: `λ₊ ≤ min(S₊, max{R_k : U_k>0})` | `coupled_dominates_upper` (parts: `lamHi_le_sepHi`, `lamHi_le_maxRatio`) | `WellFormed`, **`0 < Σ L`** (S₊ part only), nonempty filter | P |
| Thm 1.5 lower | `lamLo_eq_sepLo_iff` | `WellFormed`, `0 < sepLo`; `lamLo = sepLo ↔ ∀ k ∈ s, r k * (U k - L k) = 0` | P |
| Thm 1.5 upper | `lamHi_eq_sepHi_iff` | `WellFormed`, `0 < Σ L`; `lamHi = sepHi ↔ ∀ k ∈ s, R k * (U k - L k) = 0` | P |
| Thm 1.6 uniform ratio | `lamLo_lamHi_of_uniform` (also `Q_eq_of_uniform`, no hypotheses) | `WellFormed`, `∃ k ∈ s, 0 < U k`, `∀ k ∈ s, 0 < U k → r k = ρ ∧ R k = ρ` | P |
| Thm 1.7(a) monotonicity | `QSet_mono` (no data condition), `lamLo_mono`, `lamHi_anti`; certificate form `gLo_mono`, `gHi_anti` | coarse `WellFormed`, refined `QSet` nonempty, `L ≤ L'`, `U' ≤ U`, `r ≤ r'`, `R' ≤ R` on `s` (test form: `0 ≤ L' ≤ U'`) | P |
| Thm 1.7(b) splitting | `QSet_refine`, `lamLo_refine`, `lamHi_refine` | map `π : ι' → ι` sending `s'` into `s`; `L k ≤ Σ_{π j = k} L' j`; `Σ_{π j = k} U' j ≤ U k`; `[r' j, R' j] ⊆ [r (π j), R (π j)]` | P (any number of pieces) |
| §1.4 T1′ soundness, lower | `le_Q_of_lowerDualCert` with `IsLowerDualCert`, `IsAdmPoly` | `y ≥ 0` on `c`; `0 ≤ r k - t + Σ_j G j k y j`; `Σ_j b j y j ≤ 0` (`b` = spec's `h`) | P |
| §1.4 T1′ soundness, upper | `Q_le_of_upperDualCert` | `0 ≤ t - R k + Σ_j G j k y j`, rest as above | P |
| §1.4 completeness via LP duality | — | — | N (spec: not needed) |
| Executable checker | `checkLower`, `checkUpper`, `checkEnclosure : List (ℚ×ℚ×ℚ×ℚ) → …`; entries are `(L, U, r, R)` | gate `regionWF`: `0 ≤ L ≤ U`, `0 ≤ r ≤ R` | P |
| Checker soundness | `checkLower_sound`, `checkUpper_sound`, `checkEnclosure_sound` | `checkEnclosure regions lo hi = some bounds` ⇒ `bounds = (lo, hi)`, the admissible set is nonempty, and `Q univ α β ∈ [bounds.1, bounds.2]` for **every** real `ListAdm regions α β` | P |
| Checker exactness | `checkLower_iff`, `checkUpper_iff`, `checkEnclosure_sharp_bounds`, `checkEnclosure_eq_none_of_not_wf` | the checker accepts iff the input is well formed and `t` is a valid bound; accepted endpoints bracket `lamLo`/`lamHi`; malformed input gives `none` | P |

### §2 Separator certificates

| Spec | Lean name(s) | Exact hypotheses | Status |
|---|---|---|---|
| Lemma 2.0 | `Zsys_eq_sum` | none | P |
| Thm 2.1 | `local_enclosure`, `local_enclosure_of_bounds`; spec-literal `local_enclosure_spec` | `a ≥ 0`, `e ≥ 0`, `∀ τ, Phi a τ = 0 → Phi b τ = 0`, `0 < Zsys a e` (spec-literal: `b ≥ 0` and `¬∃ τ, Φ_A τ = 0 ∧ 0 < Φ_B τ`) | P |
| Thm 2.2, closed interval | `sharp_closed` (core: `exists_two_point`, `exists_weight_of_mem_Icc`) | `[Nonempty O]`, `(Tplus (Phi a)).Nonempty` only; the exterior is `≥ 0`, supported on ≤ 2 boundary states, `Z_A = 1`, `Z_B/Z_A = λ` | P |
| Thm 2.2, open interval with strictly positive exterior | `sharp_open_pos` (core: `exists_pos_weight`) | `[Nonempty O]`, `a ≥ 0`, `T⁺ ≠ ∅` | P |
| Thm 2.2 single-state exterior | `sharp_closed` with `O := Unit` | — | P |
| Realisability by pairwise hidden-variable exteriors | — | — | N |
| Cor 2.2′ | `instance_optimal`, `local_procedure_contains`, `separator_enclosure_tightest` | soundness assumed only for exteriors on one fixed nonempty `O` | P |
| Thm 2.3 | `separation`; exact joint image `achievable_eq`; (i) `Z_unbounded`, `freeEnergy_unbounded` | `separation`: `T⁺ ≠ ∅`, `z > 0`; `achievable_eq`: also `a ≥ 0` and support | P |
| Thm 2.4 kernel form from annulus absorption | `Phi_absorb` | shared factor `k : T → W → T' → ℝ` | P |
| Thm 2.4 mixture identity, κ a probability vector | `rho_push_eq`, `kappa_nonneg`, `kappa_sum_eq_one`, `kappa_eq_zero_of_not_pos`, `push_supp` | **support condition at the inner level**, `0 < Φ'_A(τ')` (plus `φA ≥ 0`, `K ≥ 0` for κ ≥ 0) | P (corrected) |
| Thm 2.4 nesting | `push_mem_Icc`, `hierarchy_nested`; counterexample `hierarchy_needs_support` | `φA ≥ 0`, `K ≥ 0`, support, both supports nonempty | P (corrected) |
| Thm 2.5 | `finite_termination` | `[Unique T]`, `0 < Zsys a e` (no sign hypotheses) | P |
| Thm 2.6 TV inequality | `abs_sum_sub_le_tv_mul` | `Σ_s p = Σ_s p'` (any signs), `m ≤ f ≤ M` on `s` | P |
| Thm 2.6 contraction | `hierarchy_pair_contraction`, `osc_contraction` (δ = any bound on all pairwise TVs, i.e. ≥ the Dobrushin coefficient) | `φA ≥ 0`, support, supports nonempty | P |
| Thm 2.6 products of layer coefficients | — | — | N |
| Executable checker | `checkSeparator : List (ℚ × ℚ) → Option (ℚ × ℚ)`; entries are `(Φ_A τ, Φ_B τ)` | gate: `Φ_A, Φ_B ≥ 0` and `Φ_A = 0 → Φ_B = 0`; `none` if `T⁺ = ∅` | P |
| Checker soundness | `checkSeparator_sound` (marginal Ψ), `checkSeparator_sound_sys` (any interiors with these marginals, any exterior) | `Ψ ≥ 0`, `Z_A > 0` | P |
| Checker exactness | `checkSeparator_attained`, `ratioBounds_spec` | both endpoints are attained local ratios, hence sharp by Thm 2.2 | P |
| §2.3 partially known exteriors | — (an application of T1/T1′) | — | N |
| §2.4 Thm 2.7 soft separator (optional, added later) | — | — | N |

### §3 Binding and selectivity

`Kbind ΨR ΦL ZL = (Σ_τ ΨR τ ΦL τ) / ((Σ_τ ΨR τ) ZL)`; `kloc ΦL ZL τ = ΦL τ / ZL`;
`papo ΨR τ = ΨR τ / Σ ΨR`.

| Spec | Lean name(s) | Exact hypotheses | Status |
|---|---|---|---|
| `K_L = Σ p_apo k_L` | `Kbind_eq_sum`; derived from a microscopic A-loc model in `Kbind_of_sys` | none | P |
| Thm 4.1 enclosure | `affinity_mem_Icc` | `[Nonempty T]`, `ΨR ≥ 0`, `0 < Σ ΨR` | P. *Review note: when Z_L = 0 the statement holds only through Lean's x/0 = 0 convention (K_L = 0, where the physical value is +∞). Use it with Z_L > 0.* |
| Thm 4.1 sharpness | `affinity_sharp` | `[Nonempty T]` | P |
| Thm 4.2 enclosure | `selectivity_mem_Icc` (`selectivity_eq_ratio`: `Z_R` cancels) | `ΨR ≥ 0`, `0 < Σ ΨR`, **`∀ τ, 0 < k₁ τ`** | P (corrected) |
| Thm 4.2 sharpness | `selectivity_sharp` | `∀ τ, 0 < k₁ τ` | P |
| Thm 4.2 ordering (`min k₂/k₁ > 1` ⇒ L2 tighter) | `L2_tighter_of_one_lt_min` | as for the enclosure | P |
| Thm 4.3 with two explicit ensembles | `nonidentifiable` (point masses `pointEnsemble τp`, `pointEnsemble τm`), `nonidentifiable_minmax` | `k₁ > 0` at the two states, `k₂/k₁ < 1` at `τm`, `> 1` at `τp` | P |
| Thm 4.3 with strictly positive ensembles | `nonidentifiable_pos` | `∀ τ, 0 < k₁ τ`, `min < 1 < max` | P |
| Thm 4.4 | — | — | N (not requested) |

## Deviations from the spec

1. **Theorems 1.1 and 1.2: the hypothesis `Adm ≠ ∅` is not needed.** The §1.1 data conditions
   (`WellFormed`) are enough. If `Adm = ∅`, then every `U k = 0` (`adm_nonempty_iff`), so
   `g₋ ≡ g₊ ≡ 0` and both sides of the iff hold. The spec's version is kept as
   `forall_le_Q_iff_of_nonempty`. The data conditions themselves are needed for (⇒).
   Counterexample with `L < 0`: one region `L = −1, U = 1, r = R = 1`, `t = 1/2`. Every
   admissible `Q` equals `1 ≥ 1/2`, but `g₋(1/2) = −1/2 < 0`. The soundness directions (⇐)
   hold for arbitrary data.
2. **Theorem 1.4, upper side: `λ₊ ≤ S₊` needs `Σ L > 0`.** Otherwise `S₊ = +∞`, which Lean
   encodes as `0`. The `max R` part needs no such hypothesis.
3. **Theorem 1.5:** formalised exactly as stated. The hypothesis `S₋ > 0` is necessary.
   Counterexample: regions `(L,U,r) = (0,1,1)` and `(1,1,0)` give `S₋ = λ₋ = 0`, yet
   `r₁(U₁ − L₁) = 1`.
4. **Theorem 1.6** needs `∃ k, U_k > 0`; otherwise `λ₋` does not exist. "Whatever the mass
   bounds are" holds only under this condition.
5. **Theorem 2.1:** stated with support in the form `Φ_A = 0 → Φ_B = 0`, which needs no
   `b ≥ 0`. The spec's literal form, which needs `b ≥ 0`, is `local_enclosure_spec`.
6. **Theorem 2.2:** the hypothesis `T⁺ = T` is unnecessary, for both the closed interval and
   the strictly positive open-interval claim. Only `T⁺ ≠ ∅` is used (plus `a ≥ 0` for
   `Z_A > 0` in the open case). Even when the support condition fails, every value strictly
   between `min_{T⁺} ρ` and `max_{T⁺} ρ` is realised by a strictly positive exterior.
7. **Theorem 2.4 is false without the support condition on the inner separator, which the
   spec does not state for 2.4.** `hierarchy_needs_support`: `T = Bool`,
   `Φ_A = (1, 0)`, `Φ_B = (1, 1)`, `K ≡ 1` onto one new state. Then `ρ = 1` on `T⁺` but
   `ρ' = 2`, so `[min ρ', max ρ'] = {2} ⊄ {1}`. With the condition, the identity and the nesting
   hold, and the condition propagates to the larger separator (`push_supp`). The coordinator's
   random check (`cex_nest.py`) uses strictly positive weights and therefore cannot hit this
   case.
8. **Theorem 4.2 needs `k₁(τ) > 0` for every pocket state**, which the spec leaves implicit.
   `selectivity_needs_positive_k1`: `k₁ = (1, 0)`, `k₂ = (1, 1)` and a uniform receptor give
   `S = 2`, while `k₂/k₁ ≤ 1` wherever `k₁ > 0`. Mathematically `k₂/k₁ = +∞` at the second
   state.
9. **Corollary 2.2′** is proved under a weaker premise: soundness is assumed only for
   exteriors on one fixed nonempty exterior space, for example a single state.
10. **Theorem 2.6:** the TV inequality holds for any weight vectors with equal sums, not only
    probability vectors. δ(κ) enters as any upper bound on all pairwise TV distances.

No spec statement in §1 was found to be false under the §1.1 data conditions. This agrees
with the coordinator's random check.

## Axiom audit

`FoldBounds/AxiomAudit.lean` now lists 143 declarations, 93 of them in the Coupled namespaces.
- Every Coupled report reads `depends on axioms: [propext, Classical.choice, Quot.sound]`.
- No report mentions `Lean.ofReduceBool` or `Lean.trustCompiler`.
- The pre-existing gate prints `Axiom audit: 143 declarations, standard axioms only`.
- `scripts/check_coupled_lean.py` prints `Coupled axiom audit: 93 declarations, standard axioms
  only`.

**Wrapped reports.** Lean wraps long `#print axioms` lines. The existing
`scripts/check_lean_sources.py` reads only the first line of a report, so for a wrapped report
it sees just `[propext,`. That script was not modified. The new script instead parses the
bracketed list across lines and fails on any report header it cannot parse completely. One
Coupled theorem was renamed (`nonidentifiable_minmax`) so that all Coupled reports fit on one
line, which lets the old gate see them completely as well.

The tests use `decide +kernel`: the kernel evaluates the checker on rational input. Plain
`decide` gets stuck, because the elaborator does not unfold the well-founded `Nat.gcd` used
in `ℚ` normalisation. This is not `native_decide` and adds no axioms.

## Reproduce

```bash
export PATH=$HOME/.elan/bin:$PATH
python3 scripts/check_coupled_lean.py                      # source audit (before building)
lake build                                                 # library, tests and audit file
lake env lean FoldBounds/CoupledTests.lean                 # accept/reject tests
lake env lean FoldBounds/AxiomAudit.lean > /tmp/claude-0/ax_lean.log 2>&1
python3 scripts/check_lean_sources.py /tmp/claude-0/ax_lean.log   # existing gate, unchanged
python3 scripts/check_coupled_lean.py /tmp/claude-0/ax_lean.log   # new gate
```

## What these theorems do not establish

- **Model-relative only.** Everything is relative to the abstract finite models above. Nothing
  is claimed about real molecules, experimental affinities or pKa values.
- **Checkers check, they do not search.** They verify proposed endpoints (§1) or compute
  `[min ρ, max ρ]` from given exact rational marginals (§2). No proposer, search or
  threshold-sorting routine is verified.
  - Computing the marginals `Φ_A(τ)`, `Φ_B(τ)` from a factor graph is not verified. Neither is
    the cost of enumerating boundary states, which can be exponential.
- **No efficiency claim.** Efficiency, size-independence of separator widths, and
  completeness of LP dual certificates are not formalised.
- **No link to other code.** The Python implementation, floating-point arithmetic and the 2D HP
  instantiation (§4) are not connected to these theorems.
- **Not formalised.** Theorem 1.3's sorted-threshold form, §2.3, Theorem 2.7 (soft separators),
  Theorem 4.4, realisability by pairwise exteriors, and multi-layer products of contraction
  coefficients.
- **No novelty claim.** The identities used (linear-fractional duality, the Markov
  factorisation, conformational selection) are classical.
