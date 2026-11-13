# Independent-review package: Coupled Ensemble Certificates (PR #96)

**Purpose.** This is a short brief for a reviewer deciding whether PR #96 is sound and worth
merging. It lists the strongest machine-checked results with their exact Lean statements, what
is classical and what may be original, the prior art, the failed validations, and the open
gaps. Supporting detail is in `CEC_THEORY.md` (spec), `LEAN_FORMALIZATION.md`, `PRIOR_ART.md`,
`ADVERSARIAL_REVIEW.md` and `DISCOVERY_REPORT.md`.

**Bottom line.** The PR contributes a small, fully machine-checked library of *sharp* bounds on
ratios of partition functions, Z_B/Z_A, and of *impossibility* results for what local data can
certify. It also contributes verified executable checkers and honest benchmarks. The underlying
mathematics is mostly classical. **No scientific predictive success is claimed:** the one
real-structure model tested loses to a null baseline.

## 1. How to verify (about 5 minutes with the Mathlib cache)

```
lake build && lake env lean FoldBounds/CoupledTests.lean
lake env lean FoldBounds/AxiomAudit.lean > ax.log
python3 scripts/check_lean_sources.py ax.log && python3 scripts/check_coupled_lean.py ax.log
PYTHONPATH=src python3 -m pytest -q
```

Expected output: 143 audited declarations, 93 of them new, all depending only on `propext`,
`Classical.choice` and `Quot.sound`. There is no `sorry` or `native_decide`. The CI jobs
`lean-proof` and `coupled-lean` are green on GitHub.

## 2. Strongest machine-checked results (verbatim statements)

Notation, from the Lean source. Regions are a `Finset ι`. Each region has mass bounds `L ≤ U`
and ratio bounds `r ≤ R`. `WellFormed` means `0 ≤ L ≤ U` and `0 ≤ r ≤ R`. `IsAdm α β` means
`L ≤ α ≤ U`, `rα ≤ β ≤ Rα`, `α, β ≥ 0` and `Σα > 0`. `Q = Σβ/Σα`, and
`gLo t = Σ min((r−t)L, (r−t)U)`. `lamLo` is the infimum of Q over the admissible set, and
`sepLo = Σ rL / Σ U` is the "bound each partition function separately" endpoint.

**(R1) Exact certificate for the coupled lower bound.** A single linear inequality is necessary
and sufficient:

```lean
theorem forall_le_Q_iff (hwf : WellFormed s L U r R) (t : ℝ) :
    (∀ α β, IsAdm s L U r R α β → t ≤ Q s α β) ↔ 0 ≤ gLo s L U r t
```

**(R2) When coupling strictly beats separate bounds.** For `sepLo > 0`:

```lean
theorem lamLo_eq_sepLo_iff (hwf : WellFormed s L U r R) (hS : 0 < sepLo s L U r) :
    lamLo s L U r R = sepLo s L U r ↔ ∀ k ∈ s, r k * (U k - L k) = 0
```

There is also a dominance theorem, `coupled_dominates_lower`, and an LP-dual soundness
theorem, `le_Q_of_lowerDualCert`.

**(R3) The executable checker is sound.** This is proved about the function that actually
runs:

```lean
theorem checkEnclosure_sound {lo hi : ℚ} {bounds : ℚ × ℚ}
    (h : checkEnclosure regions lo hi = some bounds) :
    bounds = (lo, hi) ∧ (∃ α β : Fin regions.length → ℝ, ListAdm regions α β) ∧
      ∀ α β : Fin regions.length → ℝ, ListAdm regions α β →
        Q univ α β ∈ Set.Icc (bounds.1 : ℝ) (bounds.2 : ℝ)
```

**Separator setting.** `Zsys x e = Σ_{(i,τ,o)} x(i,τ)·e(τ,o)`: the interior weight `a`
(system A) or `b` (system B) times a shared exterior `e`. `Phi x τ = Σ_i x(i,τ)`. `rhoMin` and
`rhoMax` are the extremes of `Φ_B/Φ_A` over the support `Tplus (Phi a)`.

**(R4) Local enclosure.**

```lean
theorem local_enclosure {e : T → O → ℝ} (ha : ∀ i τ, 0 ≤ a i τ) (he : ∀ τ o, 0 ≤ e τ o)
    (hsupp : ∀ τ, Phi a τ = 0 → Phi b τ = 0) (hZ : 0 < Zsys a e)
    (hT : (Tplus (Phi a)).Nonempty) :
    Zsys b e / Zsys a e ∈ Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT)
```

**(R5) Instance optimality.** Any interval that is sound for every exterior contains the local
enclosure:

```lean
theorem instance_optimal [Nonempty O] (hT : (Tplus (Phi a)).Nonempty) {l u : ℝ}
    (hsound : ∀ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) → 0 < Zsys a e →
      Zsys b e / Zsys a e ∈ Set.Icc l u) :
    Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT) ⊆ Set.Icc l u
```

**(R6) Absolute vs relative separation.** Over all exteriors, the realisable pairs
(Z_A, Z_B/Z_A) are *exactly* a product set. So local data cannot bound Z_A at all, yet they pin
the ratio sharply:

```lean
theorem achievable_eq [Nonempty O] (ha : ∀ i τ, 0 ≤ a i τ)
    (hsupp : ∀ τ, Phi a τ = 0 → Phi b τ = 0) (hT : (Tplus (Phi a)).Nonempty) :
    {p : ℝ × ℝ | ∃ e : T → O → ℝ, (∀ τ o, 0 ≤ e τ o) ∧ 0 < Zsys a e ∧
        p = (Zsys a e, Zsys b e / Zsys a e)} =
      Set.Ioi 0 ×ˢ Set.Icc (rhoMin (Phi a) (Phi b) hT) (rhoMax (Phi a) (Phi b) hT)
```

**Further results.**
- `hierarchy_nested`: separator enclosures are nested under a shared nonnegative kernel. It
  needs the inner support condition; `hierarchy_needs_support` is the formal counterexample
  without it.
- `checkSeparator_sound_sys`: the executable separator checker is sound for every exterior.
- `selectivity_mem_Icc`: certified selectivity, which needs `kloc Φ1 Z1 τ > 0` for every τ;
  `selectivity_needs_positive_k1` is the counterexample without it.
- `nonidentifiable`: when per-state selectivity ratios straddle 1, two explicit receptor
  ensembles reverse the ordering.

**What these theorems do NOT say.** They are about finite sums of nonnegative weights with an
explicit factorisation. They say nothing about whether a physical system has that
factorisation, whether its weights are accurate, or how expensive the extremes are to compute.

## 3. What is classical, and what may be original

| Item | Status |
|---|---|
| Z_B/Z_A = E_A[e^{−βV}] (FEP); the Markov/DLR factorisation; conformational selection K = Σ p_i K_i | **Classical** (Zwanzig 1954; Georgii's book; Monod–Wyman–Changeux 1965). The repo's own PR #91 already uses the FEP identity, with enumeration. |
| R1–R3: linear-fractional optimisation over a box, with a Dinkelbach-type test | **Classical mathematics** (Charnes & Cooper 1962; Dinkelbach 1967). New only as a *verified certificate kernel* for coupled partition-function bounds. |
| R4: convex-hull bound over Markov-blanket conditionals | **Classical** in substance (Leisink & Kappen 2003, bound propagation; Mooij & Kappen 2008). |
| R5–R6: instance optimality, and the exact achievable set separating absolute from relative | **Possibly original as stated.** Elementary, and not found stated this way in the literature reviewed. Medium novelty risk. |
| Soft-separator enclosure (T2.7, hand proof only) | Possibly original as stated; medium risk; **not formalized** |
| Local certifiability ⇔ weak spatial mixing (T2.8) | **Classical.** It is essentially the WSM definition. |
| Proton binding is an antiferromagnetic 2-spin system, which links certified titration to uniqueness/hardness results | The sign structure is known (anticooperativity). The explicit use for certification may be new; medium/low risk. |
| Machine-checked kernel and executable checkers with soundness proofs; independent fail-closed checkers | **The most defensible contribution.** It is engineering and assurance, not new mathematics. |

Closest applied prior art: OSPREY's K*/BBK*/MARK* and Toulbar2's guaranteed weighted counting
give provable partition-function approximations for protein design. To our reading of the MARK*
abstract (PubMed PMID 31855059), they bound each state separately and then combine. **This is
provisional:** the BBK*/MARK* full texts were not read.

## 4. Computational evidence

All benchmarks use exact rational checkers. Soundness failures: **0** everywhere.

| Setting | Result | Interpretation |
|---|---|---|
| 2D Ising / MRF separators | The width shrinks with radius above T_c and plateaus below (0.0007 vs about 0.48 at radius 10). Sharpness is realised exactly (12/12). | This is what the theory predicts. **Size independence holds by construction**, not as evidence, and the spanning-tree "separate" baseline was a strawman (review M2). |
| 2D HP lattice mutations | 71,859 oracle rows, all sound. The coupled width ≈ the trivial ‖V‖∞ bound for n ≥ 12, and is exactly trivial at n = 30 and 50. | **Negative.** Polymers lack a separator. The method adds a checkable rule, not a capability. |
| Real-structure protonation (PROPKA model pKa, Coulomb ε = 4r, 12 Å cutoff) | 10,000 exact-validated pairs, 0 violations. On GroEL (1,064 sites, treewidth ~31) the median certified width is 1e-4 at 8 Å; 19% of sites do not converge. | A genuine *model-computation* capability beyond exact inference. It is **not** a predictor; see §5. |

## 5. Failed or negative scientific validations (preserved)

1. **Model vs experiment (pre-registered).** On PKAD-2, the certified-midpoint RMSE is 1.59
   against 0.74 for the null model (95% CI of the difference +0.34 to +1.34). The correlation
   of predicted and experimental shifts is 0.02. **The model has no predictive skill.**
2. **MC agreement.** Under the registered criterion, only 55–65% of (site, pH) pairs fall
   inside ±3 SE. An earlier "100%" figure used an unregistered slack; it is now disclosed as a
   post-hoc deviation.
3. **HP lattice:** no gain over the trivial bound at realistic sizes.
4. **MRF "size independence"** was tautological. It was reframed after review.
5. **Spec errors.** Five spec statements were wrong as first written, and were caught by the
   Lean formalization and the adversarial review: Theorems 2.4, 4.2, 4.4, §2.3 (missing support
   or positivity conditions) and T2.8 (a sign error). All are corrected, with counterexamples.
6. Pre-existing and unchanged: the registered weighted-ANM RNase A test (in PR #91) FAILED.

## 6. Remaining gaps

- **Python ↔ Lean refinement.**
  - The production checkers (`coupled/checker.py`, `coupled_mrf/checker.py`,
    `coupled_pka/checker.py`) are independent Python, *not* refined to the Lean checkers.
  - The adversarial review cross-checked an independent T1/separator implementation against
    the Lean `#eval` output on 716 inputs and found 0 discrepancies. That is differential
    testing, not proof.
- **Unformalized:** T2.7, T2.8, §2.3, Theorem 4.4, T1′ completeness, the sorted-threshold form
  of T1.3, pairwise-exterior realisability, and multi-layer contraction products.
- **Physics.** The pKa model is deliberately crude: no desolvation, no flexibility, Coulomb
  only. The 19% non-converging sites are not yet attributed to non-uniqueness versus budget
  limits.
- **Prior art to confirm from full texts:** OSPREY BBK*/MARK* ratio handling,
  Leisink–Kappen sharpness statements, and Onufriev et al. (2001) on irregular titration.
- **Assurance gate.** The two soundness bugs in `scripts/check_lean_sources.py` are fixed in
  the separate PR #97, not here. That PR's documented coverage limits still apply.

## 7. Questions for a reviewer

1. Are R5 and R6, together with the verified checkers, a sufficient contribution for a short
   methods or verification paper? Or is the material best kept as internal infrastructure?
2. Should the certified-titration application continue only with a realistic electrostatics
   model and a head-to-head comparison against an established pKa engine?
3. Is merging the Lean library to `main` desirable on its own, separately from the
   experiments?
